'''
Created on Jun 6, 2018

@author: dgrewal
'''
import single_cell
import os
import yaml
import collections

import pipeline_config


class folded_unicode(unicode):
    pass


class literal_unicode(unicode):
    pass


def folded_unicode_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='>')


def literal_unicode_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(folded_unicode, folded_unicode_representer)
yaml.add_representer(literal_unicode, literal_unicode_representer)


def override_config(config, override):
    def update(d, u):
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        return d

    if not override:
        return config

    cfg = update(config, override)

    return cfg


def get_version():
    version = single_cell.__version__
    # strip setuptools metadata
    version = version.split("+")[0]
    return version


def get_batch_params(override=None):

    pools = {
        "standard": {
            "tasks_per_node": 4,
            "numcores": 4,
            "primary": True
        },
        "highmem": {
            "tasks_per_node": 2,
            "numcores": 4,
            "primary": False
        },
        "multicore": {
            "tasks_per_node": 1,
            "numcores": 8,
            "primary": False
        }
    }

    data = {
        "storage_container_name": "tasks-container",
        "no_delete_pool": True,
        "no_delete_job": False,
        "pools": pools,
        "version": get_version(),
        "reference": "grch37"
    }

    data = override_config(data, override)

    data["version"] = data["version"].replace('.', '_')

    return data


def generate_autoscale_formula(tasks_per_node):

    formula = (
        "tasksPerNode = {};\n"
        "numAddMax = 20;\n"
        "numDelMax = 20;\n"
        "startingNumberOfVMs = 0;\n"
        "minNumberofVMs = 0;\n"
        "maxNumberofVMs = 1000;\n"
        "pendingTaskSamplePercent = $PendingTasks.GetSamplePercent(180 * TimeInterval_Second);\n"
        "pendingTaskSamples = pendingTaskSamplePercent < 70 ? startingNumberOfVMs : avg($PendingTasks.GetSample(180 * TimeInterval_Second));\n"
        "cores = $TargetLowPriorityNodes * tasksPerNode;\n"
        "$extraVMs = (pendingTaskSamples - cores) / tasksPerNode;\n"
        "$extraVMs = $extraVMs + (tasksPerNode-1)/tasksPerNode;\n"
        "$extraVMs = min(numAddMax, $extraVMs);\n"
        "$extraVMs = max(-numDelMax, $extraVMs);\n"
        "targetVMs = ($TargetLowPriorityNodes + $extraVMs);\n"
        "$TargetLowPriorityNodes = max(minNumberofVMs,min(maxNumberofVMs, targetVMs));\n"
    )

    formula = formula.format(tasks_per_node)

    formula = literal_unicode(formula)

    return formula


def create_vm_commands():

    # GETSIZE AND MOUNT:
    # the disks are not guaranteed to appear under the same device name when starting a new
    # vm from a managed image. So we cannot be sure where the devices should be mounted,
    # this code uses the disk size to guess where to mount
    # If disk in under 40ish GB then mount it to refdata (our default reference disk is 32 GB)
    # If disk is over 900ish GB then mount it to datadrive (our default scratch space is 1TB)
    # DOCKER GROUP:
    # We do not have sudo available in tasks, docker requires sudo. Docker also comes with a
    # 'docker' user group. This user group has escalated permissions, any user in this group
    # can run docker without sudo (but still as an escalated user). Add the current user (default _azbatch)
    # to the docker user group. This cannot be done during the image creation, image capture removes all user
    # information.
    commands = (
        "if [ `sudo blockdev --getsize64 /dev/sdc` -le 40000000000 ]; then sudo mount /dev/sdc /refdata; else sudo mount /dev/sdd /refdata; fi\n"
        "if [ `sudo blockdev --getsize64 /dev/sdd` -gt 900000000000 ]; then sudo mount /dev/sdd /datadrive; else sudo mount /dev/sdc /datadrive; fi\n"
        "sudo chmod -R 777 /datadrive /refdata\n"
        "sudo gpasswd -a $USER docker"
    )

    commands = literal_unicode(commands)

    return commands


def get_vm_size_azure(numcores):
    if numcores <= 2:
        return "STANDARD_E2_V3"
    elif numcores <= 4:
        return "STANDARD_E4_V3"
    elif numcores <= 8:
        return "STANDARD_E8_V3"
    else:
        # max 16 cores
        return "STANDARD_E16_V3"
#     if numcores <= 2:
#         return "STANDARD_DS11_V2"
#     elif numcores <= 4:
#         return "STANDARD_DS12_V2"
#     elif numcores <= 8:
#         return "STANDARD_DS13_V2"
#     else:
#         # max 16 cores
#         return "STANDARD_DS14_V2"


def get_vm_image_id(version):
    subscription = os.environ.get("SUBSCRIPTION_ID", "id-missing")
    resource_group = os.environ.get("RESOURCE_GROUP", "id-missing")
    return "/subscriptions/{}/resourceGroups/{}/providers/Microsoft.Compute/images/singlecellpipeline_{}".format(
        subscription, resource_group, version)


def get_pool_def(
        tasks_per_node, reference, pool_type, version, numcores, primary=False):

    autoscale_formula = generate_autoscale_formula(tasks_per_node)

    vm_commands = create_vm_commands()

    poolname = "singlecell{}{}_{}".format(reference, pool_type, version)

    pooldata = {
        "pool_vm_size": get_vm_size_azure(numcores),
        "primary": primary,
        'node_resource_id': get_vm_image_id(version),
        'node_os_publisher': 'Canonical',
        'node_os_offer': 'UbuntuServer',
        'node_os_sku': 'batch.node.ubuntu 16.04',
        'data_disk_sizes': None,
        'max_tasks_per_node': tasks_per_node,
        'auto_scale_formula': autoscale_formula,
        'create_vm_commands': vm_commands,
        'start_resources': None
    }

    return {poolname: pooldata}


def get_compute_start_commands():
    # Azure batch starts a process on each node during startup
    # this process runs the entire time VM is up (speculation?). running sudo gasswd
    # adds our user to the docker group in node startup. but all user/group
    # permission updates require restarting the bash process (normally by
    # running newgrp or by logout-login). Since we cannot do that with batch
    # and since the start task and compute task run under same process, we need
    # to explicitly specify the user group when running docker commands. sg in
    # linux can be used to set group when executing commands
    dockeruser = os.environ["CLIENT_ID"]
    dockerpass = os.environ["SECRET_KEY"]

    commands = (
        'clean_up () {',
        '  echo "clean_up task executed"',
        '  find $AZ_BATCH_TASK_WORKING_DIR/ -xtype l -delete',
        '  exit 0',
        '}',
        'trap clean_up EXIT',
        'export PATH=/usr/local/miniconda2/bin/:$PATH',
        'mkdir -p /datadrive/$AZ_BATCH_TASK_WORKING_DIR/',
        'cd /datadrive/$AZ_BATCH_TASK_WORKING_DIR/',
        'sg docker -c "docker login singlecellcontainers.azurecr.io -u {} -p {}"'.format(
            dockeruser,
            dockerpass),
        'sg docker -c "docker pull singlecellcontainers.azurecr.io/scp/single_cell_pipeline"',
    )

    commands = '\n'.join(commands)

    commands = literal_unicode(commands)

    return {"compute_start_commands": commands}


def get_compute_run_commands():

    run_command  = ['pypeliner_delegate',
                   '$AZ_BATCH_TASK_WORKING_DIR/{input_filename}',
                   '$AZ_BATCH_TASK_WORKING_DIR/{output_filename}',
                   ]
    run_command = ' '.join(run_command)
    run_command = literal_unicode(run_command)

    docker_mounts = pipeline_config.get_docker_params(
        'azure')['docker']['mounts']

    mount_string = ['-v {}:{}'.format(mount, mount) for mount in docker_mounts]
    mount_string += ['-v /mnt:/mnt']
    mount_string = ' '.join(mount_string)

    dockerize_run_command = ['docker run -w $PWD',
                             mount_string,
                             '-v /var/run/docker.sock:/var/run/docker.sock',
                             '-v /usr/bin/docker:/usr/bin/docker',
                             'singlecellcontainers.azurecr.io/scp/single_cell_pipeline',
                             '{command}'
                            ]
    dockerize_run_command = ' '.join(dockerize_run_command)

    # wrap it up as docker group command
    dockerize_run_command = 'sg docker -c "{}"'.format(dockerize_run_command)

    # use block format
    dockerize_run_command = literal_unicode(dockerize_run_command)

    return {"compute_run_command": run_command,
            "dockerize_run_command": dockerize_run_command}


def get_compute_finish_commands():
    commands = (
        "find /datadrive/$AZ_BATCH_TASK_WORKING_DIR/ -xtype l -delete\n"
        "find /datadrive/$AZ_BATCH_TASK_WORKING_DIR/ -type f -delete\n"
    )

    commands = literal_unicode(commands)

    return {"compute_finish_commands": commands}


def get_all_pools(pool_config, reference, version):

    pooldefs = {}

    for pooltype, poolinfo in pool_config.iteritems():
        tasks_per_node = poolinfo["tasks_per_node"]
        numcores = poolinfo["numcores"]
        primary = poolinfo["primary"]

        pool_def = get_pool_def(
            tasks_per_node, reference, pooltype,
            version, numcores, primary
        )

        pooldefs.update(pool_def)

    return {"pools": pooldefs}


def write_config(params, filepath):
    with open(filepath, 'w') as outputfile:
        yaml.dump(params, outputfile, default_flow_style=False)


def get_batch_config(defaults, override=None):
    config = {}

    config.update(
        get_all_pools(
            defaults['pools'],
            defaults['reference'],
            defaults['version']))

    config.update(
        {"storage_container_name": defaults["storage_container_name"]}
    )

    config.update(get_compute_start_commands())
    config.update(get_compute_run_commands())
    config.update(get_compute_finish_commands())

    config.update({"no_delete_pool": defaults["no_delete_pool"]})
    config.update({"no_delete_job": defaults["no_delete_job"]})

    config = override_config(config, override)

    return config
