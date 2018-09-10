'''
Created on Nov 21, 2017

@author: dgrewal
'''
import os
import gzip
import pypeliner

from single_cell.utils import helpers
from single_cell.utils import bamutils


def split_bam_worker(bam, output_bam, region, config):

    region = '{}:{}-{}'.format(*region.split('-'))

    bamutils.bam_view(
        bam, output_bam, region,
        dockerize=config['dockerize'],
        mounts=config['mounts'],
        image=config['images']['samtools']['image'],
        username=config['images']['samtools']['username'],
        password=config['images']['samtools']['password'],
        server=config['images']['samtools']['server'],
    )


def index_bam_worker(bam, bai, config):

    bamutils.bam_index(
        bam, bai,
        dockerize=config['dockerize'],
        mounts=config['mounts'],
        image=config['images']['samtools']['image'],
        username=config['images']['samtools']['username'],
        password=config['images']['samtools']['password'],
        server=config['images']['samtools']['server'],
    )


def split_bam_file_one_job(bam, bai, outbam, outbai, regions, docker_config, ncores=None):

    args = [(bam, outbam(region), region, docker_config) for region in regions]

    helpers.run_in_parallel(split_bam_worker, args, ncores=ncores)

    args = [(outbam(region), outbai(region), docker_config) for region in regions]

    helpers.run_in_parallel(index_bam_worker, args, ncores=ncores)


def split_bam_file(bam, bai, outbam, outbai, interval, docker_config):

    bamutils.bam_view(
        bam, outbam, interval,
        dockerize=docker_config['dockerize'],
        mounts=docker_config['mounts'],
        image=docker_config['images']['samtools']['image'],
        username=docker_config['images']['samtools']['username'],
        password=docker_config['images']['samtools']['password'],
        server=docker_config['images']['samtools']['server'],
    )

    bamutils.bam_index(
        outbam, outbai,
        dockerize=docker_config['dockerize'],
        mounts=docker_config['mounts'],
        image=docker_config['images']['samtools']['image'],
        username=docker_config['images']['samtools']['username'],
        password=docker_config['images']['samtools']['password'],
        server=docker_config['images']['samtools']['server'],
    )


def split_bam_file_by_reads(bam, bai, outbams, outbais, tempspace, intervals, docker_config):
    # sort bam by reads and convert to sam

    kwargs = {
        'dockerize':kwargs.get('dockerize'),
        'image':kwargs.get('image'),
        'mounts':kwargs.get('mounts'),
        'username':kwargs.get("username"),
        'password':kwargs.get('password'),
        'server':kwargs.get('server'),
    }

    helpers.makedirs(tempspace)

    headerfile = os.path.join(tempspace, "bam_header.sam")

    cmd = ['samtools', 'view', '-H', bam, '-o', headerfile]
    pypeliner.commandline.execute(*cmd, **kwargs)

    collate_prefix = os.path.join(
        tempspace, os.path.basename(bam) + "_collate_temp"
    )
    collated_bam = os.path.join(tempspace, "bam_file_collated_sam_format.sam")

    cmd = [
        'samtools', 'collate', '-u', '-O', bam, collate_prefix, '|',
        'samtools', 'view', '-', '-o', collated_bam
    ]

    pypeliner.commandline.execute(*cmd, **kwargs)

    tempoutputs = [
        os.path.join(tempspace, os.path.basename(outbams[interval]) + ".split.temp")
        for interval in intervals
    ]

    split(collated_bam, tempoutputs, headerfile=headerfile)

    for inputsam, interval in zip(tempoutputs, intervals):
        outputbam = outbams[interval]

        cmd = ['samtools', 'view', '-Sb', inputsam, '-o', outputbam]

        pypeliner.commandline.execute(*cmd, **kwargs)


def get_file_handle(filename, mode="r"):
    if os.path.splitext(filename)[-1] == ".gz":
        return gzip.open(filename, mode + "b")
    else:
        return open(filename, mode)


def get_file_linecount(infile):
    def blocks(files, size=65536):
        while True:
            b = files.read(size)
            if not b:
                break
            yield b

    with get_file_handle(infile, mode="r") as inputdata:
        numlines = sum(bl.count("\n") for bl in blocks(inputdata))

    return numlines


def split(infile, outfiles, headerfile=None):

    if headerfile:
        with open(headerfile) as headerdata:
            header = headerdata.readlines()

    with get_file_handle(infile, mode="r") as inputdata:

        lines_per_file = get_file_linecount(infile) / len(outfiles)

        file_number = 0
        line_number = 0

        outputfilehandle = get_file_handle(outfiles[file_number], mode='w')
        if headerfile:
            outputfilehandle.writelines(header)

        while True:
            line = inputdata.readline()
            line_number += 1
            if not line:
                break

            outputfilehandle.write(line)

            if line_number > ((file_number + 1) * lines_per_file):
                # if next line has the same readname as last line in file
                # then write to the same file else write to next
                line2 = inputdata.readline()
                line_number += 1
                readname2 = line2.split()[0]
                readname = line.split()[0]

                if readname2 == readname:
                    outputfilehandle.write(line2)

                outputfilehandle.close()
                outputfilehandle = open(outfiles[file_number], 'w')
                if headerfile:
                    outputfilehandle.writelines(header)

                if not readname2 == readname:
                    outputfilehandle.write(line2)

                file_number += 1

        outputfilehandle.close()
