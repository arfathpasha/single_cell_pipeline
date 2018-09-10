'''
Created on Feb 19, 2018

@author: dgrewal
'''
import pypeliner
import shutil
import os
from helpers import makedirs


def produce_fastqc_report(fastq_filename, output_html, output_plots, temp_dir,
                          **kwargs):
    makedirs(temp_dir)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_dir,
        fastq_filename,
        dockerize=kwargs.get('dockerize'),
        image=kwargs.get('image'),
        mounts=kwargs.get('mounts'),
        username=kwargs.get("username"),
        password=kwargs.get('password'),
        server=kwargs.get('server'),
    )


    fastq_basename = os.path.basename(fastq_filename)
    if fastq_basename.endswith(".fastq.gz"):
        fastq_basename = fastq_basename.replace(".fastq.gz", "")
    elif fastq_basename.endswith(".fq.gz"):
        fastq_basename = fastq_basename.replace(".fq.gz", "")
    else:
        raise Exception("Unknown file type")

    output_basename = os.path.join(temp_dir, fastq_basename)

    shutil.move(output_basename + '_fastqc.zip', output_plots)
    shutil.move(output_basename + '_fastqc.html', output_html)


def bwa_mem_paired_end(fastq1, fastq2, output,
                         reference, readgroup,
                         **kwargs):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """
    pypeliner.commandline.execute(
        'bwa', 'mem', '-M', '-R', readgroup,
        reference, fastq1, fastq2,
        '>', output,
        dockerize=kwargs.get('dockerize'),
        image=kwargs.get('image'),
        mounts=kwargs.get('mounts'),
        username=kwargs.get("username"),
        password=kwargs.get('password'),
        server=kwargs.get('server'),
    )

def samtools_sam_to_bam(samfile, bamfile,
                         **kwargs):

    pypeliner.commandline.execute(
        'samtools', 'view', '-bSh', samfile,
        '>', bamfile,
        dockerize=kwargs.get('dockerize'),
        image=kwargs.get('image'),
        mounts=kwargs.get('mounts'),
        username=kwargs.get("username"),
        password=kwargs.get('password'),
        server=kwargs.get('server'),
    )



def bwa_aln_paired_end(fastq1, fastq2, output, tempdir,
                         reference, readgroup,
                         **kwargs):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    read_1_sai = os.path.join(tempdir, 'read_1.sai')
    read_2_sai = os.path.join(tempdir, 'read_2.sai')

    pypeliner.commandline.execute(
        'bwa',
        'aln',
        reference,
        fastq1,
        '>',
        read_1_sai,
        dockerize=kwargs.get('dockerize'),
        image=kwargs.get('image'),
        mounts=kwargs.get('mounts'),
        username=kwargs.get("username"),
        password=kwargs.get('password'),
        server=kwargs.get('server'),
    )

    pypeliner.commandline.execute(
        'bwa',
        'aln',
        reference,
        fastq2,
        '>',
        read_2_sai,
        dockerize=kwargs.get('dockerize'),
        image=kwargs.get('image'),
        mounts=kwargs.get('mounts'),
        username=kwargs.get("username"),
        password=kwargs.get('password'),
        server=kwargs.get('server'),
    )

    pypeliner.commandline.execute(
        'bwa', 'sampe',
        '-r', readgroup,
        reference,
        read_1_sai,
        read_2_sai,
        fastq1,
        fastq2,
        '>',
        output,
        dockerize=kwargs.get('dockerize'),
        image=kwargs.get('image'),
        mounts=kwargs.get('mounts'),
        username=kwargs.get("username"),
        password=kwargs.get('password'),
        server=kwargs.get('server'),
    )


def bam_index(infile, outfile, **kwargs):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile,
        dockerize=kwargs.get('dockerize'),
        image=kwargs.get('image'),
        mounts=kwargs.get('mounts'),
        username=kwargs.get("username"),
        password=kwargs.get('password'),
        server=kwargs.get('server'),
    )


def bam_flagstat(bam, metrics, **kwargs):

    pypeliner.commandline.execute(
        'samtools', 'flagstat',
        bam,
        '>',
        metrics,
        dockerize=kwargs.get('dockerize'),
        image=kwargs.get('image'),
        mounts=kwargs.get('mounts'),
        username=kwargs.get("username"),
        password=kwargs.get('password'),
        server=kwargs.get('server'),
    )


def bam_merge(bams, output, **kwargs):

    cmd = ['samtools', 'merge']
    if kwargs.get('region'):
        cmd.extend(['-R', kwargs.get('region')])

    cmd.append(output)
    cmd.extend(bams)

    kwargs = {
        'dockerize':kwargs.get('dockerize'),
        'image':kwargs.get('image'),
        'mounts':kwargs.get('mounts'),
        'username':kwargs.get("username"),
        'password':kwargs.get('password'),
        'server':kwargs.get('server'),
    }

    pypeliner.commandline.execute(*cmd, **kwargs)


def bam_view(bam, output, region, **kwargs):

    cmd = ['samtools', 'view', '-b', bam, '-o', output, region]
    pypeliner.commandline.execute(*cmd, **kwargs)

