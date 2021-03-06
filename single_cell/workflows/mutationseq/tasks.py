'''
Created on Jul 24, 2017

@author: dgrewal
'''
import os
import pypeliner

from scripts import ParseMuseq

from single_cell.utils import helpers
from single_cell.utils import vcfutils


def run_museq(tumour, normal, out, log, region, config, docker_kwargs={}):
    '''
    Run museq script for each chromosome

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to temporary output VCF file for the chromosome
    :param log: path to the log file
    :param config: path to the config YAML file
    :param chrom: chromosome number
    '''

    reference = config['ref_genome']

    region = '{}:{}-{}'.format(*region.split('-'))

    cmd = ['museq', 'normal:' + normal, 'tumour:' + tumour,
           'reference:' + reference, '--out', out,
           '--log', log, '--interval', region]

    pypeliner.commandline.execute(*cmd, **docker_kwargs)


def concatenate_vcfs(inputs, output):
    vcfutils.concatenate_vcf(inputs, output)


def parse_museq(infile, output):
    parser = ParseMuseq(infile=infile, tid='NA', nid='NA', output=output,
                        keep_dbsnp=True,keep_1000gen=True,
                        remove_duplicates=True)

    parser.main()
