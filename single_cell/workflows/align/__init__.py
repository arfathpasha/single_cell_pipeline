'''
Created on Jul 6, 2017

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
import single_cell
from single_cell.utils import helpers

def create_alignment_workflow(
        fastq_1_filename,
        fastq_2_filename,
        bam_filename,
        bai_filename,
        ref_genome,
        config,
        args,
        instrumentinfo,
        centerinfo,
        sample_info,
        cell_ids):

    out_dir = args['out_dir']

    merge_metrics = os.path.join(out_dir, 'metrics')

    lane_metrics = os.path.join(args['out_dir'], 'metrics_per_lane', '{lane}')

    bam_filename = dict([(cellid, bam_filename[cellid])
                         for cellid in cell_ids])

    bai_filename = dict([(cellid, bai_filename[cellid])
                         for cellid in cell_ids])

    chromosomes = config["chromosomes"]

    ctx = {'mem_retry_increment': 2, 'ncpus': 1}
    docker_ctx = helpers.get_container_ctx(config['containers'], 'single_cell_pipeline')
    ctx.update(docker_ctx)

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id', 'lane'),
        value=fastq_1_filename.keys(),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('instrument', 'cell_id', 'lane', axes_origin=[]),
        value=instrumentinfo)

    workflow.setobj(
        obj=mgd.TempOutputObj('center', 'cell_id', 'lane', axes_origin=[]),
        value=centerinfo)

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo', 'cell_id', axes_origin=[]),
        value=sample_info)


    fastqc_reports = os.path.join(
        lane_metrics,
        "fastqc",
        "{cell_id}_reports.tar.gz")
    flagstat_metrics = os.path.join(lane_metrics, 'flagstat', '{cell_id}.txt')
    workflow.transform(
        name='align_reads',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        axes=('cell_id', 'lane',),
        func="single_cell.workflows.align.tasks.align_pe",
        args=(
            mgd.InputFile(
                'fastq_1', 'cell_id', 'lane', fnames=fastq_1_filename),
            mgd.InputFile(
                'fastq_2', 'cell_id', 'lane', fnames=fastq_2_filename),
            mgd.TempOutputFile(
                'aligned_per_cell_per_lane.sorted.bam', 'cell_id', 'lane'),
            mgd.OutputFile(fastqc_reports, 'cell_id', 'lane'),
            mgd.OutputFile(flagstat_metrics, 'cell_id', 'lane'),
            mgd.TempSpace('alignment_temp', 'cell_id', 'lane'),
            ref_genome,
            mgd.TempInputObj('instrument', 'cell_id', 'lane'),
            mgd.TempInputObj('center', 'cell_id', 'lane'),
            mgd.TempInputObj('sampleinfo', 'cell_id'),
            mgd.InputInstance('cell_id'),
            mgd.InputInstance('lane'),
            args['library_id'],
            config
        )
    )

    workflow.transform(
        name='merge_bams',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        func="single_cell.workflows.align.tasks.merge_bams",
        axes=('cell_id',),
        args=(
            mgd.TempInputFile(
                'aligned_per_cell_per_lane.sorted.bam',
                'cell_id',
                'lane'),
            mgd.TempOutputFile('merged_lanes.bam', 'cell_id'),
            mgd.TempOutputFile('merged_lanes.bam.bai', 'cell_id'),
            config
        )
    )

    if args['realign']:
        workflow.transform(
            name='realignment',
            axes=('chrom',),
            ctx=dict(mem=config['memory']['high'],
                     pool_id=config['pools']['highmem'],
                     **ctx),
            func="single_cell.workflows.align.tasks.realign",
            args=(
                mgd.TempInputFile('merged_lanes.bam', 'cell_id'),
                mgd.TempInputFile('merged_lanes.bam.bai', 'cell_id'),
                mgd.TempOutputFile('realigned.bam', 'chrom', 'cell_id'),
                mgd.TempSpace('realignment_temp', 'chrom', cleanup='before'),
                config,
                mgd.InputInstance('chrom')
            )
        )

        workflow.transform(
            name='merge_realignment',
            ctx=dict(mem=config['memory']['high'],
                     pool_id=config['pools']['highmem'],
                     **ctx),
            axes=('cell_id',),
            func="single_cell.workflows.align.tasks.merge_realignment",
            args=(
                mgd.TempInputFile('realigned.bam', 'chrom', 'cell_id'),
                mgd.TempOutputFile('merged_realign.bam', 'cell_id'),
                config,
                mgd.InputInstance('cell_id')
            )
        )

    final_bam = mgd.TempInputFile('merged_lanes.bam', 'cell_id')
    if args["realign"]:
        final_bam = mgd.TempInputFile('merged_realign.bam', 'cell_id')

    markdups_metrics = os.path.join(
        merge_metrics,
        'markdups_metrics',
        '{cell_id}.markdups_metrics.txt')
    flagstat_metrics = os.path.join(
        merge_metrics,
        'flagstat_metrics',
        '{cell_id}.flagstat_metrics.txt')
    workflow.transform(
        name='postprocess_bam',
        ctx=dict(mem=config['memory']['med'],
                 pool_id=config['pools']['standard'],
                 **ctx),
        axes=('cell_id',),
        func="single_cell.workflows.align.tasks.postprocess_bam",
        args=(
            final_bam,
            mgd.OutputFile('sorted_markdups', 'cell_id', fnames=bam_filename),
            mgd.OutputFile(
                'sorted_markdups_index', 'cell_id', fnames=bai_filename),
            mgd.TempSpace('tempdir', 'cell_id'),
            config,
            mgd.OutputFile(markdups_metrics, 'cell_id'),
            mgd.OutputFile(flagstat_metrics, 'cell_id'),
        ),
    )

    return workflow
