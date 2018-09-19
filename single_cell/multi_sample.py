import pypeliner
import pypeliner.managed as mgd

import infer_haps


def get_default_ctx(config):
    singlecellimage = config['docker']['images']['single_cell_pipeline']

    ctx = {
        'mem_retry_increment': 2,
        'ncpus': 1,
        'image': singlecellimage['image'],
        'dockerize': config['docker']['dockerize'],
        'mounts': config['docker']['mounts'],
        'username': singlecellimage['username'],
        'password': singlecellimage['password'],
        'server': singlecellimage['server'],
    }

    return ctx


def multi_sample_workflow(
    normal_wgs_bam,
    tumour_cell_bams,
    config,
):
    """ Multiple sample pseudobulk workflow. """

    normal_region_bam_template
    tumour_region_bam_template
    museq_vcf_template
    strelka_snv_template
    strelka_indel_template
    snv_annotations_template
    normal_seqdata_file
    tumour_cell_seqdata_template
    haplotypes_file
    allele_counts_filename

    workflow = pypeliner.workflow.Workflow(
        default_ctx=get_default_ctx(config))

    workflow.set_filenames('normal_regions.bam', 'region', template=normal_region_bam_template)
    workflow.set_filenames('tumour_cells.bam', 'sample_id', 'cell_id', fnames=tumour_cell_bams)
    workflow.set_filenames('tumour_regions.bam', 'sample_id', 'region', template=tumour_region_bam_template)

    workflow.set_filenames('museq.vcf', 'sample_id', template=museq_vcf_template)
    workflow.set_filenames('strelka_snv.vcf', 'sample_id', template=strelka_snv_template)
    workflow.set_filenames('strelka_indel.vcf', 'sample_id', template=strelka_indel_template)
    workflow.set_filenames('snv_annotations.h5', 'sample_id', template=snv_annotations_template)
    workflow.set_filenames('cells.h5', 'sample_id', 'cell_id', template=tumour_cell_seqdata_template)

    workflow.subworkflow(
        name='split_normal',
        func=split_bams.create_split_workflow,
        args=(
            mgd.InputFile(normal_wgs_bam, extensions=['.bai']),
            mgd.OutputFile('normal_regions.bam', 'region', extensions=['.bai'], axes_origin=[])
            'region',
            config,
        ),
    )

    workflow.subworkflow(
        name='split_merge_tumour',
        func=merge_bams.create_merge_bams_workflow,
        args=(
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.OutputFile('tumour_regions.bam', 'sample_id', 'region', axes_origin=[], extensions=['.bai']),
            cellids,
            config,
            mgd.TempInputObj('region'),
        )
    )

    workflow.subworkflow(
        name='variant_calling',
        func=variant_calling_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.InputFile('tumour_regions.bam', 'sample_id', 'region', axes_origin=[], extensions=['.bai']),
            mgd.InputFile('normal_regions.bam', 'region', extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('museq.vcf', 'sample_id'),
            mgd.OutputFile('strelka_snv.vcf', 'sample_id'),
            mgd.OutputFile('strelka_indel.vcf', 'sample_id'),
            mgd.OutputFile('snv_annotations.h5', 'sample_id'),
            config,
        ),
    )

    workflow.transform(
        name='merge_museq_snvs',
        func='single_cell.utils.vcfutils.merge_vcfs',
        args=(
            mgd.OutputFile('museq.vcf', 'sample_id'),
            mgd.TempOutputFile('museq.vcf'),
        ),
    )

    workflow.transform(
        name='merge_strelka_snvs',
        func='single_cell.utils.vcfutils.merge_vcfs',
        args=(
            mgd.OutputFile('strelka_snv.vcf', 'sample_id'),
            mgd.TempOutputFile('strelka_snv.vcf'),
        ),
    )

    workflow.subworkflow(
        name='variant_counting',
        func=variant_counting_workflow,
        axes=('sample_id',),
        args=(
            [
                mgd.TempInputFile('museq.vcf'),
                mgd.TempInputFile('strelka_snv.vcf'),
            ],
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.InputFile('tumour_regions.bam', 'sample_id', 'region', axes_origin=[], extensions=['.bai']),
            mgd.InputFile('normal_regions.bam', 'region', extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('snv_counts.h5'),
            config,
        ),
    )

    workflow.subworkflow(
        name='infer_haps_from_bulk_normal',
        func=infer_haps.infer_haps_from_bulk_normal,
        args=(
            mgd.InputFile(normal_wgs_bam, extensions=['.bai']),        
            mgd.OutputFile(normal_seqdata_file),
            mgd.OutputFile(haplotypes_file),
            config,
        ),
    )

    workflow.subworkflow(
        name='extract_allele_readcounts',
        func=infer_haps.extract_allele_readcounts,
        axes=('sample_id',),
        args=(
            mgd.InputFile(haplotypes_filename),
            mgd.InputFile('tumour_cells.bam', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.OutputFile('cells.h5', 'sample_id', 'cell_id', extensions=['.bai']),
            mgd.InputFile(allele_counts_filename),
            config,
        ),
    )

