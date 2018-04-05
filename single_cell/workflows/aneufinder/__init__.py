import os
import pypeliner
import pypeliner.managed as mgd

import tasks

from single_cell.utils import csvutils

def create_aneufinder_workflow(bam_file,
                               cell_ids,
                               config,
                               aneufinder_output,
                               aneufinder_segs_filename,
                               aneufinder_reads_filename,
                               library_id):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=cell_ids,
    )

    workflow.transform(
        name='run_aneufinder_on_individual_cells',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.run_aneufinder,
        axes=('cell_id',),
        args=(
            mgd.InputFile('bam_file', 'cell_id', fnames=bam_file),
            mgd.TempSpace('working_dir', 'cell_id', fnames=bam_file),
            mgd.InputInstance('cell_id'),
            aneufinder_output,
            mgd.TempOutputFile('segments.csv', 'cell_id'),
            mgd.TempOutputFile('reads.csv', 'cell_id'),
            mgd.TempOutputFile('dnacopy.pdf', 'cell_id'),
        ),
    )

    workflow.transform(
        name='merge_segments',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=csvutils.concatenate_csv_lowmem,
        args=(
            mgd.TempInputFile('segments.csv', 'cell_id'),
            mgd.OutputFile(aneufinder_segs_filename)
        ),
    )

    workflow.transform(
        name='merge_reads',
        ctx={'mem': config["memory"]['low'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=csvutils.concatenate_csv_lowmem,
        args=(
            mgd.TempInputFile('reads.csv', 'cell_id'),
            mgd.OutputFile(aneufinder_reads_filename)
        )
    )

    dnacopy_pdf_output = os.path.join(aneufinder_output, 'plots', '{}_reads.pdf'.format(library_id))
    workflow.transform(
        name='merge_aneufinder_pdfs',
        ctx={'mem': config["memory"]['med'], 'pool_id': config['pools']['standard'], 'ncpus':1},
        func=tasks.merge_pdf,
        args=(
            [mgd.TempInputFile('dnacopy.pdf', 'cell_id')],
            [mgd.OutputFile(dnacopy_pdf_output)],
        )
    )

    return workflow
