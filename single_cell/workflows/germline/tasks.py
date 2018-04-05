import pypeliner
import vcf


def run_samtools_variant_calling(
        bam_file,
        bai_file,
        ref_genome_fasta_file,
        out_file,
        max_depth=int(1e7),
        min_bqual=0,
        min_depth=0,
        min_mqual=0,
        region=None):

    mpileup_cmd = [
        'samtools',
        'mpileup',
        '-ugf', ref_genome_fasta_file,
        '-Q', min_bqual,
        '-q', min_mqual,
        bam_file
    ]

    if region is not None:
        mpileup_cmd.extend(['-r', region])

    bcf_cmd = [
        'bcftools',
        'call',
        '-vmO', 'z',
        '-o', out_file,
    ]

    cmd = []

    cmd.extend(mpileup_cmd)
    cmd.append('|')
    cmd.extend(bcf_cmd)

    pypeliner.commandline.execute(*cmd)

    index_bcf(out_file)


nucleotides = ('A', 'C', 'G', 'T')


def parse_samtools_vcf(
        vcf_file,
        out_file):

    vcf_reader = vcf.Reader(filename=vcf_file)

    data = []

    for record in vcf_reader:
        ref_base = record.REF

        # Skip record with reference base == N
        if ref_base not in nucleotides:
            continue

        for alt_base in record.ALT:
            alt_base = str(alt_base)

            if (len(ref_base) != 1) or (len(alt_base) != 1):
                continue

            # Skip record with alt base == N
            if alt_base not in nucleotides:
                continue

            # Format output record
            out_row = {
                'chromosome': record.CHROM,
                'start': record.POS,
                'stop': record.POS,
                'ref': ref_base,
                'alt': alt_base,
            }

            data.append(out_row)

    data = pd.DataFrame(data)
    data['case_id'] = 'NA'

    data.to_csv(vcf_file, sep=',')