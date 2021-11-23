import argparse
from telomere_utils.extract_reads import extract_telomeric_reads_and_metrics
from telomere_utils.merge import merge
from telomere_utils.overlap import get_overlap
from telomere_utils.split_collated_bam import split_bam

def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()

    extract_reads = subparsers.add_parser("extract_reads")
    extract_reads.set_defaults(which='extract_reads')
    extract_reads.add_argument(
        'input',
        help='input bam file'
    )

    extract_reads.add_argument(
        'outbam',
        help='output directory'
    )
    extract_reads.add_argument(
        'outcsv',
        help='output directory'
    )
    extract_reads.add_argument(
        'sample_id',
        help='output directory'
    )
    extract_reads.add_argument(
        '--perc_threshold',
        type=float,
        default=0.85,
        help='output directory'
    )
    extract_reads.add_argument(
        '--mapping_quality',
        type=int,
        default=30,
        help='output directory'
    )
    extract_reads.add_argument(
        '--telomere_length_threshold',
        type=int,
        default=36,
        help='output directory'
    )

    merge_files = subparsers.add_parser("merge_files")
    merge_files.set_defaults(which='merge_files')
    merge_files.add_argument(
        '--inputs',
        nargs='*',
        help='input bam file'
    )
    merge_files.add_argument(
        '--output',
        help='input bam file'
    )

    get_overlap = subparsers.add_parser("get_overlap")
    get_overlap.set_defaults(which='get_overlap')
    get_overlap.add_argument(
        '--normal_bam',
        help='input bam file'
    )
    get_overlap.add_argument(
        '--normal_data',
        help='input bam file'
    )
    get_overlap.add_argument(
        '--tumour_data',
        help='input bam file'
    )
    get_overlap.add_argument(
        '--output',
        help='input bam file'
    )
    get_overlap.add_argument(
        '--bin_counts',
        help='input bam file'
    )

    split_bam = subparsers.add_parser("split_bam")
    split_bam.set_defaults(which='split_bam')
    split_bam.add_argument(
        '--infile',
        help='input bam file'
    )
    split_bam.add_argument(
        '--outdir',
        help='input bam file'
    )

    args = parser.parse_args()

    args = vars(args)
    return args



def main():
    args = parse_args()

    if args['which'] == 'extract_reads':
        extract_telomeric_reads_and_metrics(
            args['input'],
            args['outbam'],
            args['outcsv'],
            args['sample_id'],
            perc_threshold=args['perc_threshold'],
            mapping_quality=args['mapping_quality'],
            telomere_length_threshold=args['telomere_length_threshold']
        )
    if args['which'] == 'merge_files':
        merge(
            args['inputs'],
            args['output']
        )
    if args['which'] == 'get_overlap':
        get_overlap(
            args['normal_bam'],
            args['normal_data'],
            args['tumour_data'],
            args['output'],
            args['bin_counts'],
        )
    if args['which'] == 'split_bam':
        split_bam(
            args['infile'],
            args['outdir']
        )



