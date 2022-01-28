#!/usr/bin/env python

import argparse
import pandas as pd
import pysam


def _find_telomere_end_in_seq_forwards(query_sequence, kmers):
    end_pos = None
    i = 0

    while i < (len(query_sequence) - 5):
        if any(tr == query_sequence[i:(i + 6)] for tr in kmers):
            end_pos = i + 6
            i += 6
        else:
            i += 1

    return 0, end_pos


def _find_telomere_end_in_seq_backwards(query_sequence, kmers):
    end_pos = None
    i = len(query_sequence)

    while i > 5:
        if any(tr == query_sequence[(i - 6):i] for tr in kmers):
            end_pos = i - 6
            i -= 6
        else:
            i -= 1

    return end_pos, len(query_sequence)


def _keep_telomere(start_pos, end_pos, perc_threshold=0.85, length_threshold=36):
    # no telomere found
    if end_pos is None or start_pos is None:
        return False

    telo_len = end_pos - start_pos

    # telomere not long enough
    if telo_len < length_threshold:
        return False

    percentage = telo_len / end_pos
    # most of the telomeric section is not telomeric enough
    if percentage < perc_threshold:
        return False

    return True


def _find_all_telomeres(read, perc_threshold=0.85, mapping_quality=30, length_threshold=36):
    if read.mapping_quality > mapping_quality:
        return []

    kmers = ['TTAGGG', 'TCAGGG', 'TGAGGG', 'TTGGGG']
    kmers_revcomp = ['CCCTAA', 'CCCTGA', 'CCCTCA', 'CCCCAA']

    query_sequence = read.query_sequence

    positions = []

    telo_start, telo_end = _find_telomere_end_in_seq_forwards(query_sequence, kmers)
    if _keep_telomere(telo_start, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append(
            {'start': telo_start, 'end': telo_end, 'kmers_revcomp': False, 'read_rev': False, 'forwards': True}
        )

    telo_start, telo_end = _find_telomere_end_in_seq_forwards(query_sequence, kmers_revcomp)
    if _keep_telomere(telo_start, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append(
            {'start': telo_start, 'end': telo_end, 'kmers_revcomp': True, 'read_rev': False, 'forwards': True}
        )

    telo_start, telo_end = _find_telomere_end_in_seq_forwards(query_sequence[::-1], kmers)
    if _keep_telomere(telo_start, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append(
            {'start': telo_start, 'end': telo_end, 'kmers_revcomp': False, 'read_rev': True, 'forwards': True}
        )

    telo_start, telo_end = _find_telomere_end_in_seq_forwards(query_sequence[::-1], kmers_revcomp)
    if _keep_telomere(telo_start, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append(
            {'start': telo_start, 'end': telo_end, 'kmers_revcomp': True, 'read_rev': True, 'forwards': True}
        )

    ## backwards
    telo_start, telo_end = _find_telomere_end_in_seq_backwards(query_sequence, kmers)
    if _keep_telomere(telo_start, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append(
            {'start': telo_start, 'end': telo_end, 'kmers_revcomp': False, 'read_rev': False, 'forwards': False}
        )

    telo_start, telo_end = _find_telomere_end_in_seq_backwards(query_sequence, kmers_revcomp)
    if _keep_telomere(telo_start, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append(
            {'start': telo_start, 'end': telo_end, 'kmers_revcomp': True, 'read_rev': False, 'forwards': False}
        )

    telo_start, telo_end = _find_telomere_end_in_seq_backwards(query_sequence[::-1], kmers)
    if _keep_telomere(telo_start, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append(
            {'start': telo_start, 'end': telo_end, 'kmers_revcomp': False, 'read_rev': True, 'forwards': False}
        )

    telo_start, telo_end = _find_telomere_end_in_seq_backwards(query_sequence[::-1], kmers_revcomp)
    if _keep_telomere(telo_start, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append(
            {'start': telo_start, 'end': telo_end, 'kmers_revcomp': True, 'read_rev': True, 'forwards': False}
        )

    return positions


def _get_read_pairs(bamfile):
    read_data = {}

    for al in bamfile.fetch(until_eof=True):

        # must be primary read alignment, not secondary or supplementary
        if al.is_secondary or al.flag & 2048 == 2048:
            continue

        # skip unpaired reads
        if not al.is_paired:
            continue

        # add read name to dictionary if not already there
        key = al.qname
        if key not in read_data:
            read_data.setdefault(key, al)
        # print matched read pairs
        else:
            yield al, read_data[key]


def _get_csv_data(read, telomere, sample_id, bamfile):
    if bamfile.get_reference_name(read.reference_id) is None:
        return

    if bamfile.get_reference_name(read.reference_id) is None:
        return

    data = {
        'read_id': read.query_name,
        'sample_id': sample_id,
        'strand': '-' if read.is_reverse else '+',
        'chromosome': bamfile.get_reference_name(read.reference_id),
        'start': read.reference_start,
        'end': read.reference_end,
        'telomere_start': telomere['start'],
        'telomere_end': telomere['end'],
        'reverse_complemented_kmers': telomere['kmers_revcomp'],
        'reversed_read': telomere['read_rev'],
        'match_beginning_of_read': telomere['forwards'],
        'readend': '1' if read.is_read1 else '2'
    }

    return data


def extract_telomeric_reads_and_metrics(
        infile, outbam, outcsv, sample_id, perc_threshold=0.85, mapping_quality=30,
        length_threshold=36
):
    bamfile = pysam.AlignmentFile(infile, "rb")
    outfile = pysam.AlignmentFile(outbam, "wb", template=bamfile)

    csvdata = []

    for r1, r2 in _get_read_pairs(bamfile):
        telo_r1 = _find_all_telomeres(
            r1, perc_threshold=perc_threshold, mapping_quality=mapping_quality,
            length_threshold=length_threshold
        )
        telo_r2 = _find_all_telomeres(
            r2, perc_threshold=perc_threshold, mapping_quality=mapping_quality,
            length_threshold=length_threshold
        )

        if len(telo_r1) == 0 and len(telo_r2) == 0:
            continue

        outfile.write(r1)
        outfile.write(r2)

        for telomere in telo_r1:
            assert None not in telomere
            read_pair_data = _get_csv_data(r1, telomere, sample_id, bamfile)
            if read_pair_data is not None:
                csvdata.append(read_pair_data)

        for telomere in telo_r2:
            assert None not in telomere
            read_pair_data = _get_csv_data(r2, telomere, sample_id, bamfile)
            if read_pair_data is not None:
                csvdata.append(read_pair_data)

    df = pd.DataFrame(csvdata)

    if df.empty:
        columns = [
            'read_id', 'sample_id', 'strand', 'chromosome',
            'start', 'end', 'telomere_start', 'telomere_end',
            'reverse_complement', 'readend'
        ]
        df = pd.DataFrame(columns=columns)

    df.to_csv(outcsv, index=False, na_rep='NA')


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input',
        help='input bam file'
    )

    parser.add_argument(
        'outbam',
        help='output directory'
    )
    parser.add_argument(
        'outcsv',
        help='output directory'
    )
    parser.add_argument(
        'sample_id',
        help='output directory'
    )
    parser.add_argument(
        '--perc_threshold',
        type=float,
        default=0.85,
        help='output directory'
    )
    parser.add_argument(
        '--mapping_quality',
        type=int,
        default=30,
        help='output directory'
    )
    parser.add_argument(
        '--telomere_length_threshold',
        type=int,
        default=36,
        help='output directory'
    )

    args = parser.parse_args()

    args = vars(args)
    return args


def main():
    args = parse_args()
    extract_telomeric_reads_and_metrics(
        args['input'],
        args['outbam'],
        args['outcsv'],
        args['sample_id'],
        perc_threshold=args['perc_threshold'],
        mapping_quality=args['mapping_quality'],
        length_threshold=args['telomere_length_threshold']
    )


if __name__ == '__main__':
    main()
