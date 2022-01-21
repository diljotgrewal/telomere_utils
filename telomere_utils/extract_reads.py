#!/usr/bin/env python

import argparse
import pandas as pd
import pysam


def _find_telomere_end_in_seq(query_sequence, kmers):
    tot_len = 0
    end_pos = None
    i = 0
    while i < (len(query_sequence) - 5):
        if any(tr == query_sequence[i:(i + 6)] for tr in kmers):
            tot_len += 6
            end_pos = i + 6
            i += 6
        else:
            i += 1

    return tot_len, end_pos


def _keep_telomere(telo_len, end_pos, perc_threshold=0.85, length_threshold=36):
    # no telomere found
    if end_pos is None:
        return False

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

    telo_len, telo_end = _find_telomere_end_in_seq(query_sequence, kmers)
    if _keep_telomere(telo_len, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append((0, telo_end, False))

    telo_len, telo_end = _find_telomere_end_in_seq(query_sequence, kmers_revcomp)
    if _keep_telomere(telo_len, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append((0, telo_end, True))

    telo_len, telo_end = _find_telomere_end_in_seq(query_sequence[::-1], kmers)
    if _keep_telomere(telo_len, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append((len(query_sequence) - telo_end, len(query_sequence), False))

    telo_len, telo_end = _find_telomere_end_in_seq(query_sequence[::-1], kmers_revcomp)
    if _keep_telomere(telo_len, telo_end, perc_threshold=perc_threshold, length_threshold=length_threshold):
        positions.append((len(query_sequence) - telo_end, len(query_sequence), True))

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
        'telomere_start': telomere[0],
        'telomere_end': telomere[1],
        'reverse_complement': telomere[2],
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
