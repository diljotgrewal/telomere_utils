#!/usr/bin/env python

import argparse
import pandas as pd
import pysam


def _revcomp(seq):
    seq1 = seq.translate(str.maketrans("AGCTagct", "TCGAtcga"))
    seq2 = seq1[::-1]
    return seq2


def _find_telomere(query_sequence):
    kmers = ['TTAGGG', 'TCAGGG', 'TGAGGG', 'TTGGGG']

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


def _get_telomeric_pos(read, perc_threshold=0.85, mapping_quality=30, telomere_length_threshold=36):
    # good mapping qual reads cannot be telomeric
    if read.mapping_quality > mapping_quality:
        return

    query_sequence = read.query_sequence
    if read.is_reverse:
        query_sequence = _revcomp(query_sequence)

    tot_len, telomere_end = _find_telomere(query_sequence)

    # no telomere found
    if telomere_end is None:
        return

    # telomere not long enough
    if tot_len < telomere_length_threshold:
        return

    percentage = tot_len / telomere_end

    # most of the telomeric section is not telomeric enough
    if percentage < perc_threshold:
        return

    # we revcomped the read when calculating positions. so count it in reverse here
    if read.is_reverse:
        end = len(query_sequence)
        start = end - telomere_end
    else:
        start = 0
        end = telomere_end

    return [start, end]


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


def _get_csv_data(r1, r2, telo_r1, telo_r2, sample_id, bamfile):
    if r1.is_read1:
        rend1 = r1
        rend1_telo = telo_r1
        rend2 = r2
        rend2_telo = telo_r2
    else:
        assert r2.is_read1
        rend1 = r2
        rend1_telo = telo_r2
        rend2 = r1
        rend2_telo = telo_r1

    if bamfile.get_reference_name(rend1.reference_id) is None:
        return

    if bamfile.get_reference_name(rend2.reference_id) is None:
        return

    data = {
        'read_id': rend1.query_name,
        'sample_id': sample_id,

        'strand_1': '-' if rend1.is_reverse else '+',
        'chromosome_1': bamfile.get_reference_name(rend1.reference_id),
        'start_1': rend1.reference_start,
        'end_1': rend1.reference_end,
        'telomere_start_1': None if rend1_telo is None else rend1_telo[0],
        'telomere_end_1': None if rend1_telo is None else rend1_telo[1],

        'strand_2': '-' if rend2.is_reverse else '+',
        'chromosome_2': bamfile.get_reference_name(rend2.reference_id),
        'start_2': rend2.reference_start,
        'end_2': rend2.reference_end,
        'telomere_start_2': None if rend2_telo is None else rend2_telo[0],
        'telomere_end_2': None if rend2_telo is None else rend2_telo[1],
    }

    return data


def extract_telomeric_reads_and_metrics(
        infile, outbam, outcsv, sample_id, perc_threshold=0.85, mapping_quality=30,
        telomere_length_threshold=36
):
    bamfile = pysam.AlignmentFile(infile, "rb")
    outfile = pysam.AlignmentFile(outbam, "wb", template=bamfile)

    csvdata = []

    for r1, r2 in _get_read_pairs(bamfile):
        telo_r1 = _get_telomeric_pos(
            r1, perc_threshold=perc_threshold, mapping_quality=mapping_quality,
            telomere_length_threshold=telomere_length_threshold
        )
        telo_r2 = _get_telomeric_pos(
            r2, perc_threshold=perc_threshold, mapping_quality=mapping_quality,
            telomere_length_threshold=telomere_length_threshold
        )

        if telo_r1 is None and telo_r2 is None:
            continue

        outfile.write(r1)
        outfile.write(r2)

        read_pair_data = _get_csv_data(r1, r2, telo_r1, telo_r2, sample_id, bamfile)

        if read_pair_data is not None:
            csvdata.append(read_pair_data)

    df = pd.DataFrame(csvdata)

    if df.empty:
        columns = [
            'read_id', 'sample_id',
            'strand_1', 'chromosome_1', 'start_1', 'end_1', 'telomere_start_1', 'telomere_end_1',
            'strand_2', 'chromosome_2', 'start_2', 'end_2', 'telomere_start_2', 'telomere_end_2',
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
        telomere_length_threshold=args['telomere_length_threshold']
    )


if __name__ == '__main__':
    main()