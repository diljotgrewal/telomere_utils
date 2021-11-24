from collections import defaultdict
from operator import itemgetter

import pandas as pd
import pysam


def get_chrom_lengths(bamfile):
    chromosomes = bamfile.references
    lengths = bamfile.lengths
    lengths = {name: length for name, length in zip(chromosomes, lengths)}
    return lengths


def get_overlapping_bin(bamfile, chrom, pos, binsize=1000):
    chr_lengths = get_chrom_lengths(bamfile)

    start = (pos // binsize) * binsize
    end = start + binsize
    if end > chr_lengths[chrom]:
        end = chr_lengths[chrom]

    return (chrom, int(start), int(end))


def build_germline_ref(bamfile, telomeres, binsize=1000):
    """
    telomere length max is 150 (read length)
    so check the bin that covers start and end of telomere
    keeping a boolean per bin as proof of telomere
    """

    telomeres = pd.read_csv(telomeres)

    outdata = defaultdict(bool)

    for i, row in telomeres.iterrows():

        if not pd.isna(row['telomere_start_1']):
            chromosome_1 = row['chromosome_1']
            telomere_start_1 = int(row['start_1'] + row['telomere_start_1'])
            telomere_end_1 = int(row['start_1'] + row['telomere_end_1'])
            binval_start = get_overlapping_bin(bamfile, chromosome_1, telomere_start_1, binsize=binsize)
            binval_end = get_overlapping_bin(bamfile, chromosome_1, telomere_end_1, binsize=binsize)

            outdata[binval_start] = 1
            outdata[binval_end] = 1
        else:
            assert not pd.isna(row['telomere_start_2'])
            chromosome_2 = row['chromosome_2']
            telomere_start_2 = int(row['start_2'] + row['telomere_start_2'])
            telomere_end_2 = int(row['start_2'] + row['telomere_end_2'])

            binval_start = get_overlapping_bin(bamfile, chromosome_2, telomere_start_2, binsize=binsize)
            binval_end = get_overlapping_bin(bamfile, chromosome_2, telomere_end_2, binsize=binsize)

            outdata[binval_start] = 1
            outdata[binval_end] = 1

    return outdata


def build_counts(bamfile, telomeres, output, binsize=1000):
    """
    telomere length max is 150 (read length)
    so check the bin that covers start and end of telomere
    keeping a boolean per bin as proof of telomere
    """
    chr_lengths = get_chrom_lengths(bamfile)

    telomeres = pd.read_csv(telomeres)

    telomeres = telomeres[telomeres['germline'] == False]

    outdata = defaultdict(int)

    for i, row in telomeres.iterrows():

        if not pd.isna(row['telomere_start_1']):
            chromosome_1 = row['chromosome_1']
            telomere_start_1 = int(row['start_1'] + row['telomere_start_1'])
            telomere_end_1 = int(row['start_1'] + row['telomere_end_1'])
            binval_start = get_overlapping_bin(bamfile, chromosome_1, telomere_start_1, binsize=binsize)
            binval_end = get_overlapping_bin(bamfile, chromosome_1, telomere_end_1, binsize=binsize)
            if binval_start == binval_end:
                outdata[binval_start] += 1
            else:
                outdata[binval_start] += 1
                outdata[binval_end] += 1
        else:
            assert not pd.isna(row['telomere_start_2'])
            chromosome_2 = row['chromosome_2']
            telomere_start_2 = int(row['start_2'] + row['telomere_start_2'])
            telomere_end_2 = int(row['start_2'] + row['telomere_end_2'])
            binval_start = get_overlapping_bin(bamfile, chromosome_2, telomere_start_2, binsize=binsize)
            binval_end = get_overlapping_bin(bamfile, chromosome_2, telomere_end_2, binsize=binsize)
            if binval_start == binval_end:
                outdata[binval_start] += 1
            else:
                outdata[binval_start] += 1
                outdata[binval_end] += 1

    sorted_counts = sorted(outdata.items(), key=itemgetter(1))
    outdata = []
    for binval, count in sorted_counts:
        chromosome, start, end = binval
        data = {
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'count': count,
            'at_chromosome_end': False
        }
        if start < 10000 or end > (chr_lengths[chromosome] - 10000):
            data['at_chromosome_end'] = True

        outdata.append(data)

    outdata = pd.DataFrame(outdata)
    outdata.to_csv(output, index=False, na_rep='NA')


def filter_telomeres(tumour_telomeres, normal_data, output, bamfile, binsize=1000):
    tumour_telomeres = pd.read_csv(tumour_telomeres)

    alldata = []

    for i, row in tumour_telomeres.iterrows():
        row['germline'] = False

        if not pd.isna(row['telomere_start_1']):
            chromosome_1 = row['chromosome_1']
            telomere_start_1 = int(row['start_1'] + row['telomere_start_1'])
            telomere_end_1 = int(row['start_1'] + row['telomere_end_1'])

            binval_start = get_overlapping_bin(bamfile, chromosome_1, telomere_start_1, binsize=binsize)
            binval_end = get_overlapping_bin(bamfile, chromosome_1, telomere_end_1, binsize=binsize)

            if normal_data[binval_start] or normal_data[binval_end]:
                row['germline'] = True

        else:
            assert not pd.isna(row['telomere_start_2'])
            chromosome_2 = row['chromosome_2']
            telomere_start_2 = int(row['start_2'] + row['telomere_start_2'])
            telomere_end_2 = int(row['start_2'] + row['telomere_end_2'])

            binval_start = get_overlapping_bin(bamfile, chromosome_2, telomere_start_2, binsize=binsize)
            binval_end = get_overlapping_bin(bamfile, chromosome_2, telomere_end_2, binsize=binsize)
            if normal_data[binval_start] or normal_data[binval_end]:
                row['germline'] = True

        alldata.append(row)

    alldata = pd.DataFrame(alldata)
    alldata.to_csv(output, index=False, na_rep="NA")


def get_overlap(bamfile, normal_telomeres, tumour_telomeres, outfile, bin_counts, binsize=1000):
    bamfile = pysam.AlignmentFile(bamfile)

    normal_data = build_germline_ref(bamfile, normal_telomeres, binsize=binsize)

    filter_telomeres(tumour_telomeres, normal_data, outfile, bamfile, binsize=binsize)

    build_counts(bamfile, outfile, bin_counts, binsize=binsize)
