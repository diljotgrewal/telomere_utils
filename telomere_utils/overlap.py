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


def build_germline_ref(bamfile, telomeres, chromosome=None, binsize=1000):
    """
    telomere length max is 150 (read length)
    so check the bin that covers start and end of telomere
    keeping a boolean per bin as proof of telomere
    """

    telomeres = pd.read_csv(telomeres)

    if chromosome:
        telomeres = telomeres[telomeres['chromosome'] == chromosome]

    outdata = defaultdict(bool)

    for i, row in telomeres.iterrows():
        chrom = row['chromosome']
        telomere_start = int(row['start'] + row['telomere_start'])
        telomere_end = int(row['start'] + row['telomere_end'])
        binval_start = get_overlapping_bin(bamfile, chrom, telomere_start, binsize=binsize)
        binval_end = get_overlapping_bin(bamfile, chrom, telomere_end, binsize=binsize)

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

    telomeres = pd.read_csv(telomeres, dtype={"chromosome": "str"})

    telomeres = telomeres[telomeres['germline'] == False]

    outdata = defaultdict(int)

    for i, row in telomeres.iterrows():

        chromosome = row['chromosome']
        telomere_start = int(row['start'] + row['telomere_start'])
        telomere_end = int(row['start'] + row['telomere_end'])
        binval_start = get_overlapping_bin(bamfile, chromosome, telomere_start, binsize=binsize)
        binval_end = get_overlapping_bin(bamfile, chromosome, telomere_end, binsize=binsize)
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

    if outdata.empty:
        columns = ['chromosome', 'start', 'end', 'count', 'at_chromosome_end']
        outdata = pd.DataFrame(columns=columns)

    outdata.to_csv(output, index=False, na_rep='NA')


def filter_telomeres(tumour_telomeres, normal_data, output, bamfile, chromosome=None, binsize=1000):
    tumour_telomeres = pd.read_csv(tumour_telomeres)

    if chromosome:
        tumour_telomeres = tumour_telomeres[tumour_telomeres['chromosome'] == chromosome]

    alldata = []

    for i, row in tumour_telomeres.iterrows():
        row['germline'] = False

        chrom = row['chromosome']

        telomere_start = int(row['start'] + row['telomere_start'])
        telomere_end = int(row['start'] + row['telomere_end'])

        binval_start = get_overlapping_bin(bamfile, chrom, telomere_start, binsize=binsize)
        binval_end = get_overlapping_bin(bamfile, chrom, telomere_end, binsize=binsize)

        if normal_data[binval_start] or normal_data[binval_end]:
            row['germline'] = True

        alldata.append(row)

    alldata = pd.DataFrame(alldata)

    if alldata.empty:
        columns = [
            'read_id', 'sample_id', 'strand', 'chromosome', 'start', 'end',
            'telomere_start', 'telomere_end', 'reverse_complement', 'readend',
            'germline'
        ]
        alldata = pd.DataFrame(columns=columns)

    alldata.to_csv(output, index=False, na_rep="NA")


def get_overlap(bamfile, normal_telomeres, tumour_telomeres, outfile, bin_counts, chromosome, binsize=1000):
    bamfile = pysam.AlignmentFile(bamfile)

    normal_data = build_germline_ref(bamfile, normal_telomeres, binsize=binsize, chromosome=chromosome)

    filter_telomeres(tumour_telomeres, normal_data, outfile, bamfile, binsize=binsize, chromosome=chromosome)

    build_counts(bamfile, outfile, bin_counts, binsize=binsize)
