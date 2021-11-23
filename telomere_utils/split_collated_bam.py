import os


def getheader(bam):
    header = []

    for line in open(bam, 'rt'):
        if line.startswith("@"):
            header.append(line)
        else:
            break
    return header


def getfile(outdir, filenumber, header):
    filepath = os.path.join(outdir, '{}.sam'.format(filenumber))
    curr_file = open(filepath, 'wt')
    for line in header:
        curr_file.write(line)

    return curr_file


def split_bam(filepath, outdir, lines_per_file=1e6):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    header = getheader(filepath)

    bamreader = open(filepath, 'rt')

    filenumber = 0
    lines_written = 0

    currfile = getfile('outdir', filenumber, header)

    prev_read_id = None

    for line in bamreader:
        if line.startswith('@'):
            continue

        if prev_read_id is None:
            prev_read_id = line.split()[0]

        curr_read_id = line.split()[0]

        if curr_read_id == prev_read_id:
            currfile.write(line)
            lines_written += 1
            continue

        if lines_written < lines_per_file:
            currfile.write(line)
            lines_written += 1
            prev_read_id = curr_read_id
        else:
            filenumber += 1
            lines_written = 0
            currfile.close()
            currfile = getfile('outdir', filenumber, header)
            currfile.write(line)
            lines_written += 1
            prev_read_id = curr_read_id
