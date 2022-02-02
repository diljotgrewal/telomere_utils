import pandas as pd


def merge(infiles, outfile):
    infiles = [pd.read_csv(filepath) for filepath in infiles]
    infiles = pd.concat(infiles)

    infiles['start'] = infiles['start'].astype("Int64")
    infiles['end'] = infiles['end'].astype("Int64")

    if 'telomere_start' in list(infiles.columns.values):
        infiles['telomere_start'] = infiles['telomere_start'].astype("Int64")
    if 'telomere_end' in list(infiles.columns.values):
        infiles['telomere_end'] = infiles['telomere_end'].astype("Int64")


    infiles.to_csv(outfile, index=False, na_rep='NA')
