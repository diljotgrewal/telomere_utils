import pandas as pd


def merge(infiles, outfile):
    infiles = [pd.read_csv(filepath) for filepath in infiles]
    infiles = pd.concat(infiles)

    infiles['start_1'] = infiles['start_1'].astype("Int64")
    infiles['end_1'] = infiles['end_1'].astype("Int64")
    infiles['telomere_start_1'] = infiles['telomere_start_1'].astype("Int64")
    infiles['telomere_end_1'] = infiles['telomere_end_1'].astype("Int64")

    infiles['start_2'] = infiles['start_2'].astype("Int64")
    infiles['end_2'] = infiles['end_2'].astype("Int64")
    infiles['telomere_start_2'] = infiles['telomere_start_2'].astype("Int64")
    infiles['telomere_end_2'] = infiles['telomere_end_2'].astype("Int64")

    infiles.to_csv(outfile, index=False, na_rep='NA')
