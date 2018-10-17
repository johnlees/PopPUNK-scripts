#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from math import floor
import tqdm

def main():

    parser = argparse.ArgumentParser(description='Pairwise distances between MLST alleles')
    parser.add_argument('infile', type=str, help="Tab separated file containing alleles")
    parser.add_argument('outfile', type=str, help="Name for output file")
    args = parser.parse_args()

    alleles_in = pd.read_csv(args.infile, sep="\t", header=0, index_col=0, dtype=str)
    num_samples = len(alleles_in.index)
    dists = np.zeros((num_samples, num_samples), dtype=int)

    with tqdm.tqdm(total = int(0.5*num_samples*(num_samples-1))) as pbar:
        for i in range(num_samples):
            row1 = alleles_in.iloc[i, :].values
            for j in range(i+1, num_samples):
                row2 = alleles_in.iloc[j, :].values
                diffs = np.sum(np.not_equal(row1, row2))
                dists[i,j] = diffs
                dists[j,i] = diffs
                pbar.update(1)

    dists_out = pd.DataFrame(dists, index=alleles_in.index, columns=alleles_in.index)
    dists_out.to_csv(args.outfile)
    print("Done\n")

if __name__ == "__main__":
    main()


