# PopPUNK-scripts
Scripts used to help produce extra analysis for the [PopPUNK](https://github.com/johnlees/PopPUNK) manuscript

Some more general scripts can be found in the /scripts directory of the PopPUNK repository, and are documented on [readthedocs](https://poppunk.readthedocs.io/en/latest/scripts.html).

## allele_dists.py

Calculate the pairwise number of allele differences between isolates, using a tab-separated file of MLST-like calls. First row is a header with gene names, first column is the sample names. Entries are allele identifiers e.g.:
```
Sample  adk     fumC    gyrB    icd     mdh     purA    recA
11679_5#87      4       26      2       25      5       5       19
11679_4#35      53      40      47      13      36      28      29
11679_4#85      36      24      9       13      17      11      25
11791_3#27      13      24      19      14      23      1       10
11657_6#42      13      43      9       37      17      37      25
```

Calculate distances with:
```
python allele_dists.py mlst_calls.tsv mlst.dists.csv
```
where `mlst.dists.csv` is the output. This output can be used in `mlst_clusters.R` to define clusters using the PopPUNK network structure.

## mlst_clusters.R

Using a pairwise distance matrix of number of allele changes (from `allele_dists.py`) define clusters in the same way as PopPUNK, using a manually specified cutoff.

## tree_clusters.R

1) Compare clusters to a phylogeny to count polyphyly
2) Compare clusters to a pairwise distance matrix (of SNPs) to calculate and plot between- and within-cluster distances.

For 2) the SNP distance matrices were calculated with https://github.com/gtonkinhill/pairsnp

## sketch_time.R

Plots time and memory taken for sketches of different sizes, using the file `sketch_times.txt`. This was produced using gnu `time`.
