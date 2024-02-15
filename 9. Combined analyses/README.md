# Investigations of DE genes, gene duplication and TE insertions

Here is the pipeline for all combinatory analyses performed on the output of several previous analyses together.

## Distribution of DE genes

Stats on DE gene distribution across both genomes were calculated in `GetDistributionStats_DEgenes.R`. Statistical test results are in `Test_statistics.DE_distribution.txt`.

## Gene duplication in DE genes

The proportion of DE orthogroups belonging to different relative size classes was obtained with `01_get_duplication_stats.sh`. The significance of the differences in proportion between all orthogroups and DE orthogroups was studied in `02_Significance_duplication.R`. Gene family expansion in families underlying enriched GO terms were visualized in `03_Script_GOterm_figure.R` using [GOplot](https://wencke.github.io/).

## TE insertions in DE Genes

We calculated the number of genes with a TE insertion in introns, both for the whole genome and for DE genes.
As input, we used the EDTA output file "*.EDTA.TEanno.gff3" and removed all unknown TEs by `grep -v "unknown" [file] > *.EDTA.TEanno.no-unknown.gff3`. We also used the gene annotation file and ran bedtools intersect with the script `05_bedtools_intersect.sh`. This creates two files, one containing the number of insertions for each gene, and another one containing the name of each TE insertion in each gene.

The significance of TE insertions in DE genes compared to the whole genome was then calculated in `06_Script_test_significance_TE_insertions_DE-nonDE.R`.
