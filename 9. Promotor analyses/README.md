# Promotor Analyses

This section lists files and scripts used to study circadian-clock related transcription factor binding motifs in upstream regions of DE and CAM-related genes.

## Count motifs

The script `count_motifs.py` will extract the upstream regions (length is user-defined, we worked with 2-kb lengths) of genes using a gff file and fasta file as input. It will then search for motifs in these region and save the counts for each gene. It also calculates the observed and expected frequency of specific motifs in the genome based on GC content.

## Test significance of motif frequency

The script `significance.py` takes a per-gene list of motif counts for background genes (not-DE) and for targeted genes (DE). It will test for normality of these counts and use a parametric or non-parametric test to see if the motif frequency is significantly different between background and target genes. It will also perform a chi-square test with a 2x2 contingency table to test if target genes are enriched for presence of a motif compared to the background genes.
