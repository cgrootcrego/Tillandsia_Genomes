# Visualizing gene, repetitive and GC content in *T. fasciculata* and *T. leiboldiana* genomes

Here are all the steps to reproduce the circular figure showing window-based repetitive, gene and GC content for each scaffold of both species' genome assemblies, and comparisons of GC, TE and gene content in syntenic scaffolds.

## Calculating gene content

Gene content per 1 MB window was calculated both as the number of genes starting per window and the proportion of bases being in genic regions. This was done for all orthologous genes, single-copy and multicopy genes for both species in the Rscript `Calculate_Gene_content.R`.

## Calculating repetitive content

Using the masked version of the assemblies (after detailed annotation, see 2a. TE annotation), the percentage of masked bases per 1 MB windows was calculated with the python script `script_calculate_repetitive_content_perwindow_maskedfasta.py`:

````
    python calculate_repetitive_content_perwindow_maskedfasta.py Tillandsia_leiboldiana_26_scaffolds.fasta.masked 1000000
````

## Calculating GC content

With the script `script_calculate_GC_content_perwindow.py`, which is very similar to the one for repetitive content, I calculated GC content in 1 MB windows for each genome assembly:

````
python script_calculate_GC_content_perwindow.py Tillandsia_fasciculata_25_scaffolds.fasta 1000000
````

## Visualizing per-scaffold gene, GC and repetitive content

The circular figures showing per-scaffold genic and repetitive content were made using [circlize](https://jokergoo.github.io/circlize_book/book/) in R with the script `Script_Circular_figure_Gene-TE-GC_content.R`.

## Calculating abundances of TE classes for the whole genome and per scaffold

Per-scaffold TE content was calculated in the Rscript `Analyze_per-scaffold_TE_abundances.R`. This showed that while large scaffolds had on average 65 % TE content in T. fasciculata and 77 % in T. leiboldiana, small scaffolds consisted on average of 94 % TEs. 35 % of all TE content can be found on small scaffolds, which make up 25 % of the assembly. This, together with evidence from orthofinder showing that almost all orthologous genes are on the main scaffolds, shows that small scaffolds in the assembly are mostly TE remnants that were not incorporated in the main assembly.

## Calculating per scaffold genic / repetitive content ratios

Since visually it is clear that TE content is larger in T. leiboldiana, we computed genic-to-repetitive content ratios for each contig for both species, to have a numeric comparison. We already have the number of repetitive basepairs per scaffold from the in-depth analyses of TE content (see above). To obtain total numbers of genic basepairs per scaffold, I counted basepairs marked as exon from the gff files of both genomes with the script `script_genic_proportion_perscaff.py`. The final per-scaffold proportions of exonic and TE content can be found in `Tillandsia_fasciculata.perscaff.genic_and_TE_content.txt` and `Tillandsia_leiboldiana.perscaff.genic_and_TE_content.txt`. The significance of the difference in TE to exonic content ratio between *T- fasciculata* and *T. leiboldiana* was calculated in `Sign_TE_to_exonic_perScaffold.R`.

## Correlations between Gene, GC and TE content across the genome

The correlation between gene, TE and GC content was calculated in `Correlations_gene_te_gc_content.R`

## Comparisons of gene, TE and GC content in syntenic scaffold triplets

See 5a.
