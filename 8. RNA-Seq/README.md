# timecourse RNA-seq experiment between CAM and C3 Tillandsia

This folder reports the steps undertaken to infer differential gene expression across 6 time points within a 24 hour period between a CAM Tillandsia (*Tillandsia fasciculata*) and a C3 plant (*Tillandsia leiboldiana*).

## Raw data processing

Data trimming was performed using [AdapterRemoval v.2.3.1](https://adapterremoval.readthedocs.io/en/stable/installation.html) with the bash script `01_trim.sh`. Post-trimming, samples contained 58 - 97 million pairs, with three outlier samples containing about 150 million pairs.

## Choosing a reference genome for mapping

We mapped a subset of our data to both Tillandsia genome assemblies and also to *A. comosus* to infer possible mapping bias and choose the better reference genome for this analysis.

We subsampled samples A,C and D for Tlei (D is individual used for the assembly) and A,C and F for Tfas (F is individual for the assembly), at times 01, 05, 13 and 17. This results in a total of 24 samples used in the test: 2 species x 3 samples x 4 timepoints.

I extracted random reads from each sample to 10 % of its original size using [seqkit](https://bioinf.shenwei.me/seqkit/) sample in the bash script `02_subsample.sh`. These were then mapped to each genome using [STAR](https://github.com/alexdobin/STAR) with `03_map.sh`. Mapping statistics were collected from each log file with the bash script `04_collect_mapping_stats.sh` and assessed with the R script `05_Assessment_mapping_bias.R`. The assessment showed that RNA data maps far better to the Tillandsia genomes than to *A. comosus* (around 90 % to around 50 %). Mapping bias seems much stronger when mapping to *T. leiboldiana* than to *T. fasciculata*. This seems to be rather due to unmapped reads than due to multimapper reads. Based on these results, we decided to map the full dataset to *T. fasciculata*. The raw mapping results of the test can be found in the Results folder as `mapping_stats.testset.txt`.

## Mapping full dataset

All 72 samples were  mapped to the *T. fasciculata* genome with the same script as above. General mapping stats were obtained with MultiQC and also with the script above. Further assessment can be found in `05_Assessment_mapping_bias.R`. The mapping stats are available in Results as `mapping_stats.full.txt`.

## Obtaining counts per gene

Counts per exon were computed using [featureCounts](https://subread.sourceforge.net/) with the bash script `06_counts_exon.sh`. Per-exon counts were summed up over all exons per gene with the python script `07_script_sum-up_featurecounts_exon.py`. The count data was inspected with Principal Component Analysis in the R sheet `08_Inspect_count_data.R`. The same procedure was repeated for reads mapped to the *T. leiboldiana* genome.
The raw count tables for both species are available in Results as `counts.mapped_to_Tfas.exons.txt` and `counts.mapped_to_Tlei.exons.txt`.

## Co-expression analysis with maSigPro

We ran [maSigPro](https://www.bioconductor.org/packages/devel/bioc/vignettes/maSigPro/inst/doc/maSigProUsersGuide.pdf) in the script `09_maSigPro_script_Tfas-Tlei.R`.

The steps for differential gene expression are as follows: we normalize the data in EdgeR and remove all genes with a mean(cpm) < 1. Then, we create the design matrix, containing the time points, replicates and the experimental groups (species). MaSigPro then identifies differentially expressed (DE) genes across time and species.

GO term enrichment was performed on all DE genes and each cluster separately with the automated script in 6. The raw GO term enrichment results for both species are available in Results. Expression curves of DE genes were drawn using `10_Script_Expression_curves_modules_maSigPro.R`. Heatmap was plotted with script `11_Expression_Heatmap_FamilySpecific.z-scores.R`.
