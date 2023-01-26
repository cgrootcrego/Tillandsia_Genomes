
# Calling orthogroups between *T. fasciculata*, *T. leiboldiana* and *A.comosus*

We called orthogroups between three bromeliad gene model sets using [OrthoFinder](https://github.com/davidemms/OrthoFinder).

# First run: all T. fasciulata and T. leiboldiana gene models included, regardless of location

The first orthofinder run involved all called gene models and was run to see the distribution of genes across the assembly, making distinctions between one-to-one orthogroups and orthogroups with other relationships. The idea behind this was that, if most gene models, especially gene models with one-to-one relationships are on main scaffolds, it confirms that the main scaffolds represent the chromosomes and we can limit downstream analyses to these sequences. Alternatively, if many genemodels are also present on smaller scaffolds, this may indicate fragmentation in our assembly and an anchoring step may be necessary.

**1) Input sequences**

The following peptide sequence datasets were provided for running orthofinder:

  - T. fasciculata: MAKER Protein sequences (See 2b. Gene annotation). Contains all gene models, including those that didn't blast in functional annotation (34,886 sequences)
  - T.leiboldiana: MAKER Protein sequences. Again, the full set was used including sequences that didn't blast (38,180).
  - Ananas comosus: [Acomosus F153 peptide sequences](http://www.life.illinois.edu/ming/LabWebPage/Downloads.html) (27,024).

**2) Running Orthofinder**

Orthofinder was run with version 2.4.1. with the following command:

    /apps/orthofinder/2.4.1/orthofinder \
	  -f [run_directory] \
	  -t 48

 **3) Compiling orthology and annotation results**

Orthology and annotation info per gene was compiled from the N0.tsv file and gff file with the scripts  `script_make_per_gene_og_table_calculate_counts.py` and`script_compile_gff_info_og_table_all_species.py`.

**4) Analyzing the spatial and length distribution of global orthogroups**

Genes of different species were placed in separate files:

    grep "Tlei" orthogroups_Tfas-Tlei-Acom.full-set.no_Acom_specific_og.all_info.txt > Tlei_orthology_info.txt
    grep "Tfas"  orthogroups_Tfas-Tlei-Acom.full-set.no_Acom_specific_og.all_info.txt > Tfas_orthology_info.txt

Scaffolds lengths were added with `add_scaffold_lengths.py`.

The actual analysis was done in the Rscript `Analyze_global_orthogroups.R`.

# Second run: Only gene models on main scaffolds (> 1 Mb) for T.fas and T. lei, all gene models of A.comosus

With the aim of studying synteny and gene family evolution, we reran Orthofinder for genes on the main scaffolds of the Tfas and Tlei assemblies.

**1) Input sequences**

I filtered out shorter isoforms and also peptide sequences < 40 amino acids.
The longest isoform was selected using the script `select_longest_isoform.py`.
Peptide sequences < 40 AA were filtered out with the following lines of code:

	# select peptides longer than 40 AA
	awk '$2 > 40 {print $0}' \
	Tillandsia_fasciculata_v1.25chrom.longest_isoforms.proteins.fasta.fai \
	| cut -f 1 > proteins.w.min40AA

	seqkit grep -f proteins.w.min.40.AA \
	Tillandsia_fasciculata_v1.25chrom.longest_isoforms.proteins.fasta > \
	Tillandsia_fasciculata_v1.25chrom.longest_isoforms.min40AA.proteins.fasta

**2) Running OrthoFinder**

Orthofinder was run with the same command as above.

**3) Compiling orthology with gff info and functional annotations, plus a bit more filtering**

Remove orthogroups unique to *A. comosus*:

`awk '!($3 == 0 && $4 == 0) {print $0}' Orthogroups.GeneCount.tsv > orthogroup_counts_no_Acom_specific_og.txt`

Gff info and orthogroup info were compiled with `script_compile_gff_info_with_og_table_all_species.py`:

All orthogroups containing genes with TE description were removed, by searching for the words:

    transposable
    transposon
    transposase
    Transposon
    Gag-Pol polyprotein
    Pro-Pol polyprotein
    virus

I obtained the IDs of these orthogroups and then removed them from the table :

	grep -f tes orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.txt | cut -f 7 | \
	sort | uniq > orthogroups_containing_TES
    	grep -v -f orthogroups_containing_TES \
	orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.txt > orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt

This removed a total of 6042 genes.
The file `orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt` is the final set of curated orthogroups that has been used in downstream analysis. It contains 70,963 genes and 19,101 orthogroups (26325 for Tfas, 23584 for Tlei, 21055 for Aco). Additional statistics on the orthology analysis can be found [here](https://docs.google.com/spreadsheets/d/1hv_Oe6MvV1fBuBkGIhvM1ihDwYviuY-LcYUGYKnIwNQ/edit#gid=324678736).

# Re-evaluation of multi-copy families

I decided to re-evaluate our curated orthogroups based on other criteria than the ones previously used (annotation, coverage), to further ascertain to validity of multi-copy genes. For this, I decided to calculate per gene the number of exons, the length, absence of start and stop codon, and expression in both species. Then, I made comparative measures by taking the ratio of the largest member in the gene family (of the same species) to each gene both for full exonic length and the number of exons. This should give us an idea how many multi-copy genes are "incorrect" or incomplete and are perhaps transposed or falsely assembled genes.

The measures were taken from many different input files (gff, fasta sequence, orthology and expression) and compiled in the script `script_checklist_curated_ogs.py`. I later updated this script to distinguish between one stopcodon, multiple stopcodons or no stopcodons. I also incorporated a metric of the percentage of exons expressed per gene.

Later even, I decided to be more stringent with expression rates. previously, an exon that had 1 read in only one sample would already count as expressed. This time, I decided to take a harsher threshold. I made this choice after looking at the distribution of mean CPM per exon for both genomes and species. The median was always around 0.6, meaning that 50 % of the exons has less than 600,000 reads mapped on average across samples in both genomes and species. Our cutoff for DGE analysis is quite stringent (we removed about half of the genes at the time), but here I just want to see if there is some kind of expression. Expression can vary for many reasons, so I didn't want to choose such a high threshold to decide if a gene is dead or not. In the end, I decided that all exons with an average CPM across samples lower than 0,001 (1000 counts in total) would be marked as not expressed. This will approximately remove between 10-20 % of exons. The new changes were applied in `script_checklist_curated_ogs.updated062022.py`.

I then assessed the checklist in `Script_Checklist_ortholog_assessment.R`. It demonstrates that most genes with an absent stop/startcodon are in fact fully expressed. So, though especially in T. leiboldiana a large proportion of gene models (~ 40 %) are missing either a start/stopocodon or both, they seem to be functional. If applying a filter of not-expressed + absence of start/stopcodon, only 3000 genes would be remove din Tlei and ~ 1000 in Tfas.

In the end we decided to add a tag to genes declaring if they were ROBUST or not. A gene is regarded as ROBUST if all exons are expressed, OR if it has start/stopcodon and a minimum length of 50 amino-acids. The tag was added to the final full genome GFF file with the script `Script_add_robustness_expression_gff.py`

----
I also evaluated whether extra copies in T. fasciculata are also expressed in T. leiboldiana, as this seemed to be the case for some of the DE families. I did this using the checklist and the list of Tfasciculata gene families (so all orthologs, including the differentially expressed gene) that are differentially expressed and larger in Tfas. I counted the number of genes that were expressed in Tfas but not in Tlei, this was 62 genes. I also counted the total difference in counts between Tlei and Tfas in this family, this was 532. So, it seems that a lot of copies that supposedly don't exist in T. leiboldiana are still mapping to T. leiboldiana RNA-seq.
