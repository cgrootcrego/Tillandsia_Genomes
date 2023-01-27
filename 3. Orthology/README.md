
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

The resulting per-gene orthology table is available as `OrthologyInfo.per-gene.AcomTfasTlei.txt`
