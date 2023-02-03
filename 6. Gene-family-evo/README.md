# Inferring Gene Family Evolution between *T. fasciculata* and *T. leiboldiana*

Using the inferred gene models and orthogroups between the Tillandsia genomes, we studied gene family size changes between *T. fasciculata* and *T. leiboldiana*.

## Detecting faulty gene models resulting from haplotigs in *T. fasciculata*

*T. fasciculata* is a wide-spread species, often regarded as a species complex, with very high levels of heterozygosity. Unsurprisingly, the accession used for our *T. fasciculata* reference genome was also more heterozygous than desired for genome assembly. This is especially an issue for inferring gene family sizes, since haplotigs will be considered separate copies.
Using 50x whole genome sequencing data of our reference genome accession, looked at coverage across the genome, and more importantly across gene models. The idea here is that legitimate gene copies should have similar coverage to the rest of the genome. However, erroneous duplications that were inferred for example due to the existence of haplotigs, will have a lower coverage than the rest of the genome, as reads from one gene will be distributed across two genic regions in the genome.

**1) Alignment of 50x WGS data to *T. fasciculata* reference:**

Alignment of the raw data to the *T. fasciculata* genome was done using [Bowtie2](https://github.com/BenLangmead/bowtie2). We chose Bowtie2 because its manual explicitly mentions that reads aligning equally well in two different locations will be assigned to one of these two locations at random. This is important, since we wanted to infer a decrease in coverage when gene duplicates are false, and this would be more difficult if reads are mapped at multiple locations. BWA doesn't explicitly specify this.

    # Genome indexing:
    bowtie2-build -f Tillandsia_fasciculata_25_scaffolds.fasta Tillandsia_fasciculata_25_scaffolds
    # alignment
    bowtie2 --very-sensitive-local \
	-x /proj/grootcrego/Genome_assemblies/fasciculata/final/Tillandsia_fasciculata_25_scaffolds.fasta \
	-1 /proj/grootcrego/Genome_assemblies/fasciculata/0_raw_data/illumina/Tfas_Illumina_50x_trimmed_pair1.fq \
	-2 /proj/grootcrego/Genome_assemblies/fasciculata/0_raw_data/illumina/Tfas_Illumina_50x_trimmed_pair2.fq \
	-S Tfas_50x_illumina_to_Tfas25chrom.sam -p 24
    # Transforming to bam and sorting
    samtools view -Sb Tfas_50x_illumina_to_Tfas25chrom.sam > Tfas_50x_illumina_to_Tfas25chrom.bam
    samtools sort Tfas_50x_illumina_to_Tfas25chrom.bam -o Tfas_50x_illumina_to_Tfas25chrom.sorted.bam

**2) Calculating per-gene median and average coverage**

Per-base coverage in genic regions of *T. fasciculata* was calculated using samtools depth. We first obtained a bedfile with the coordinates of ortholohous gene models:

    # Select only Tfas genes from file containing curated gene models
	grep "Tfasc_v1"  orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt > curated_Tfas_orthologues.txt
    # Obtain only gene ID
	cut -f 1 curated_Tfas_orthologues.txt > curated_Tfas_orthologues.IDonly.txt
    # Select curated genes from annotation GFF file
	grep -f curated_Tfas_orthologues.IDonly.txt /proj/grootcrego/Genome_assemblies/fasciculata/4_final_assembly/Tillandsia_fasciculata_v1.2.edited_allfeatures.gff > Tillandsia_fasciculata_v1.2.curated_orthologues_only.gff
    # Limit curated gff file only to mRNA entries (no exon, UTR, or CDS...)
	awk '$3 == "mRNA" {print $0}' Tillandsia_fasciculata_v1.2.curated_orthologues_only.gff > Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff
    # This curated GFF file was then converted to a BED format using gff2bed from bedops v. 2.4.37:  
    gff2bed < Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff > Tfas_ortholog_regions.gff.bed

Using this bedfile, we calculated per-base coverage with samtools. Coverage of all bases included in the bedfile was reported, including areas of zero coverage (-a):

    samtools depth -a -b $bedfile $bamfile -o $output

The name of the gene was then added at each position in the samtools depth output file using the python script `script_add_gene_names_to_cov_file.py`. Mean and median coverage were then computed per gene and compiled into a table together with orthology information with the script `script_calculate_mediancov_add_info.py`. The resulting tables are available as `Tfas_pergene_mediancov_and_orthoinfo.txt` and `Tlei_pergene_mediancov_and_orthoinfo.txt` for gene models of *T. fasciculata* and *T. leiboldiana* respectively.

## Correcting gene family sizes based on average coverage

Density plots of the average coverage per gene for different categories of genes were made using the Rscript `Assessing_multicopy_genemodels_cov.R`. These density plots showed that multi-copy gene models in *T.fasciculata* have a bimodal distribution of mean coverage - meaning that about half of the multicopy genes have a mean coverage around the genome-wide average, while the other half has a much lower mean coverage.

The average coverage for the full genome was calculated by running samtools depth -a without specifying regions and running the awk one-liner `awk '{ total += $3; count++ } END { print total/count }`. To obtain the median coverage over the full genome, I ran `cut -f 3 CovDepth_perbase_fullgenome.txt | sort -n | awk -f median.awk`. Median.awk was the following short script:

    #/usr/bin/env awk
    {
    count[NR] = $1;
    }
    END {
    if (NR % 2) {
	    print count[(NR + 1) / 2];
    } else {
	    print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;
    }
    }

We only applied gene family size correction to multi-copy gene families which did not belong to plastid, mitochondrial or ribosome classes. We removed these genes for correction, since their coverage doesn't necessarily reflect their copy number (multiple plastids / mitochondria in a cell, ribosomal genes have very complicated gene family evolution). For our strategy to remove such genes, please see 6a. Ientifying plastid, mitochondrial and ribosomal genes.

Gene family sizes were corrected by obtaining the total mean coverage in the orthogroup (sum of mean coverage of all *T.fasciculata* genes in the OG) and dividing this by the "expected mean coverage" (single-copy genes mean coverage x number of *T.fasciculata* genes in the orthogroup). I called this ratio the correction factor, which I multiplied by the original number of *T.fasciculata* genes in the orthogroup to obtain a corrected family size.

Note: Size corrections were made not based on the full genome average but on the average coverage of ancestral single-copy genes. The idea is that this category represents better the expected coverage in genic regions since the full genome coverage is extremely variable and largely influenced by non-genic regions.

In cases where the total mean coverage of the orthogroup didn't reach the single-copy mean average (our expectation of coverage for 1 gene), the correction usually brought the gene family size down to 1, but sometimes also to 0. In those cases, since the gene has been properly annotated and has an orthologous sequence in at least *T. leiboldiana* or *A.comosus*, we assumed the gene must exist and therefore corrected the size back to 1.

Size correction was performed on all multi-copy genes of *T. fasciculata* and *T. leiboldiana* in `Assessing_multicopy_genemodels_cov.R`.

## Studying gene family size differences in *T. fasciculata* and *T. leiboldiana*

After gene family size correction, the relationship of gene family size between species was explored in `Gene_family_evolution.R`.

To better understand which levels of orthogroup size changes between species were "substantial" or "deviating from expectations", I calculated Log ratios of the family sizes of *T. fasciculata* versus *T. leiboldiana* and ranked them. This showed that by selecting the top 2 % of observations, already the full multi-copy gene set was included. This is mostly because we have a very large proportion of orthogroups where family sizes have not changed between both species. My interpretation of this is that we have no real power to distinguish between multicopy families which are fast or slowly evolving. Therefore, it may be justified to analyze all multicopy groups as done below.

For downstream multicopy analyses (GO terms and genes of interest) I made a selection of non-unique gene families that had a family size > 1 in at least one of the Tillandsia species.

## GO term enrichment of multi-copy gene families

GO term enrichment was run on (i) all multicopy genes, (ii) a subset of genes belonging to families where copy number was larger in Tfas and (iii) on a subset of genes belonging to families where copy number was larger in Tlei. This can give us an idea of the direction in which certain functions have had gene family expansion. This was done with the automated script `script_GO_term)enrichment.R` and [TopGo](https://bioconductor.org/packages/release/bioc/html/topGO.html).

The top 100 GO terms in each category (BP, MF and CC) were stored alongside with their p-value and underlying orthogroups and genes, as `Enriched_GO_terms_Tfas-Tlei_ALL_multicopy_genes.txt`, `Enriched_GO_terms_multicopy_genes_larger-in-Tfas.txt`, and `Enriched_GO_terms_multicopy_genes_larger-in-Tlei.txt`

## Detecting "genes of interest" in the multicopy gene family pool

A reverse approach, where we search for genes of interest in our set of multicopy gene families, may give us a better idea on gene family evolution and key innovation traits. I made use of the ~ 700 genes of interest identified in previous work by [De la Harpe et al. 2020](https://onlinelibrary.wiley.com/doi/epdf/10.1111/pce.13847) and [Yardeni et al. 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13523) to see if a proportion of these genes were to be found in our mutli-copy gene family set.

Genes were selected from the Bromeliad1776 kit, which relies on the pineapple annotation (so these are A. comosus genes). The genes are subdivided in many different categories (see table S1 as supplementary material of Yardeni et al. 2021), and I chose the following to select genes of interest related to CAM:

Differentially expressed in CAM / C3 experiment (186)
Positive selection in CAM / C3 shifts (79)
Gene families associated with CAM (22)
CAM-related Acomosus Ming2015 (29)
stomata function (48)
aquaporin regulation (24)
drought resistance (61)
circadian metabolism (47)
malate transferase (28)
circadian clock (3)

Selection was done with following code:

`awk 'BEGIN{FS="\t"} ($6=="yes")||($7=="yes")||($8=="yes")||($10=="yes")||($11=="yes")||($12=="yes")||($13=="yes")||($16=="yes")||($17=="yes")||($21=="yes") {print $0}' Bromeliad1776_gene_list.csv > Genes_of_interest_CAM_related.csv`

Using the Acomosus IDs of this subset, I searched the multi-copy gene families in our orthology table. For CAM related genes, 76 out of 500 or 15 % of genes were in multicopy families. Using the script `script_make_gene_family_categories.py` I separated these genes into categories based on the gene counts per species. The scripts adds the category to the original subset file and puts out a summary of counts of orthogroups per categories as such (for CAM related genes):

    Types of CNV changes:
    Nr. gene families that remained equal: 19
    Nr. gene families larger in A.comosus: 11
    Nr. gene families larger in Tillandsia: 2
    Nr. gene families larger in T. fasciculata: 8
    Nr. gene families larger in T. leiboldiana: 2
    Nr. gene families with changes in more than one species: 16
