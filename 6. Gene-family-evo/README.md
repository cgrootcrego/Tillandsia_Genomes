# Inferring Gene Family Evolution between *T. fasciculata* and *T. leiboldiana*

Using the inferred gene models of both reference genomes, we can get a first idea of which gene families have changed in size between the two species. However, before looking at the families with greatest changes, some filtering of faulty gene models has to be done first, especially for the *T. fasciculata* annotation.

# Detecting faulty gene models resulting from haplotigs in *T. fasciculata*

*T. fasciculata* is a wide-spread species, often regarded as a species complex, with very high levels of heterozygosity. Unsurprisingly, the accession used for our *T. fasciculata* reference genome was also more heterozygous than desired for genome assembly. This is especially an issue for inferring gene family sizes, since haplotig genes will be considered separate copies.
Using 50x whole genome sequencing data of our reference genome accession, we can look at coverage across the genome, and more importantly across gene models. The idea here is that legitimate duplications should have similar coverage to the rest of the genome. However, erroneous duplications that occurred for example due to the existence of haplotigs, will have a lower coverage than the rest of the genome, as reads from one gene will be distributed across two genic regions in the genome.

**1) Alignment of 50x WGS data to T.fas reference:**

Alignment of the raw data to T.fas genome was done using bowtie2. The choice of aligner is because the manual of bowtie2 explicitly mentions that reads aligning equally well in two different locations will be assigned to one of these two locations at random. This is important, since we want to infer a decrease in coverage when gene duplicates are false, and this would be more difficult if reads are mapped at multiple locations. BWA doesn't explicitly specify this.

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

Per-base coverage in genic regions of T. fasciculata was calculated using samtools depth. I first obtained a bedfile with the coordinates of only our curated gene models (see orthofinder analysis):

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

Using this bedfile, I then calculated per-base coverage with samtools. Coverage of all bases included in the bedfile was reported, including areas of zero coverage (-a):

    samtools depth -a -b $bedfile $bamfile -o $output

I then appended the name of the gene at each position in the samtools depth output file using the python script `script_add_gene_names_to_cov_file.py`.

Note: This script was updated later on to include overlapping genes. See [here](https://github.com/cgrootcrego/Tillandsia-compgenomics/tree/main/Annotation/Gene_curation/overlapping_genes) for more info on these overlapping genes.

Mean and median coverage were then computed per gene and compiled into a table containing also orthology information with the script `script_calculate_mediancov_add_info.py`:

    python script_calculate_mediancov_add_ortho-info.py Tfas_genemodel_assessment/CovDepth_Tfas_ortholog_regions.fromgff.edited.txt Tfas_genemodel_assessment/curated_Tfas_orthologues.txt > computed_lengths_cov.txt

The file `computed_lengths_cov.txt` records the length of the vector of coverage entries for each gene. This was a way for me to make sure the script was recording the full length of overlapping genes.

# Correcting gene family sizes based on average coverage

I then made density plots of the average coverage per gene for different categories of genes using the Rscript `Assessing_multicopy_genemodels_cov.R`. These density plots showed that multi-copy gene models in *T.fasciculata* have a bimodal distribution of mean coverage - meaning that about half of the multicopy genes have a mean coverage around the genome-wide average, while the other half has a much lower mean coverage. I then determined a coverage threshold using a finite mixture model (FMM) with the package [cutoff].(http://marcchoisy.free.fr/fmm/index.html). This cutoff represents an estimated "split" of the two peaks in the bimodal distribution, in other words, the point at which, if we would split the bimodal distribution into two unimodal ones, a datapoint is unlikely to belong to both distributions. So, I used this threshold to separate "true" gene models from "faulty" ones. The threshold was determined at a median and mean coverage of 35. I then isolated all multicopy genes with a mean coverage under the threshold.

Note: the average full genome coverage was calculated by running samtools depth without specifying regions and running the awk one-liner `awk '{ total += $3; count++ } END { print total/count }`. To obtain the median coverage over the full genome, I ran `cut -f 3 CovDepth_perbase_fullgenome.txt | sort -n | awk -f median.awk`. Median.awk was the following short script:

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

I then decided to rescale gene family sizes for the orthogroups which contained faulty multi-copy genes. There were 1981 such orthogroups, out of 17,641 (11 %), containing 6861 genes (26 %) I corrected gene family sizes by obtaining the total mean coverage in the orthogroup (sum of mean coverage of all Tfas genes in the OG) and dividing this by the "expected mean coverage" (full genome mean coverage x number of Tfas genes in the group). I called this ratio the correction factor, which I multiplied by the original number of   Tfas genes in the orthogroup to obtain a corrected family size.

For 63 orthogroups, the correction factor was > 1, meaning that the total average coverage of the orthogroup was larger than expected. In these cases, coverage is non-informative regarding the validity of genome sizes, so I decided not to correct for these families.

for 803 orthogroups, the total mean coverage of the orthogroup didn't even reach 46 (what we would expect for 1 gene). In this case, the correction usually brought the gene family size down to 1, but sometimes also to 0. In those cases, since the gene has been properly annotated and has an orthologous sequence in at least *T. leiboldiana* or *A.comosus*, I assumed the gene must exist and therefore corrected the size back up to 1.

I ran the same corrections for T. leiboldiana, although I couldn't separate faulty genes due to the fact the distribution was unimodal with just a small shoulder. However, while in the case of T.fas 1719 orthogroups were corrected, only 424 orthogroups were size corrected in T. leiboldiana.

Next, I integrated the new family sizes into my general per-gene orthology table, which we will feed back into gene family evolution analysis, with the python script `script_insert_corrected_sizes.py`.

# Changes to workflow July 27th

I made a few changes to the workflow above:
 - orthogroups containing mito / plastid / ribosomal information were removed just before size correction, as it probably doesn't make sense to correct sizes for these genes since their coverage doesn't exactly reflect their presence in the genome (multiple plastids / mitochondria in a cell, ribosomal genes have very complicated gene family evolution)
 - Size corrections were made for all multi-copy orthogroups in Tfas, instead of only on the orthogroups with genes < 34.5x. This only expanded the number of orthogroups run through size corrections from 1981 to 2632 (25 % more) and from 6861 genes to 8400 genes (19 % more).
 - Size corrections were also made for upward corrections in orthogroups that did not belong to chloroplasts.
 - Size corrections were made not based on the full genome average but on the average coverage of ancestral single-copy genes. The idea is that this category represents better the sort of coverage expected in genic regions (full genome coverage is extremely variable).

Testing different forms of size corrections:
 - Accounting for "expected" coverage variability: I designed a size correction method where orthogroups with coverage lying within an interval from the expected coverage were not corrected, under the assumption that this is "expected" coverage deviation from the mean. However, it was hard to find proper thresholds for this interval since the coverage of ancestral-single-copy genes has a relatively wide distribution and taking strict cutoff such as the 1st and 99th quantile would lead to a much too large interval. I tested it by taking the 25th and 27th interval, but it is difficult to argue that this interval really accounts for "expected" coverage variability, it is quite arbitrary and it also didn't impact many orthogroups. So we decided to drop this avenue.
 - Designing a per-gene size correction strategy: I decided to generate correction factors on a per-gene basis, by taking the average ancestral-single-copy average as a correction factor 1 and scaling all orther coverages against it. Then I would sum up the weighted gene values to obtain a new size. This is in fact the same as doing a per-orthogroup approach. So we can continue with the per-orthogroup approach as I find it more intuitive.

 After final size corrections, the total gene counts for Tillandsia fasciculata went down from 8139 to 6763 (1376 or 17 % difference). For Tillandsia leiboldiana, gene counts went from 4098 to 4026 (72 or 2 % difference).

# Studying gene family size differences in Tfas and Tlei

I first selected only orthogroups that were not mitochondrial, ribosomal or plastid. The procedure for this combined both gene name searches and blasting to existing plastid and mitochondrial sequences. The details are described in `Removing_plastid_genes.md`

This resulted in 69,229 genes and 18,697 orthogroups.

I then obtained a table of per-orhogroup gene counts by selecting the orthogroup and count fields of the orthology table:

    cut -f 7,8,9,10 orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo2.txt | sort -u > orthogroups_Tfas_Tlei_Acom.counts.no_TEs.size_corrections.no_plastid-mito-ribo2.txt

The relationship of gene family size between species was explored in `Gene_family_evolution_GO_term_enrichment.R`

To better understand which levels of orthogroup size changes between species were "substantial" or "deviating from expectations", I calculated Log ratios of the family sizes and ranked them. This showed that by selecting the top 2 % of observations, already the full multicopy gene set was included. This is mostly because we have a very large set of orthogroups where family sizes have not changed between both species. My interpretation of this is that we have no real power to distinguish between multicopy families which are fast or slow evolving, and that in any case any form of family size changes is unexpected. Therefore, it may be justified to analyze all multicopy groups as done below.

For downstream multicopy analyses (GO terms and genes of interest) I made a selection of genes that were not unique to one species and where at least one of the Tillandsia species had a gene count > 1.This set contained 2061 orthogroups.

I still made an additional selection of fast-changing families by selecting within the multicopy dataset the top 2 % of orthogroups by log ratio, in both directions. This resulted in 53 orthogroups where family size is larger in T. fasciculata and 45 orthogroups where size is larger in T. leiboldiana. This translated to 1182 genes for the former selection (1297 after correction, 9.7 % difference, 755 genes have no GO terms, 64 %) and 480 genes for the latter selection (501 after correction, 4.3 % difference, 295 genes have no GO terms, 59 %). GO enrichment was done as explained below.

# GO term enrichment of multi-copy gene families

To have an idea of what sort of gene families occur in multicopy in either species, I used the above set of orthogroups,  where at least one species has a gene count > 2 and no species-unique groups are included. Using this list of orthogroups (obtained from the R script) I extracted the underlying genes (uncorrected sizes):

    grep -f orthogroup_selection_multicopy_for_GO_term_ID.txt \
    orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch.txt > \
    dup_genes.txt

Then I removed all pineapple genes (`grep -v "Aco" dup_genes.txt`) to obtain the list of genes for GO term enrichment, which was done using TopGo.

The subset of multicopy genes contained 6642 genes in T. fasciculata and 4820 genes in T. leiboldiana. After size correction, the sum of reported sizes for these two subsets are 6261 genes and 4693 genes. In other words, the difference in gene numbers is not very large (594 in Tfas and 127 in Tlei)- perhaps a correction in GO terms is therefore not necessary for the big scale. Of the selected multi-copy genes, 1968 (29 %) in Tfas and 1169 (24 %) have no GO terms at all.

GO term enrichment was run on all multicopy genes (2061 orthogroups), on a subset of genes belonging to families where copy number was larger in Tfas (916 orthogroups, 5721 genes) and on a subset of genes belonging to families where copy number was larger in Tlei (583 orthogroups, 2975 genes). This can give us an idea of the direction in which certain functions have had gene family expansion.

Important: somehow some neuronal go terms have ended up in our annotation, we need to clean this out??

GO terms were studied by eye and interesting terms were studied on a gene-level. A number of genes putatively involved in CAM photosynthesis belong to multicopy gene families. However, this was only a small subset of multicopy genes and broad patterns of gene duplication related to the evolution of key innovative traits remained unclear in both species. This is mostly because GO terms poorly reflect the traits we are interested in and remain quite difficult to interpret.

The top 100 GO terms in each category (BP, MF and CC) were stored alongside with their p-value and underlying orthogroups and genes, for further analysis in combination with RNA-seq data.

# Detecting "genes of interest" in the multicopy gene family pool

A reverse approach, where we search for genes of interest in our set of multicopy gene families, may give us a better idea on gene family evolution and key innovation traits. I made use of the ~ 700 genes of interest identified in previous work by Marylaure and Gil to see if a proportion of these genes were to be found in our mutlicopy gene family set.

Genes were selected from the Bromeliad1776 kit, which relies on the pineapple annotation (so these are A. comosus genes). The genes are subdivided in many different categories (see table S1 as supplementary material of Gil's target capture paper), and I chose the following to select genes of interest related to CAM:

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

This resulted in 500 genes.

A similar selection was made for flower-related genes:
Anthocyanin and self incompatibility
Flavonoid and anthocyanin pathways in Pitcairnia

`awk 'BEGIN{FS="\t"} ($15=="yes")||!($18=="NA") {print $0}' Bromeliad1776_gene_list.csv > Genes_of_interest_flower_colour_related.csv`

This resulted in 57 genes.
Using the Acomosus IDs of this subset, I searched the multicopy gene set. For CAM related genes, 76 out of 500 or 15 % of genes were in multicopy families. For flowering genes, 6 out of 57 (11 %) genes are in multicopy families. Using the script `script_make_gene_family_categories.py` I separated these genes into categories based on the gene counts per species. The scripts adds the category to the original subset file and puts out a summary of counts of orthogroups per categories as such (for CAM related genes):

    Types of CNV changes:
    Nr. gene families that remained equal: 19
    Nr. gene families larger in A.comosus: 11
    Nr. gene families larger in Tillandsia: 2
    Nr. gene families larger in T. fasciculata: 8
    Nr. gene families larger in T. leiboldiana: 2
    Nr. gene families with changes in more than one species: 16
