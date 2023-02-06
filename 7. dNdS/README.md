# Inferring signatures of selection with dN/dS ratios between *T. fasciculata* and *T. leiboldiana*

In this section, we infer dN / dS ratios for all single-copy pairs of orthologous genes between *T. fasciculata* and *T. leiboldiana*. We obtained dN/dS ratios by aligning the orthologous sequences using a codon-aware aligner ([MACSE](https://github.com/ranwez/MACSE_V2_PIPELINES)) and inferring ratios with [PAML](http://web.mit.edu/6.891/www/lab/paml.html).

## Selecting single-copy orthogroups and compiling per-orthogroup fasta files

We compiled the gene sequences of each one-to-one orthogroup into individual fasta files. We selected all orthogroups with a single-copy gene in *T. fasciculata* and *T. leiboldiana* (regardless of the number of copies in *A. comosus*).

To extract the features belonging only to one-one orthologous genes (here only shown for *T. fasciculata* genes):

	# Select all features of single-copy genes:
	grep -f one-to-one_orthogroups_Tfas_Tlei.IDonly.txt \ Tillandsia_fasciculata_v1.2.edited_allfeatures.gff > \
	Tillandsia_fasciculata_v1.2.edited_allfeatures.one-to-one_orthologs.gff
	# Retain only "CDS" features and modify the 9th field to allow the next step to work
    	awk '$3 == "CDS" {print $0}' [gff] | sed 's/ID=/Name=/g' > [gff.CDSonly]

Fasta sequences for each ortholog was extracted using with `01_script_CutSeq_modified.py`:

    while read line; do python2 cutSeqGff_mod.py  \
    Tfas_assembly/per_chrom_fasta/Tfas_$line.fasta \
    Tillandsia_fasciculata_v1.2.one-to-one_orthologs.CDSonly.gff \
    $line CDS; done < chr_Tfas.txt

Note: this script needs per-scaffold fastafiles of the genome assembly to work efficiently, and a list of these scaffold names (chr_Tfas.txt).

A list of orthologous gene pairs per orthogroup was compiled with the python script `02_script_compile_per-orthogroup_gene_list.py`. This list wass used to join the fasta sequences originating from the two genome assemblies in single per-orthogroup fasta files. The following bash script was run to compile orthologous fasta sequences into one fasta file:

    cat one-to-one_orthogroups_Tfas_Tlei.perOG.txt | while read line ;
    do
     Orthogroup=`echo "$line"|awk '{print $1}'`;
     echo $Name;
     Tfas=`echo "$line"|awk '{print $2}'`;
     Tlei=`echo "$line"|awk '{print $3}'`;
     cat ../Tfas_seq/$Tfas:cds.fst ../Tlei_seq/$Tlei:cds.fst > $Orthogroup.fasta
    done

## Checking for completeness and length differences between orthologous pairs

Since dN/dS ratios are very sensitive to alignment quality, we first wanted to have an idea of the completeness of our gene pairs (i.e. presence of a start and stop codon) and their relative lengths. I therefore made a "checklist" containing this information for each gene using the script `03_script_make_checklist_completeness_length_orthologous_genes.py`. We also collected exon numbers per gene with `04_run_count_exon_per_gene.sh`.
The checklist was modified to a per-orthogroup format and the absolute difference in sequence length between pairs was also calculated. We calculated relative differences for each species (diff/length_Tfas and diff/length_Tlei) and noted completeness in the orthogroup (both complete / one complete / zero complete).

## Testing filtering options of alignments based on length differences and completeness

With the info obtained above, we filtered orthogroups based on completeness and length using `05_completeness_filtering.R` and `06_codeforfiltering.txt`.We created four subsets with different degrees of filtering:

- stringent filtering: both alignments are complete (they have a start and stop codon), they have an absolute length difference of maximum 200 basepairs and they have a relative length difference that is smaller than 0.1 (meaning that the difference in length is no more than 10 % of the total length of that gene).
- relaxed completeness: same requirements as above regarding length differences, but only one gene is complete.
- relaxed relative difference: both genes have to be complete and the same requirement applies for absolute length differences as above. However, the relative length difference ranges from 0.1 to 0.2.
- relaxed completeness and relative difference: absolute length requirement remains the same, but relative differences range from 0.1 to 0.2 and only one gene is complete.

The total number of orthogroups qualifying for each filtering level was:
6208  stringent filtering (47 %)
493   relaxed relative difference (4 %)
2569  relaxed completeness (19 %)
316   relaxed relative difference and completeness (2 %)

We then randomly selected 200 orthogroups for each filtering setting in `06_codeforfiltering.txt` and performed alignments:

    # Extensive alignment
    for i in *; do
     java -jar macse_v2.05.jar -prog alignSequences -seq $i -local_realign_init 1 -local_realign_dec 1 ;
    done

The alignments were then converted into .axt files using a modified form of the script
03_ConvertFasta.py from the {AlignmentProcessor}(https://github.com/WilsonSayresLab/AlignmentProcessor/blob/master/bin/03_ConvertFasta.py) software. Preliminary dN/dS ratios were then calculated with [KaKs calculator2](https://github.com/kullrich/kakscalculator2) and the MYN model using a modified script of 04_CallKaKs.py from ALignmentProcessor.

Distribution of dN/dS ratios across filtering subsets were assessed in `07_Assess_filtering_alignments.R`. While relaxing length difference didn't change dN/dS ranges, relaxed completeness did cause a few outliers to show up, though most estimates were still in the same range as with stringent filtering. We therefore decided to keep the filtering loose at this point, and only remove alignments of genes where both are incomplete and where length differences are > 0.2 in both genes. This kept 10,362 orthologous pairs for pairwise alignment, which is almost 80 % of the original set.

Using the checklist, we selected all genes meeting quality requirements:

	grep -v "0 complete" checklist_with_diff_completeness.txt | \
	awk '$5 > -0.2 && $5 < 0.2 && $6 > -0.2 && $6 < 0.2' | cut -f 1 | \
	sed 's/"//g' > orthogroups_final_subset_0.2_lengthdiff_1_complete.txt

## Calculating pairwise dN / dS ratios

Alignment was done with optimization parameters for all orthologous pairs with the bash script `08_run_macse.sh`. The selected alignments were converted to PHYLIP format using `09_03_ConvertFastatoAXTorPhylip_modified.py`. I used modified scripts from ALignmentProcessor `10_04_CallCodeML_modified.py` and `10a_parallelcodeML.py` to run codeML automatically for all orthologous pairs.

Note: Phylip files were slightly modified because codeML doesn't accept "!" in the alignments. MACSE introduces these whenever there is a change in frameshift. In other words, when an entire codon is deleted, this will be shown as "---" but when there is a gap < 3, it will show as "!!A" or something of the like. I replaced all "!" by "-" so that codeML wouldn't throw errors.

codeML was run twice for a null model and alternative hypothesis. The settings are summarized in the two control files `11_codeml_null.ctl` and `12_codeml_alt.ctl`. Concretely, we used the site-specific model MO (Nssties = 0) for both runs. In the null model, omega (dN/dS) was fixed to 1, whereas in the alternative, it was estimated from the data starting from 1. The reason to run codeML twice is that we gather likelihoods for the null and alternative model and can compute p-values, to see how significant our inferences are.
Additional choices in the settings of codeML were to set codon frequencies to be inferred from the data at all three positions in the codon (F3x4) (CodonFreq = 2) and kappa (Ts/Tv ratio) was set to 3, but is estimated from data (kappa = 3, fix_kappa = 0). The results of both runs were compiled with the script `13_script_compile_codeml_LRT.py`. This script also performs a likelihood ratio test (chisquare). I called the script on all files as follows:

    null=(codeML/null/*.out)
    alt=(codeML/alt/*.out)
    n=0
    for ((i=0;i<=${#null[@]};i++))
    do
      python script_compile_codeml_LRT.py "${null[i]}" "${alt[i]}"
      n=$((n+1))
      echo $n
    done

Mulitple testing correction was performed with `14_script_multiple_testing_correction.py`.+

The dN/dS results were analyzed with `15_Assessment_dnds_values.R`. Alignments of significant genes were checked with AliView. Additionally, RNA-seq data used for annotation was mapped back to the main scaffolds and visualized in Jbrowse to further assess alignments.

# Obtaining per-chromosome dN/dS statistics

To study the general dN/dS distribution across scaffolds, we decided to only work with alignments that had at least 5 variant sites. We obtained this information with [AMAS](https://github.com/marekborowiec/AMAS) and created a list of orthogroups to keep. The dN/dS distribution across all scaffolds and dedicated study of rearranged scaffolds was performed in `16_Script_dNdS_perchrom.R`.

# Running dN/dS for 1:1:2 and 1:2:1 paralogs

Pairwise tests were also run for duplicated orthogroups that are either 1:1:2 (duplicated *T. leiboldiana*) or 1:2:1 (duplicated in *T. fasciculata*). We decided to be very stringent with this analysis and only work with orthogroups that did not change after size correction (the group is reported as 1:1:2 or 1:2:1 both before and after size correction). This resulted in 108 groups in 1:1:2 and 190 groups in 1:2:1.

Obtaining fasta sequences per orthogroup was done similarly as above but with `17_script_compile_per-orthogroup_gene_list.paralog.py`. Concatenation of the fastasequences into two files happened with following loop:

	cat ../orthologs-121.perOG.txt | while read line ;
	do
 		Orthogroup=`echo "$line"|awk '{print $1}'`;
 		echo $Orthogroup;
 		Tfas1=`echo "$line"|awk '{print $2}'`;
 		Tfas2=`echo "$line"|awk '{print $3}'`;
 		Tlei=`echo "$line"|awk '{print $4}'`;
 		cat ../Tfas_seq/$Tfas1:cds.fst ../Tlei_seq/$Tlei:cds.fst > ${Orthogroup}_Tfas-copy1.fasta
 		cat ../Tfas_seq/$Tfas2:cds.fst ../Tlei_seq/$Tlei:cds.fst > ${Orthogroup}_Tfas-copy2.fasta
	done

MACSE alignments were performed with `18_run_macse_paralogs121.sh` and `19_run_macse_paralogs112.sh`. Alignments were converted to phylip format with the above modified python script:

    python 03_ConvertFastatoAXTorphylip_modfied.py -i alignments_121/ -o phylip_121/ --phylip  

Tests were run between each duplicate copy of the species with paralogs and the single-copy ortholog of the species with a single-copy gene with a null and alt model in codeML. I compiled the results and performed LRT with `20_script_compile_codeml_LRT.paralogs112.py`.
