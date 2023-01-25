## Gene annotation in Tillandsia

# Overview

Gene annotation steps:
1. Prepare input for braker: splice information and busco protein sequences
2. Run braker on the assembly using splice-info from mapping and busco proteins
3. Prepare input for maker:
4. Run maker 1 / 2 times
5. Calculate number of genes, average length, busco score on transcripts, TE content, map RNAseq back to gene models, do feature counts ...
6. If necessary, make subset based on TE content
7. Rename maker names in gff
8. run functional annotation in Blast2Go, compare blast2go runs with different DB


# 1. Preparing files for Braker

The input for braker is a soft-masked version of the assembly (See 2a. TE annotation), the splice information from RNA-seq mapping and single copy proteins predicted during BUSCO.

## 1.1. Splice information extraction:

Map RNA-seq to assembly:

	# Genome indexing
	/apps/star/2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
	--genomeDir [genome_directory] \
	--genomeFastaFiles [hardmasked_genome] \
	--runThreadN 48
	# Mapping
	/apps/star/2.5.3a/bin/Linux_x86_64/STAR --genomeDir [genome_directory] \
    --readFilesIn [RNA_fastq_R1] [RNA_fastq_R2] \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix rna_to_assembly_fasc_knTEs \
    --limitBAMsortRAM 10240000000 --runThreadN 48 \

Extract splice info:

	/apps/augustus/3.3.3/bin/bam2hints --in [File.Aligned.sortedByCoord.out.bam] --out [splice_info_from_rna.gff]

## 1.2. BUSCO protein sequences:

  busco -i [hardmasked_genome] -o [output_name] -l liliopsida -m genome --long --augustus_parameters='--progress=true' -c 24

The protein files in busco/sequences/single_copy_sequences (.faa) were concatenated with cat:

 	cat *.faa > single_copy_protein_hints_busco.faa


# 2. Running Braker

Braker was run with the following command (example of *T. leiboldiana*):

	source /home/fs71400/grootcrego/.bashrc
	conda activate braker-env
	braker.pl --species=[Tleiboldiana] --genome=tillandsia_leiboldiana_assembly.Softmasked.fasta --softmasking --hints=splice_info_from_rna_leiboldiana.gff --prot_seq=single_copy_protein_hints_busco_leiboldiana.faa --prg=gth --gth2traingenes --cores=32 --gff3

BUSCO was run on the augustus.hints.codingseq fastafile with following command:

	busco -i augustus.hints.codingseq -o busco_Tleiboldiana_braker_pred -l /home/fs71400/grootcrego/busco_downloads/lineages/liliopsida_odb10 -m transcriptome -c 16

# 3. Preparing files for Maker

We used the following input files for Maker:
- EST library, masked for complex repeats
- gff with complex repeats from EDTA
- Braker gene models
- Ananas and monocot curated protein sequences
- TE library from EDTA + Giri Monocot DB


## 3.1 EST library

Transcriptome assembly was done previously with Trinity (see Notebook entry "transcriptome assembly" from 06.02.2020). Yet, we have learned from T. fasciculata that these may contain active TE sequences, and therefore we will pass it through repeatmasker once. We only mask complex repeats, as simple repeats can be part of genic regions. The command used to generate a masked assembly:
  RepeatMasker -pa 4 -e ncbi -dir . -gff -nolow -lib /proj/grootcrego/Genome_assemblies/leiboldiana/3_masked_assembly/TE_families_EDTA_T.leiboldiana_KnownOnly.fasta leiboldiana_transcriptome_assembly.fasta

  - Complex repeat GFF:
By providing maker a list of complex repeats in the genome and their positions, maker can hardmask these areas to prevent TEs from being annotated as genes. To generate this list, we used the repeatmasker util rmOuttoGFF3 to create a gff file based on the masked assembly with only known TEs (EDTA gff file contains all TEs):
    /apps/repeatmasker/4.0.7/util/rmOutToGFF3.pl /proj/grootcrego/Genome_assemblies/leiboldiana/3_masked_assembly/Hardmask_T.leiboldiana_EDTA_knownOnly/tillandsia_leiboldiana_assembly.pilon.upper.sorted.fasta.out > /scratch/grootcrego/maker/leiboldiana/tillandsia_leiboldiana_EDTA_knownTEsonly.fasta.mask.gff3
Important: Because we only ran repeatmasker for complex repeats (--nolow), there's no need to filter out simple repeats in this gff, as we did in T.fasciculata. Technically, this step doesn't have to be performed as repeatmasker generated a gff file. But since I see slight differences in the output, I decided to still perform this step for the sake of consistency between T.fasciculata and T. leiboldiana.
Afterwards, we need to reformat the gff3 to make it compatible with maker:
  cat tillandsia_leiboldiana_EDTA_knownTEsonly.fasta.mask.gff3 | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > tillandsia_leiboldiana_EDTA_knownTEsonly.fasta.mask.complex.reformat.gff3

  - Protein database:
To the already existing file containing Ananas and monocot swissprot proteins, I added the protein sequences of all T.fasciculata genes that had a blast result (31 574 sequences). The total file had 96 322 sequences and can be found in the maker folder on the vsc4.

  - TE database:
I fed into maker the file TE_families_EDTA_T.leiboldiana_KnownOnly_AND_Giri_monocot_TE-DB.fasta.


2.4. Running maker
---
Maker was run on the VSC4 with two jobs, using all 4 available private nodes in mpi mode. Options can be found in /gpfs/data/fs71400/grootcrego/maker_Tlei/maker_R1_Tlei_opts.ctl and the running command was:
  mpirun -np 48 /home/fs71400/grootcrego/software/maker/bin/maker -base leiboldiana_R1 maker_R1_Tlei_opts.ctl /home/fs71400/grootcrego/software/maker/maker_bopts.ctl /home/fs71400/grootcrego/software/maker/maker_exe.ctl

The run took about 5-6 days. I had a couple of issues which delayed the run, for example, I first ran maker in a directory that only has 30 GB, which caused maker to crash. I then moved the entire directory to another partition with sufficient memory, but when restarting the run, maker somehow didn't recognize about 2500 scaffolds. I fixed this by running maker with the -dsindex option, which simply restructures the datalog of that specific run (add the correct basename). This then included all 10433 scaffolds. Of these, about 200 had been tried 2 times in the previous run and were therefore not being restarted. I eliminated their entries in the datalog files and also increased the number of tries in the opts file. This then worked for all but one scaffold. No matter how many times I retried it, it didn't work. Finally I resolved to eliminate it's folder in the datastore directory, removing all entries in the datalog file and rerunning it.


2.5. Examining maker gene models
---
I obtained the final maker gff file and fasta sequences with following commands:
  ~/software/maker/bin/gff3_merge -n -s -d leiboldiana_R1_master_datastore_index.log > leiboldiana_R1.all.maker.gff
  ~/software/maker/bin/fasta_merge -d leiboldiana_R1_master_datastore_index.log

To count the number of gene models and average length I ran:
  cat leiboldiana_R1.all.maker.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'

I also ran BUSCO on the transcript sequences with following command:
  busco -i leiboldiana_R1.all.maker.transcripts.fasta -o busco_maker_R1_Tleiboldiana -l /home/fs71400/grootcrego/busco_downloads/lineages/liliopsida_odb10 -m transcriptome -c 24

And I looked at the percentage of gene models with an AED score < 0.5:
perl ~/software/maker/bin/AED_generator.perl -b 0.025 leiboldiana_R1.all.maker.gff

I mapped raw RNAseq data to the transcripts too, using STAR. For this, it is important to figure out the total number of bases in the fastafile, as STAR needs an extra option for "small" genomes when indexing. To count the total number of bases:
  grep -v ">" leiboldiana_R1.all.maker.transcripts.fasta | wc | awk '{print $3-$1}'
This gave 62,536,087 bp for Round 1. The default --genomeSAindexNbases is 14, but has to be scaled down for smaller genomes. The calculation is as follows: log2(genome-length)/2 - 1. With the above mentioned length this is: 11.9 (12). The STAR command used was:
  /apps/star/2.7.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir /scratch/grootcrego/maker_tlei --genomeFastaFiles /scratch/grootcrego/maker_tlei/leiboldiana_R1.all.maker.transcripts.fasta --runThreadN 24
  /apps/star/2.7.3a/bin/Linux_x86_64/STAR --genomeDir /scratch/grootcrego/maker_tlei --readFilesIn /proj/grootcrego/RNAseq/leiboldiana/all/leiboldiana_RNA_all_trimmed.p1.fq /proj/grootcrego/RNAseq/leiboldiana/all/leiboldiana_RNA_all_trimmed.p2.fq --outSAMtype BAM Unsorted --outFileNamePrefix rna_to_maker_transcripts_Tlei_R1. --limitBAMsortRAM 10240000000 --runThreadN 24

The same STAR command was run on full gene models (exons and introns). To obtain these sequences, I used bedtools:
  awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' leiboldiana_R1.all.maker.gff | bedtools sort | bedtools merge -d 100 > leiboldiana_R1.all.maker.bed
  bedtools getfasta -fi ~/Tlei_assembly/tillandsia_leiboldiana_assembly.pilon.upper.sorted.fasta -fo  leiboldiana_R1.all.maker.full_genemodels.fasta -bed leiboldiana_R1.all.maker.bed
Take note: I changed the --genomeSAindexNbases from 12 to 13 due to size differences.

I decided to run feature counts to see genes with multimapping RNA-seq (may be TEs that I could filter out). The following command was used:
  featureCounts -a leiboldiana_R1.all.maker.gff -o ../counts_Tlei_RNAseq.txt -g ID -t gene -G /proj/grootcrego/Genome_assemblies/leiboldiana/2_polished_assembly/tillandsia_leiboldiana_assembly.pilon.upper.sorted.fasta -T 8 -p -s 2 -M --tmpDir tmp/ /scratch/grootcrego/STAR/Tlei_mapping_to_unmasked_assembly/leiboldiana_RNA_to_unmasked_assembly_trimmed.Aligned.sortedByCoord.out.bam
Note the -M option which includes also multi-read counting.

I also investigated the extent of genes lying in masked areas. To do this I followed Thibaults guidelines set out in "Filtering and selecting gene models" from 02.04.2020. Genes with high TE content can be filtered out when training maker in a second round. This may improve the annotation.

A summary of this examination can be found on Google drive in the sheet "gene annotation procedure".


2.6. Training predictors for second round
---
In the second round of maker, the predictors are further trained on the genemodels estimated in round 1. We know from our T.fasciculata annotation rounds that maker's second round usually yields just as good or even worse gene models as the first round. Our first round of T. leiboldiana gene prediction yielded more genes than one would expect (> 50 thousand, when 30-35 thousand would be expected). Most of these genes lie completely outside masked content (around 40 thousand), but the remaining ones could be TEs which cause an overestimation of gene models. In this second round, I want to make sure the gene predictors are being trained by genes which are 100 % outside masked regions, additional with the usual prerequisites such as a length > 50 amino acids and an AED-score under 0.25 (meaning they have good support from raw evidence).

2.6.1. Training SNAP
------
Training was done in directory /gpfs/data/fs71400/grootcrego/maker_Tlei/snap/ in the vsc4.
9252
8783
To only allow training with genes containing 0% masked content, I isolated the names of these genemodels from the TE content summary produced in 2.5:
  awk '$10 == 0 {print $1}' summary_TE_GCcontent_Leiboldiana.txt > genemodels_no_te_content.txt
I removed all the mRNA-1, etc., endings to these names to select the gene feature. I added the line "contig" so that scaffold entries would also be included in the subset when grepping. Then, I created a subset of the gff file containing only these genes:
  cat leiboldiana_R1.all.maker.gff | parallel -j 4 --pipe --block 10M grep -f genemodels_no_te_content.txt > leiboldiana_R1.all.maker.only_fully_unmasked_genes.gff
Using this newly produced gff file, I created the zff file for snap:
Selecting genemodels with AED < 0.25 and > 50 amino acids:
  ~/software/maker/bin/maker2zff -x 0.25 -l 50 ../../leiboldiana_R1.maker.output/leiboldiana_R1.all.maker.only_fully_unmasked_genes.gff
This resulted in errors which also occurred when running the full gff file from maker. In other words, something seems to be wrong with this command when using gff files. Because I can't filter out genes with the datastore log, I decided not to train SNAP.


2.6.2. Training Augustus
------
Obtaining training sequences for Augustus is more straightforward and allows us to filter out genes with TE content:
First, I selected the genes which had no TE content:
  awk '$10 == 0 {print $1}' summary_TE_GCcontent_Leiboldiana.txt > mRNA_no_te_content.txt (40330 mRNAs)
Then, I made a subset gff containing only these transcripts:
  cat leiboldiana_R1.all.maker.gff | parallel -j 4 --pipe --block 10M grep -f mRNA_no_te_content.txt > leiboldiana_R1.all.maker.only_fully_unmasked_mRNAs.gff
From this gff, I printed all mRNA entries and their start and end point. Then, I added 1000 bp on both ends to get the surrounding area. Then, using bedtools, I obtained the fasta sequences for these regions:
  awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' leiboldiana_R1.all.maker.only_fully_unmasked_mRNAs.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi ../../Tlei_assembly/tillandsia_leiboldiana_assembly.pilon.upper.sorted.fasta -bed - -fo Bcon_rnd1.all.maker.transcripts1000.fasta
This resulted in 40026 genes.

I then ran busco (--long version) on these sequences to train augustus. I used maize as a species:
  busco -i leiboldiana_R1.all.maker.transcripts1000.fasta -o busco_for training_Augustus_R1_Tleiboldiana -l /home/fs71400/grootcrego/busco_downloads/lineages/liliopsida_odb10 -m genome --long -sp maize -c 24

The busco results of this subset on liliopsida was:
  C:84.5%[S:79.7%,D:4.8%],F:7.4%,M:8.1%,n:3236
Which is quite close to the busco score of the full genemodel set (85 %, 6.2 % duplicated). This suggests that our training subset of 40 000 genes is almost as complete in real genes as the full set.

Inside the busco output folder, the retrained parameters can be found in:
  /gpfs/data/fs71400/grootcrego/maker_Tlei/leiboldiana_R1.maker.output/busco_for_training_Augustus_R1_Tleiboldiana/run_liliopsida_odb10/augustus_output/retraining_parameters/BUSCO_busco_for_training_Augustus_R1_Tleiboldiana

Here, I renamed the files to Tillandsia_leiboldiana:
rename BUSCO_busco_for_training_Augustus_R1_Tleiboldiana Tillandsia_leiboldiana_R1 *
And modified two of the configuration files:
sed -i 's/BUSCO_busco_for_training_Augustus_R1_Tleiboldiana/Tillandsia_leiboldiana_R1/g' Tillandsia_leiboldiana_R1_parameters.cfg
sed -i 's/BUSCO_busco_for_training_Augustus_R1_Tleiboldiana/Tillandsia_leiboldiana_R1/g' Tillandsia_leiboldiana_R1_parameters.cfg.orig1

Then I created a species folder in the Augustus configuration and moved all these files to that folder:
mkdir $AUGUSTUS_CONFIG_PATH/species/Tillandsia_leiboldiana_R1
cp Tillandsia_leiboldiana_R1*  $AUGUSTUS_CONFIG_PATH/species/Tillandsia_leiboldiana_R1/

Keep in mind that we now have 2 species folders for leiboldiana in our Augustus configuration: Tleiboldiana/ which contains parameters from training in Braker, and Tillandsia_leiboldiana_R1 which contains parameters from training in maker round 1.


2.7. Running Maker round 2
---
To keep matters simple I didn't create separate gff files for our predictions of round1, as I had done previously, but simply referred to the maker gff file within the ctl file:
maker_gff=/gpfs/data/fs71400/grootcrego/maker_Tlei/leiboldiana_R1.maker.output/leiboldiana_R1.all.maker.gff

After round2 finished running I computed the same summary statistics as above. This showed that the genemodels had actually deteriorated (see google sheet).


2.8. Rerunning Maker round 1 without Tfas protein models and extra filtering of Acomosus models
---
Next, I decided to rerun maker from scratch with exactly the same input as resulted in the best genemodels in T. fasciculata. In other words, I removed all Tfas protein sequences. Additionally, I filtered out everything in the protein file which was a transposon by doing:
grep -v "transpos" [protein_file] > Acomosus_F135_CB5_swissprot_monocot_proteins.fasta
This resulted in a file with 64,519 sequences. Apart from this nothing was changed to the maker round.
This resulted in much fewer genes (around 37,000) and a slightly improved BUSCO score. I decided to run a second round using snap and augustus despite having strong suspicions it would turn out worse.

Round 2 was done in a standard way, without filtering out potential TE candidates. Augustus was trained with 37933 sequences, resulting in following results: C:86.6%[S:81.3%,D:5.3%],F:7.1%,M:6.3%. For SNAP I decided to be more stringent with the AED score and only take models with a score < 0.1 (in stead of 0.25). This selected 5127 genes.

Finally, the statistics for round 2 were worse than for round 1 (see google sheet),
meaning that we will proceed to functional annotation with the alternative round 1 genemodels.

2.9. Renaming maker genemodels
---
As was done previously for T.fasciculata gene models, we renamed our final set using maker tools:
# Create a "name map":
~/software/maker/bin/maker_map_ids --prefix Tlei_v1. --justify 5  leiboldiana_R1_alt.all.maker.gff > leiboldiana.maker.name.map
# Change names in GFF:
~/software/maker/bin/map_gff_ids leiboldiana.maker.name.map leiboldiana_R1_alt.all.maker.gff
mv leiboldiana_R1_alt.all.maker.gff Tillandsia_leiboldiana_v1.gff

#For renaming fasta files:
awk '$3=="mRNA" {print $0}' fasciculata.all.maker.gff | cut -f 9 | sed "s/;/\t/g" | cut -f 1,3 | sed 's/=/\t/g' | cut -f 2,4 > mrna_names
cat mrna_names | while read line ;do   Name=`echo "$line"|awk '{print $1}'`;   echo $Name;   Replace=`echo "$line"|awk '{print $2}'`;   echo $Replace;   sed -i "s/$Name/$Replace/g" leiboldiana.maker.name.map2; done
map_fasta_ids leiboldiana.maker.name.map2 leiboldiana.maker.transcripts.fasta
map_fasta_ids leiboldiana.maker.name.map2 leiboldiana.maker.proteins.fasta

2.10. Functional annotation in Blast2Go
----
Like with T.fas, we loaded all CDS information onto blast2go and ran all steps with default settings (viridiplantae DB for blast)
