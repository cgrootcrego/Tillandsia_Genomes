## Gene annotation in Tillandsia

Gene annotation was run with [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and [MAKER](https://github.com/Yandell-Lab/maker).

# Overview

Gene annotation steps:
1. Prepare input for braker
2. Run braker
3. Prepare input for MAKER
4. Run MAKER 1 / 2 times with training
5. Evaluate gene models
6. Rename MAKER names in gff
7. Functional annotation in Blast2Go

# 1. Preparing files for Braker

The input for braker is a soft-masked version of the assembly (See 2a. TE annotation), the splice information from RNA-seq mapping and single copy proteins predicted by [BUSCO](https://busco.ezlab.org/).

## 1.1. Splice information extraction:

Map RNA-seq to assembly using [STAR](https://github.com/alexdobin/STAR):

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

# 2. Running BRAKER

BRAKER was run with the following command (example of *T. leiboldiana*):

	source /home/fs71400/grootcrego/.bashrc
	conda activate braker-env
	braker.pl --species=[Tleiboldiana] --genome=tillandsia_leiboldiana_assembly.Softmasked.fasta --softmasking --hints=splice_info_from_rna_leiboldiana.gff --prot_seq=single_copy_protein_hints_busco_leiboldiana.faa --prg=gth --gth2traingenes --cores=32 --gff3

BUSCO was run on the augustus.hints.codingseq fastafile with following command:

	busco -i augustus.hints.codingseq -o busco_Tleiboldiana_braker_pred -l /home/fs71400/grootcrego/busco_downloads/lineages/liliopsida_odb10 -m transcriptome -c 16

# 3. Preparing files for MAKER

We used the following input files for MAKER:
- EST library, masked for complex repeats
- gff with complex repeats from EDTA
- BRAKER gene models
- Ananas and monocot curated protein sequences
- TE library from EDTA + Giri Monocot DB


## 3.1. EST library

Transcriptome assembly with [Trinity](https://github.com/trinityrnaseq/trinityrnaseq):

	/apps/trinityrnaseq/2.8.4/Trinity --seqType fq \
	--max_memory 512G --trimmomatic --left [RNA_fastq_R1] \
	--right [RNA_fastq_R2] --CPU 48 \
	--SS_lib_type FR #type of library, here paired end \
	--output [output_dir]

 Our RNA-seq data contained TE sequences, therefore we passed it through Repeatmasker once. We only mask complex repeats, as simple repeats can be part of genic regions:

 	RepeatMasker -pa 4 -e ncbi -dir . -gff -nolow -lib [TE_families_EDTA.fasta] [transcriptome_assembly.fasta]

## 3.2. Complex repeat GFF

We used the repeatmasker util rmOuttoGFF3 to create a gff file based on the masked assembly with only known TEs:

	/apps/repeatmasker/4.0.7/util/rmOutToGFF3.pl [hardmasked_genome] > [complex_repeats.gff3]

We had to reformat the gff3 to make it compatible with MAKER:

  cat [complex_repeats.gff3] | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > [complex_repeats.reformat.gff3]

## 3.3. Protein database

Ananas and monocot swissprot proteins were compiled into one file.

## 3.4. TE database

This is the TE library (fasta) resulting from TE annotation (See 2a. TE annotation)

# 4. Running MAKER

Our annotation strategy with MAKER  is largely based on the steps explained in [this](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2) tutorial
MAKER was run in mpi mode. Options can be found in `maker_R1.ctl` and the running command was:

	mpirun -np 48 /maker/bin/maker -base maker_R1.ctl /maker/maker_bopts.ctl maker/maker_exe.ctl

MAKER was run in consecutive rounds including training of Augustus and SNAP (see tutorial). However, since our final gene sets derive from the first round of MAKER, we won't expand on training here.

# 5. Evaluating gene models

Obtain the final MAKER gff file and fasta sequences:

	~/software/maker/bin/gff3_merge -n -s -d master_datastore_index.log > maker.gff
  	~/software/maker/bin/fasta_merge -d master_datastore_index.log

Count the number of gene models and average length:

	cat leiboldiana_R1.all.maker.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'

Run BUSCO on transcript sequences:

	busco -i maker.transcripts.fasta -o [output_name] -l liliopsida_odb10 -m transcriptome -c 24

Obtain percentage of gene models with an AED score < 0.5:

	perl ~/software/maker/bin/AED_generator.perl -b 0.025 maker.gff

Raw RNAseq data was mapped to the transcripts using STAR. Important: STAR needs an extra option for "small" genomes when indexing. To count the total number of bases:

	grep -v ">" maker.transcripts.fasta | wc | awk '{print $3-$1}'

The default --genomeSAindexNbases is 14, but has to be scaled down for smaller genomes. The calculation is as follows: log2(genome-length)/2 - 1.

Map RNAseq to MAKER transcripts:

	# Genome indexing
	/apps/star/2.7.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir [genome_dir] --genomeFastaFiles maker.transcripts.fasta --runThreadN 24
	# Mapping
  	/apps/star/2.7.3a/bin/Linux_x86_64/STAR --genomeDir [genome_dir] --readFilesIn [RNA_fastq_R1] [RNA_fastq_R2] --outSAMtype BAM Unsorted --outFileNamePrefix [out] --limitBAMsortRAM 10240000000 --runThreadN 24

# 6. Renaming MAKER genemodels

We renamed our final set using MAKER tools:

	# Create a "name map":
	~/software/maker/bin/maker_map_ids --prefix Tlei_v1. --justify 5 maker.gff > maker.name.map
	# Change names in GFF:
	~/software/maker/bin/map_gff_ids maker.name.map maker.gff
	mv maker.gff Tillandsia_leiboldiana_v1.gff

	# For renaming fasta files, taking into account BRAKER names as well:
	awk '$3=="mRNA" {print $0}' maker.gff | cut -f 9 | sed "s/;/\t/g" | cut -f 1,3 | sed 's/=/\t/g' | cut -f 2,4 > mrna_names
	cat mrna_names | while read line ;do   Name=`echo "$line"|awk '{print $1}'`;   echo $Name;   Replace=`echo "$line"|awk '{print $2}'`;   echo $Replace;   sed -i "s/$Name/$Replace/g" maker.name.map; done
	map_fasta_ids maker.name.map maker.transcripts.fasta
	map_fasta_ids maker.name.map maker.proteins.fasta

# 7. Functional annotation in Blast2Go

We loaded all CDS information onto Blast2Go and ran all steps with default settings (viridiplantae DB for blast).
