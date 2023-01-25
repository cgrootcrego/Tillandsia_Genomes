# Annotation of repetitive content in *T. fasciculata* and *T. leiboldiana*

TE annotation was performed using the Extensive *de novo* TE Annotator ([EDTA](https://github.com/oushujun/EDTA)) for both species.

# Original annotation, used for masking during gene prediction

Originally, EDTA was run in version 1.7.8 for Tillandsia fasciculata and version 1.8.5 for Tillandsia leiboldiana. The original annotation didn't include whole genome annotation as these were not available in those versions. Annotation was run with the command:

    perl EDTA.pl --genome tillandsia_leiboldiana_assembly.pilon.upper.sorted.fasta \
	--cds pineapple.20150427.cds.wo_transposons.fasta \
	--sensitive --threads 16

*A. comosus* coding sequences were used to filter out genes falsely identified as TEs in the final step of annotation.

# Detailed annotation of the full assembly

To report the repetitive content of both assemblies, annotation was rerun with extra steps to obtain the most precise result. The EDTA version used here was 1.9.6:

	perl $edta --genome $genome --cds $cds \
	--sensitive 1 --anno 1 --threads 24

Anno will simply create a whole genome annotation file (GFF3) and basic statistics which were then used for analyses of genic content. This version of the annotation was *not* used for downstream gene annotation.

# Masking genomes

Repeatmasker was run with the following command (softmasking):

`RepeatMasker -dir $dir -gff -pa 4 -e ncbi -nolow -xsmall -lib $telib $fasta`

For hardmasking, remove option `-xsmall`
