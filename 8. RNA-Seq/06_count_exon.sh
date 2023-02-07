#!/bin/bash

gff=/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/Tillandsia_fasciculata_v1.2.edited_allfeatures.25chrom.for-exon_featurecounts.gff
output=/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/counts.Tfas_Tlei_6_timepoints.exons.txt
reference=/REFERENCES_TILLANDSIA/Tfas_assembly/assembly_25_scaffolds/Tillandsia_fasciculata_25_scaffolds.fasta
input=/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/2_mapped/*.bam

mkdir tmp

featureCounts -a $gff -o $output -g 'ID' -t exon -G $reference -T 48 -p -s 2 --tmpDir tmp/ $input
