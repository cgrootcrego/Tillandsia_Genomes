#!/bin/bash

list_of_orthologs="/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/dnds/dnds_paralogs/orthologs-112.perOG.txt"
macse="/home/fs71400/grootcrego/software/macse_v2.05.jar"
fasta_dir="/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/dnds/dnds_paralogs/fasta_seq_per_OG_112"

cd alignments_112/
cat $list_of_orthologs | while read line; do
	orthogroup=`echo "$line"|awk '{print $1}'`
	java -jar $macse -prog alignSequences -seq $fasta_dir/${orthogroup}_Tlei-copy1.fasta -local_realign_init 1 -local_realign_dec 1 ;
	java -jar $macse -prog alignSequences -seq $fasta_dir/${orthogroup}_Tlei-copy2.fasta -local_realign_init 1 -local_realign_dec 1 ;
done
