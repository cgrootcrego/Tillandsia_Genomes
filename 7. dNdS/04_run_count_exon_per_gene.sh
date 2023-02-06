#!/bin/bash

rm exon_number_per_gene_Tfas.txt
cat list_of_gene_names_Tfas.txt | while read line; do
 gene_name=$(echo "$line" | awk '{print $1}')
 exon_number=$(grep "$gene_name" Tillandsia_fasciculata_v1.2.one-to-one_orthologs.CDSonly.gff | wc -l)
 echo "$gene_name $exon_number" >> exon_number_per_gene_Tfas.txt
done
