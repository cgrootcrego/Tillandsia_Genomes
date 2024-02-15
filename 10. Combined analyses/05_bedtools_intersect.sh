te_gff="/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/EDTA/fasciulata/25_scaffolds/tillandsia_fasciculata_assembly.sorted.25_scaffolds.fasta.mod.EDTA.TEanno.no-unknown.gff3"
gff="/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tfas_assembly/assembly_25_scaffolds/Tillandsia_fasciculata_v1.2.edited_allfeatures.25chrom.mRNA.gff"

bedtools intersect -a $gff -b $te_gff -wa -c > Tillandsia_fasciculata_GENE-TE-intersection.ALLgenes.counts.txt
bedtools intersect -a $gff -b $te_gff -wa -wb > Tillandsia_fasciculata_GENE-TE-intersection.ALLgenes.txt
