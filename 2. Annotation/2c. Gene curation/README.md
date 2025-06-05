# Gene curation

Using orthology and expression information (see 3. Orthology, 8. RNA-seq), we added tags to our gene annotation to flag potentially unreliable genes.

We calculated the number of exons, length, absence of start and stop codon (completeness), and expression per gene in both species. This was done with `Script_add_robustness_expression_gff`, which will output the resulting information of orthology, completeness and expression per gene in a list, and also add this information in tags to the gff file.

Rule: a gene is ROBUST, when:
- The gene is expressed across all exons
OR
- The gene has a start and stopcodon

Tags included the orthology of the gene, with "NO_ORTHOLOGY" when the gene was not included in orthology analyses or assigned to an orthogroup. Expression information is shown as a proportion of exons expressed. Lowly expressed exons were considered not expressed.

The curation summary displaying the per-gene information on orthology, completeness and expression can be found in `Robustness_checklist_ALLGenes.Tillandsia_fasciculata.txt` and `Robustness_checklist_ALLGenes.Tillandsia_leiboldiana.txt`. The GFF files containing curation tags are `Tillandsia_fasciculata_v1.gff` and `Tillandsia_leiboldiana_v1.gff`.

To match the chromosome names in the gff file with the NCBI fasta headers for the top 25 /top 26 scaffolds (main assemblies), use `Tillandsia_fasciculara_v1.Top25.Chromosome_names.txt` and `Tillandsia_leiboldiana_v1.Top26.Chromosome_names.txt`. Secondary scaffolds retain their original scaffolds names ("Scaffold_N") in the fasta headers, e.g.: >JAQOTM010000045.1 Tillandsia fasciculata voucher WU:0013642 isolate MHJB-B1840 *Scaffold_1318*, whole genome shotgun sequence.

The fasta sequences are the two reference genomes with fasta headers that match the GFF files.
