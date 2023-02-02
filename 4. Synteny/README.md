# Studying synteny and rearrangements between *T. fasciculata* and *T. leiboldiana*

Using orthologous gene pairs, whole genome alignment and local alignments, we investigated synteny and identified rearrangements between both species and also *Ananas comosus*.

# Visualizing synteny with GENESPACE


# Whole genome alignment between Tfas, Tlei and A. comosus

We ran a whole-genome alignment between *A. comosus*, *T. fasciculata* and *T. leiboldiana* (only main scaffolds) to obtain a second source of evidence for rearrangements. See `run_nucmer.sh`. This alignment was then visualized with dotplots on [Dot](https://github.com/marianattestad/dot).

# Local alignments between Tfas and Tlei using LastZ

We performed LastZ alignments on hard-masked genomes (see 2a. TE annotation) to investigate rearrangement breakpoints more in-depth.

Lastz alignments were made to each main scaffold and parallelized. Create fasta files for each main scaffold:
````
    cat tillandsia_fasciculata_chrnames.txt | while read line ; do  
     Name=`echo "$line"|awk '{print $1}'`;  
     echo $Name;  
     Replace=`echo "$line"|awk '{print $2}'`;  
     echo $Replace;  
     seqkit grep -p $Name Tillandsia_fasciculata_25_scaffolds.fasta.masked > Tfas_$Replace.fasta.masked;  
     sed  -i "s/$Name/Tfas_$Replace/g" Tfas_$Replace.fasta.masked;
    done
````
To replace chromosome names in the full T.lei genome, I ran:
````
    cat Tillandsia_leiboldiana_26_scaffolds_chrnames.txt | while read line ; do  
     Name=`echo "$line"|awk '{print $1}'`;  
     echo $Name;  
     Replace=`echo "$line"|awk '{print $2}'`;  
     echo $Replace;  
     sed  -i "s/$Name/Tlei_$Replace/g" Tillandsia_leiboldiana_26_scaffolds.fasta.masked;
    done
````
LastZ alignments were run with a slurm array, see `run_lastz.sh`.

The alignments were then processed with [scripts](https://shorturl.at/xLS15) from Leroy et al. 2021:
````
    for i in Tlei_vs_Tfas_chr*; do echo $i;  
     python2 convertMaftoCoordinates.py $i > ${i%.maf}.coord;
    done
    bash script_generate_matrix_alignblastz_results.sh test_out Tillandsia_leiboldiana_26_scaffolds.fasta.masked
````

The script `script_filter_alignments_by_uniqness_threshold.py` will eliminate all alignments which overlap with an alignment from a different chromosome above a certain threshold. This script was applied both on the length-filtered and non-filtered coord file:
````
    python2 script_filter_alignments_by_uniqness_threshold.py Tlei_vs_Tfas_allchrom_lastz.coord 0.9
````
Visualization of the alignments was done with `Visualization_aln_lastz.R`
The LastZ alignments seem to support the rearrangements observed in dotplot and orthology analyses. In some cases, the breakpoints are clear-cut, in others not so much. I decided to zoom in on the clear breakpoints to find out a range in which the breakpoint may exist, and then look at the pacbio alignments in that region.

# dN/dS rates across rearranged chromosomes

To see what sort of effect rearrangements may have had on nearby genes, we decided to map the dN/dS rates obtained in 6. to the chromosomal location. We also use this opportunity to obtain some general stats and distributions of dN/dS.

First, I filtered all alignments 5 or less variant sites (1286 out of 10,363) to avoid skew in these distributions. I calculated the number of variant sites using [AMAS](https://github.com/marekborowiec/AMAS): `python AMAS.py summary -f phylip -d dna -i *.phylip` which created summary.txt. Then I selected every orthogroup with more than 5 variant sites: `awk 'BEGIN { FS = "\t" } ; $7 > 5 {print $0}' summary.txt | cut -f 1 | sed 's/_NT.phylip//g' > orthogroups_to_keep.5variant-sites.txt`.
