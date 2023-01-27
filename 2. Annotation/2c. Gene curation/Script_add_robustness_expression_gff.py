#!/usr/bin/env python
# coding: utf-8

import sys
from Bio import SeqIO
import re
import pandas as pd


gff = open(sys.argv[1]) # gff file to add tags
genes = open(sys.argv[2]) # a list of all genes to evaluate
expression = pd.read_table(sys.argv[4]) # read count info from RNA-seq mapping per transcript
expression_exon = pd.read_table(sys.argv[5]) # read count info from RNA-seq mapping per exon
directory_fasta = sys.argv[6] # directory with per-CDS fasta files
outputfilename=sys.argv[1].replace(".gff",".with_robustness_orthology.gff")
output=open(outputfilename,'w')
outputfilename2="Robustness_checklist_ALLGenes.txt" # summary list
output2=open(outputfilename2,'w')

def has_startcodon(seq):
    present = seq.startswith("ATG")
    return(present)

def multiple_stopcodon(seq):
    x = len(seq)/3
    dividable = (x - int(x) == 0)
    if dividable == True:
        n = 3 # chunk length
        codons = [seq[i:i+n] for i in range(0, len(seq), n)]
        from collections import Counter
        counts = Counter(codons)
        stopcodon_counts = counts["TAG"]+counts["TAA"]+counts["TGA"]
        if stopcodon_counts == 1:
            multiple_stopcodons = "one_stopcodon"
        elif stopcodon_counts > 1:
            multiple_stopcodons = "multiple_stopcodons"
        elif stopcodon_counts == 0:
            multiple_stopcodons = "no_stopcodon"
    else:
        multiple_stopcodons = "not dividable by 3"
    return(multiple_stopcodons)

#1st step: find orthogroup and determine robustness
# Rule: a gene is ROBUST, when:
# - The gene is expressed across all exons
# OR
# - The gene has a start and stopcodon

first_line = "Gene_id\torthology\tCDS_length\tstartcodon\tstopcodon\texpression\trobust\n"
output2.write(first_line)

gene_dict = {}
for line1 in genes.readlines():
    gene = line1.replace('\n','') # rm the return carriage
    print(gene)
    ortho_info = open(sys.argv[3])
    i = 0
    for line in ortho_info:
        if re.search(gene, line):
            i = 1
            og = line.split('\t')[6]
    if i == 0:
        og = "NO_ORTHOLOGY"
        # Extract info on length, start and stop codon from fasta sequence
    fasta_sequences = SeqIO.parse(open(directory_fasta+gene+":cds.fst"),'fasta')
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
        startcodon = has_startcodon(sequence)
        multiple_stop = multiple_stopcodon(sequence)
        length = len(sequence)
    # Extract info on exon number from gff
    exon_count = 0
    gff_Tfas = open(sys.argv[1])
    for line in gff_Tfas:
        if re.search(gene, line):
            line = line.replace('\n','') # rm the return carriage
            splitted_line = line.split('\t')
            if splitted_line[2] == "exon":
                exon_count=exon_count+1
# Extract info on expression, exons will be marked as not expressed when the average CPM across saples is < 0.001 (1000 counts in total)
    line = expression[expression['Geneid'].str.contains(gene)]
    if line.empty:
        expressed = "not_expressed"
    else:
        sum_Tfas = int(line.iloc[:, 1:37].sum(axis=1))
        if sum_Tfas > 0:
            exons = expression_exon[expression_exon['Geneid'].str.contains(gene)]
            expressed_exons = 0
            for row in exons.iterrows():
                avg = int(line.iloc[:, 1:37].mean(axis=1))
                if avg >= 0.001:
                    expressed_exons = expressed_exons + 1
            expressed = int(expressed_exons)/exon_count
            if expressed == 0:
                expressed = "not_expressed"
        else:
            expressed = "not_expressed"
    if expressed == 1:
        robust = "ROBUST"
    elif startcodon == True and multiple_stop == "one_stopcodon" and length >= 150:
        robust = "ROBUST"
    else:
        robust = "NOT_ROBUST"

    gene_list = [og, length, startcodon, multiple_stop, expressed, robust]
    gene_dict[gene] = gene_list
    to_write = gene+"\t"+og+"\t"+str(length)+"\t"+str(startcodon)+"\t"+multiple_stop+"\t"+str(expressed)+"\t"+robust+"\n"
    output2.write(to_write)

#2nd step, iterate over each line of the gff and check if the ID is in the description dictionnary. If yes, print line + functional descriptions, if no just print the original line
for line2 in gff.readlines():
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t') # split the line regarding the tabulations
    annotation=splitted_line2[8]  # note that in python, it starts from 0, so the 9th field is [8]
    annotation2 = re.split(';|:', annotation)
    ID=annotation2[0].split('=')[1]
    print(splitted_line2[2])
    if splitted_line2[2] == "gene":
        filtered_dict = dict(filter(lambda item: ID in item[0], gene_dict.items()))
        key = list(filtered_dict.keys())[0]
        stats = filtered_dict[key]
        print(stats)
        line_to_write=line2+";Orthogroup="+stats[0]+";Expression="+str(stats[4])+";"+stats[5]+'\n'
        output.write(line_to_write)
    else:
        stats = gene_dict[ID]
        line_to_write=line2+";Orthogroup="+stats[0]+";Expression="+str(stats[4])+";"+stats[5]+'\n'
        output.write(line_to_write)
