import sys
from Bio import SeqIO

fasta_ref = SeqIO.parse(open(sys.argv[1]),'fasta')
window_size = int(sys.argv[2])
outputfilename="GC_content_per{}window_py.txt".format(window_size)
output=open(outputfilename,'w')

def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break

firstline = "chrom\tstart_window\tend_window\tperc_GC\n"
output.write(firstline)

for fasta in fasta_ref:
    name = fasta.id
    sequence = str(fasta.seq)
    start = 0
    end = window_size
    for subseq in chunks(sequence, window_size, window_size):
        c_count = subseq.count("C")
        g_count = subseq.count("G")
        totalcount = c_count + g_count
        perc = (totalcount/window_size)*100
        line_to_write = name+"\t"+str(start)+"\t"+str(end)+"\t"+str(perc)+"\n"
        output.write(line_to_write)
        start = start + window_size
        end = end + window_size
