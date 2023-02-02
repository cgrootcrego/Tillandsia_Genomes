#!/usr/bin/perl

# file: gc_content.pl
# Calculate the GC content in a sliding window

use Getopt::Long;

my $WINDOW;

# parse command-line parameters
die "Usage: $0 [--window XX] [file1 file2...]" unless GetOptions('window=i' => \$WINDOW);

$WINDOW ||= 100000;  # default to a window size of 10 unless otherwise specified

while (<>) {
  chomp;
  my ($id,$sequence) = split /\t+/;    # split ID from sequence
  print "= $id =\n";                   # little ad hoc header
  gc_content($sequence,$WINDOW);       # gc_content is a local function
  print "\n";                          # extra line for neatness
}


sub gc_content {
  my $seq = shift;                        # sequence
  my $win = shift;                        # window
  for (my $i = 0; $i <= length($seq) - $win; $i+=100000) { # slide across sequence one bp at a time
    my $segment = substr($seq,$i,$win);  # fetch out a segment of the sequence $win bp long starting at $i
#    my $n_segment = substr($seq,$i,$win);  # fetch out a segment of the sequence $win bp long starting at $i
#    my $atgc_segment = substr($seq,$i,$win);  # fetch out a segment of the sequence $win bp long starting at $i
    my $n_count = $segment =~ tr/Nn/Nn/;
    my $atgcTE_count = $segment =~ tr/atgc/atgc/;  # trick alert -- see manual entry for tr////
    my $atgcnonTE_count = $segment =~ tr/ATGC/ATGC/;  # trick alert -- see manual entry for tr////
    my $atgc_count = $segment =~ tr/ATGCatgc/ATGCatgc/;  # trick alert -- see manual entry for tr////
    my $gc_count = $segment =~ tr/GCgc/GCgc/;  # trick alert -- see manual entry for tr////
    my $gcTE_count = $segment =~ tr/gc/gc/;  # trick alert -- see manual entry for tr////
    my $gcnonTE_count = $segment =~ tr/GC/GC/;  # trick alert -- see manual entry for tr////
#    my $others_count = $segment =~ tr/MRWSYKmrwsyk/MRWSYKmrwsyk/;  # trick alert -- see manual entry for tr////
    my $n_content = $n_count/length($segment);
    if ($n_count<49999) {
      my $gc_content = $gc_count/$atgc_count;
      my $TE_content = $atgcTE_count/$atgc_count;
      my $gcTE_content = $gcTE_count/$atgcTE_count;
      my $gcnonTE_content = $gcnonTE_count/$atgcnonTE_count;
      print $i+1,"\t",$gc_count,"\t",$atgc_count,"\t",$gcTE_count,"\t",$atgcTE_count,"\t",$n_count,"\t",$TE_content,"\t",$gc_content,"\t",$gcTE_content,"\t",$gcnonTE_content,"\t",$n_content,"\n"}
    if ($n_count>49999) {
      print $i+1,"\t",$gc_count,"\t",$atgc_count,"\t",$gcTE_count,"\t",$atgcTE_count,"\t",$n_count,"\t","NA","\t","NA","\t","NA","\t","NA","\t",$n_content,"\n"}
}
}
