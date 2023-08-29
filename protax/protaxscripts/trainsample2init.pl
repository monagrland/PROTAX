#!/usr/bin/perl

if (scalar(@ARGV) != 2) {
    die "ERROR (trainsample2init.pl): two arguments required startnode and seqlist\n";
}

$startnode=shift;
while (<>) {
    ($weight,$onode,$priprob,$nodeinfo,$rnode,$trainseq)=split;
    print "$nodeinfo $rnode $trainseq $startnode 0\n";
}
