#!/usr/bin/perl

if (scalar(@ARGV) != 2) {
    die "ERROR (testsample2init.pl): two arguments required startnode and seqlist\n";
}

$startnode=shift;
while (<>) {
    ($id)=split;
    print "$id $startnode 0\n";
}
