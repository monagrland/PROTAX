#!/usr/bin/perl

while (<>) {
    ($nid,$pid,$level,$name,$prior)=split;
    if ($name eq "Fungi") {
	$prior = 1.0;
    }
    print "$nid\t$pid\t$level\t$name\t$prior\n" if ($name ne "nonFungi");
}
