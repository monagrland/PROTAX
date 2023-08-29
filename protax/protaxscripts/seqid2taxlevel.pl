#!/usr/bin/perl

$targetlevel=shift;

# use full length names for taxa
while (<>) {
    ($id,$level,$nimi) = split;
    @a=split(/,/,$nimi);
    if ($level >= $targetlevel) {
	$tax=join(',',@a[0 .. ($targetlevel-1)]);
	print "$id\t$tax\n";
    }
}
