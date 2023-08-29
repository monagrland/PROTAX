#!/usr/bin/perl

$taxonomyfile=shift;
open(FD,$taxonomyfile);
while (<FD>) {
    ($nid,$pid,$level,$name)=split;
    $nid2tname{$nid} = $name;
}
close(FD);

while (<>) {
    @a=split;
    $seqid = shift(@a);
    @out = ($seqid);
    $n=scalar(@a);
    for ($i=0; $i<$n; $i=$i+2) {
	$nid=$a[$i];
	$prob=sprintf("%.3g",exp($a[$i+1]));
	$unk="";
	if ($nid =~ /,unk$/) {
	    $nid =~ s/,unk$//;
	    $unk=",unk";
	}
	push(@out,"$nid2tname{$nid}$unk");
	push(@out,$prob);
    }
    $s = join("\t",@out);
    print "$s\n";
}
