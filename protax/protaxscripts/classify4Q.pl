#!/usr/bin/perl

use warnings;
use strict;

sub readTaxonomyFile {
    # taxonomy file format: id(integer) pid(integer) level(integer) name(string)
    
    my ($file) = @_;
    my %node2parent = ();
    my %node2children = ();
    my %node2level = ();    
    my %name2node = ();    
    my ($id, $pid, $level, $name, $prior);
    my $maxid=0;
    my $numlevels=0;
    my %node2prior = ();
    open(FD,$file) or die "ERROR (readTaxonomyFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($id,$pid,$level,$name,$prior) = split;
	next if ($id eq ""); # skip empty lines
	$node2parent{$id} = $pid;
	push(@{$node2children{$pid}},$id) if ($id != $pid);
	$node2level{$id} = $level;
	$name2node{$name} = $id;
	$maxid = $id if ($id > $maxid);
	$numlevels = $level if ($level > $numlevels);
	die "ERROR (readTaxonomyFile): cannot find prior (5th column) file '$file'. $!\n" if ($prior eq "");
	$node2prior{$id} = $prior;
    } 
    close(FD);
    return (\%node2parent,\%node2children,\%node2level, \%name2node, $maxid, $numlevels, \%node2prior);
}

sub readSequence2TaxonomyFile {
    my ($file, $taxname2node) = @_;
    my %seqid2taxonomy = ();
    my %node2seqids = ();

    my ($seqid,$taxname);
    open(FD, $file) or die "ERROR (readSequence2TaxonomyFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($seqid,$taxname) = split;
	next if ($seqid eq ""); # skip empty lines
	die "ERROR (readSequence2TaxonomyFile): cannot find taxonomy node for seqid: '$seqid' taxname: '$taxname'.\n" if (!exists($taxname2node->{$taxname}));
	$seqid2taxonomy{$seqid} = $taxname2node->{$taxname};
	push(@{$node2seqids{ $taxname2node->{$taxname} }}, $seqid);
    }
    close(FD);

    return (\%seqid2taxonomy,\%node2seqids);
}

sub readTaxonomyRseqIdFile {
    my ($filename, $rseq_separator) = @_;
    my %node2rseq = ();
    my %node2rseqweight = ();
    
    open(FD, $filename) or die "ERROR (readTaxonomyRseqIdFile): cannot read file '$filename'. $!\n";
    while (<FD>) {
	my ($node, $rseqs, $weights) = split;
	@{$node2rseq{$node}} = split(/$rseq_separator/, $rseqs);
	@{$node2rseqweight{$node}} = split(/$rseq_separator/, $weights) if ($weights ne "");
    }
    close(FD);
    
    return (\%node2rseq, \%node2rseqweight);
}

sub readSeqidFile {
    my ($file) = @_;
    my @seqids = ();
    my $id;
    open(FD,$file) or die "ERROR (readSeqidFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($id)=split;
	push(@seqids,$id);
    }
    close(FD);
    return (\@seqids);
}

sub readSeqSimilarityFile {
    my ($file) = @_;
    my %sim = ();
    my ($sid1,$sid2,$seqsim);

    open(FD,$file) or die "ERROR (readSeqSimilarityFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($sid1,$sid2,$seqsim) = split;
#	$sim{$sid1}{$sid2} = $seqsim;
	$sim{$sid1}{$sid2} = $seqsim/100 - 0.7;	
	$sim{$sid1}{$sid2} = 0 if ($sim{$sid1}{$sid2}<0);
    }
    close(FD);

    return (\%sim);
}

sub readParams {
    my ($file) = @_;
    my @params = ();
    my @postprob = ();

    open (FD,$file) or die "ERROR (readParams): cannot read file '$file'. $!\n";
    my @a;
    while (<FD>) {
	# each line: posterior param1 param2 ... paramM
	@a = split;
	push(@postprob,shift(@a));
	push(@params, [@a]);
    }
    close(FD);

    return (\@params, \@postprob);
}

sub getParamVector {
    my ($params_mcmc, $pp, $mode) = @_;
    $mode = lc($mode);
    
    my @params = ();

    if ($mode eq "map") {
	# MCMC sample with largest pp
	my $i;
	my $mi=0;
	my $max=$pp->[$mi];
	for ($i=0; $i<scalar(@$pp); $i++) {
	    $mi=$i, $max=$pp->[$i] if ($pp->[$i] > $max);
	}
	@params = @{$params_mcmc->[$mi]};
    }
    elsif ($mode eq "mean") {
	# mean
	my $numparams = scalar(@{$params_mcmc->[0]});
	my $samplesize = scalar(@$params_mcmc);
	@params = (0) x $numparams;
	my ($i,$j);
	for ($i=0; $i<$samplesize; $i++) {
	    for ($j=0; $j<$numparams; $j++) {
		$params[$j] += $params_mcmc->[$i]->[$j];
	    }
	}
	for ($j=0; $j<$numparams; $j++) {
	    $params[$j] /= $samplesize;
	}
    }
    elsif ($mode eq "random") {
	my $samplesize = scalar(@$params_mcmc);
	my $i=int(rand($samplesize));
	@params = @{$params_mcmc->[$i]};
    }
    else {
	die "ERROR (getParamVector): paramsmode must be either 'map', 'mean', or 'random' (now was '$mode').\n";
    }

    return (\@params);
}

sub simMaxEml {
    my ($qseqid, $tseqidref, $tseqweightref, $seqsimref, $seq_dont_use_ref) = @_;
    my $max=0;
    my $eml=0;
    my ($sim, $w, $i, $q97);
    my $n=0;

    my @s=();
    foreach my $tseqid (@$tseqidref) {
	next if (exists($seq_dont_use_ref->{$tseqid}));
	$sim=0;
	if (exists($seqsimref->{$qseqid}->{$tseqid})) {
	    $sim = $seqsimref->{$qseqid}->{$tseqid};
	    push(@s,$sim);
	    $max = $sim if ($sim > $max);
	}
	$n++;
    }
    if ($n>1) {
	# weight(for evidence for mislabeling) based on number of refseqs
	$w = 1-$n**(-0.2);
	$i = 1 + int(0.03*$n + 0.5);
	$i = 10 if ($i>10);
	if ($i >= scalar(@s)) {
	    $q97 = 0;
	}
	else {
	    @s=sort {$b <=> $a} @s;	
	    $q97 = $s[$i];
	}
	$eml = $w * ($max - $q97);
    }
    return ($max, $eml);
}

sub createX {
    my ($node, $n2cref, $n2rseqref, $n2rseqweightref, $qseqid, $seqsimref, $n2prior, $unkprior, $seq_dont_use, $debug_mode) = @_;

    my ($i,$cnode,$max,$eml,$wref,$p);

    # store all matrix X elements in one array row by row

    my @X = ([0, 0, 0, 0]);
    # unk weight, here assuming that original priors sum to 1 (as they should)
    #my @xweight = ($unkprior/(1-$unkprior)*$n2prior->{$node});
    my @xweight = (1);

    if (exists($n2cref->{$node})) {
	for ($i=0; $i<(scalar(@{$n2cref->{$node}})); $i++) {
	    $cnode = $n2cref->{$node}->[$i];
	    #push(@xweight, $n2prior->{$cnode});
	    push(@xweight, 1);
	    
	    if (exists($n2rseqref->{$cnode})) {
		if (exists($n2rseqweightref->{$cnode})) {
		    $wref = $n2rseqweightref->{$cnode};
		}
		else {
		    $wref = [];
		}
		# calculate max similarity and eml between query seq and child node's refseqs
		($max, $eml) = simMaxEml($qseqid, $n2rseqref->{$cnode}, $wref, $seqsimref, $seq_dont_use);
		push(@X, [0, 1, $max, $eml]);
	    }
	    else {
		# known taxon without refseqs
		push(@X, [1, 0, 0, 0]);
	    }
	}
    }
    return (\@X, \@xweight);
}

sub logProbMultinomLogistic {
    my ($mref, $mweight, $vref, $vstartind, $dim) = @_;
    my ($i,$j,$maxz,$ezsum);

    # matrix m: ([row1], [row2], ..., [rowN])

    my $nrows = scalar(@$mref);
    my @logp = (0) x $nrows;

    # use logp to first store values of z: 1) z=Xb 2) scale z 3) w*ez/sum(w*ez)
    # z=Xb
    for ($i=0; $i<$nrows; $i++) {
	$logp[$i] = 0;
	for ($j=0; $j<$dim; $j++) {
	    $logp[$i] += $mref->[$i][$j] * $vref->[$j+$vstartind];
	}
	$maxz = $logp[$i] if ($i==0 || $logp[$i]>$maxz);
    }
    # scale z = z - max{z}
    $ezsum = 0;
    for ($i=0; $i<$nrows; $i++) {
	$logp[$i] = $mweight->[$i] * exp($logp[$i] - $maxz) + 1e-10;
	$ezsum += $logp[$i];
    }
    # p = w*ez/sum(w*ez)
    for ($i=0; $i<$nrows; $i++) {
	$logp[$i] = log($logp[$i]/$ezsum); 
    }

    return(\@logp);
}

sub printX {
    my ($mref, $dim, $logprobref) = @_;
    my ($i,$j);

    my $nrows = scalar(@$mref);

    for ($i=0; $i<$nrows; $i++) {
	print "[";
	for ($j=0; $j<$dim; $j++) {
	    printf(" %.3f", $mref->[$i][$j])
	}
	printf(" ]  logprob %.3f\n",$logprobref->[$i]);
    }
    return(0);
}

sub computeNodeProb {
    my ($qseqid, $node, $n2logprob, $n2cref, $n2level, $n2rseqid, $n2rseqweight, $params, $betadim, $prevlevel_logprob, $logprob_th, $seqsimref, $n2prior, $unkprior, $seq_dont_use, $debug_mode) = @_;

    my ($i, $cnode);

    print "computeNodeProb '$node' level $n2level->{$node} parentLogprob $prevlevel_logprob\n" if ($debug_mode > 1);

    my ($X, $Xweight) = createX($node, $n2cref, $n2rseqid, $n2rseqweight, $qseqid, $seqsimref, $n2prior, $unkprior, $seq_dont_use, $debug_mode);
    # my $params_startind = 1 + $n2level->{$node} * $betadim;
    # single layer parameters:
    my $params_startind = 1;
    my $logprobref = logProbMultinomLogistic($X, $Xweight, $params, $params_startind, $betadim);
    # this level's logprobs are conditional on previous level, combine
    for ($i=0; $i<scalar(@$logprobref); $i++) {
	$logprobref->[$i] += $prevlevel_logprob;
    }

    printX($X,$betadim,$logprobref) if ($debug_mode>2);

    # store logprobs of this node
    $n2logprob->{$node} = $logprobref;
    # dim(logprob) == nrow(X), first row for unk 
    for ($i=1; $i<scalar(@$logprobref); $i++) {
	# continue if child node prob is above threshold and child has child nodes
	# -1 because first X row for unk and child node index starts from 0 
	$cnode = $n2cref->{$node}->[$i-1];
	if (($logprobref->[$i] >= $logprob_th) && exists($n2cref->{$cnode})) {
	    computeNodeProb($qseqid, $cnode, $n2logprob, $n2cref, $n2level, $n2rseqid, $n2rseqweight, $params, $betadim, $logprobref->[$i], $logprob_th, $seqsimref, $n2prior, $unkprior, $seq_dont_use, $debug_mode);
	}
    }
    return (0);
}

sub readSeqidFileWithStartProb {
    my ($file) = @_;
    my @seqids = ();
    my @startnodes = ();
    my ($id,$startnode,$startprob,$n,$i,@a);
    open(FD,$file) or die "ERROR (readSeqidFileWithStartProb): cannot read file '$file'. $!\n";
    my $linenum=0;
    while (<FD>) {	
	@a=split;
	$linenum++;
	$id = shift(@a);
	push(@seqids,$id);
	$n = scalar(@a)/2;
	if (2*$n != scalar(@a)) {
	    die ("ERROR (readSeqidFileWithStartProb): extra/missing columns (seqid, 0 or more [startnode,startprob] pairs), line $linenum, file '$file'.\n");
	}
	my @snodes = ();
	for($i=0; $i<$n; $i++) {
	    $startnode=$a[2*$i];
	    $startprob=$a[2*$i +1];
	    push(@snodes,[$startnode, $startprob]);
	}
	push(@startnodes, \@snodes);
    }
    close(FD);
    return (\@seqids, \@startnodes);
}

if (scalar(@ARGV) < 11) {
    die "usage: classify5sa.pl unkprior seqid_file(seqid,start node probs) taxonomyfile seqid2taxonomy_file taxonomy2rseqid_file params_file params_mode seqsimfile max_out_num(0 for all) prob_th outfile debug_mode\n";
}

my $unkprior = shift;
my $seqidfile = shift;
my $taxonomyfile=shift;
my $seqid2taxfile = shift;
my $tax2ridfile=shift;
my $paramsfile = shift;
my $paramsmode = shift;
my $seqsimfile=shift;
my $maxoutnum = shift;
my $prob_th = shift;
my $outfile = shift;
my $debug_mode = shift;
$debug_mode = 0 if (!defined($debug_mode));

my ($qidref,$qidstart_taxid_probref) = readSeqidFileWithStartProb($seqidfile);
my ($n2pref,$n2cref,$n2lref,$tname2nref,$maxid,$numlevels,$n2prior) = readTaxonomyFile($taxonomyfile);
my ($s2tref, $n2sidref) = readSequence2TaxonomyFile($seqid2taxfile,$tname2nref);
my ($n2ridref, $n2rweightref) = readTaxonomyRseqIdFile($tax2ridfile, ',');
my ($params_mcmc,$pp) = readParams($paramsfile);
my $paramvec = getParamVector($params_mcmc,$pp,$paramsmode);
my $seqsimref = readSeqSimilarityFile($seqsimfile);

$numlevels=1;
my $logprob_th = log($prob_th + 1e-100);
my $betadim = (scalar(@$paramvec) - 1)/$numlevels;

open(FD,">$outfile") or die "ERROR: cannot write to '$outfile'. $!\n";
my ($ap,$node,$startlogprob);
my ($k,$numseqs,$qseqid, %seq_dont_use);

$numseqs = scalar(@$qidref);
for ($k=0; $k<$numseqs; $k++) {
    $qseqid = $qidref->[$k];

    print "classify '$qseqid'\n" if ($debug_mode);
    if (scalar(@{$qidstart_taxid_probref->[$k]}) == 0) {
	print FD "$qidref->[$k]\n";
	next;
    }

    my %n2logprob = ();    
    %seq_dont_use = ();

    foreach $ap (@{$qidstart_taxid_probref->[$k]}) {
	# $ap->[0] is node id
	next if (($ap->[0] =~ /,unk$/) or !exists($n2lref->{$ap->[0]}));
#	if (!exists($n2lref->{$ap->[0]})) {
#	    die ("ERROR: cannot find taxonomy node '$ap->[0]' for qseq '$qseqid'\n");
#	}
	print " start ($ap->[0], $ap->[1])\n" if ($debug_mode > 1);
	
	$node = $ap->[0];
	$startlogprob = $ap->[1];
	computeNodeProb($qseqid, $node, \%n2logprob, $n2cref, $n2lref, $n2ridref, $n2rweightref, $paramvec, $betadim, $startlogprob, $logprob_th, $seqsimref, $n2prior, $unkprior, \%seq_dont_use, $debug_mode);
    }
    
    my @result = ();
    my ($i, $cnode);
    # output sorted leaf nodes with logprob > threshold
    foreach $node (keys %n2logprob) {

	# combine unk and branches without refseqs

	my $norefsum = exp($n2logprob{$node}->[0]);

	# sum probabilities of child nodes without refseqs

	for ($i=1; $i<scalar(@{$n2logprob{$node}}); $i++) {
	    # -1 because first logprob element for unk and child node index starts from 0 
	    $cnode = $n2cref->{$node}->[$i-1];
	    if (!exists($n2ridref->{$cnode})) {
		$norefsum += exp($n2logprob{$node}->[$i]);
	    }
	}
	if ($norefsum >= $prob_th) {
	    push(@result, [log($norefsum), "$node,unk"]);
	}

	for ($i=1; $i<scalar(@{$n2logprob{$node}}); $i++) {
	    # -1 because first logprob element for unk and child node index starts from 0 
	    $cnode = $n2cref->{$node}->[$i-1];
	    if (exists($n2ridref->{$cnode}) && !exists($n2cref->{$cnode})) {
		# cnode has refseqs and is leaf node
		if ($n2logprob{$node}->[$i] >= $logprob_th) {
		    push (@result, [$n2logprob{$node}->[$i], $cnode]);
		}
	    }
	}
    }
    # sort results based on logprob in descending order
    print FD "$qidref->[$k]";
    my $count=0;
    foreach $i (sort {$b->[0] <=> $a->[0]} @result) {
	if (($count < $maxoutnum) || ($maxoutnum==0)) {
	    # my $p = exp($i->[0]);
	    # print FD " $i->[1] " . sprintf("%.6g",$p);
	    print FD " $i->[1] " . sprintf("%.6g",$i->[0]);
	}
	$count++;
    }
    print FD "\n";

    # print all nodes with logprob > threshold
    if ($debug_mode > 1) {
	foreach $node (sort {$a <=> $b} keys %n2logprob) {
	    if ($n2logprob{$node}->[0] >= $logprob_th) {
		print "   UNK $node $n2logprob{$node}->[0]\n";
	    }
	    for ($i=1; $i<scalar(@{$n2logprob{$node}}); $i++) {
		# -1 because first logprob element for unk and child node index starts from 0 
		$cnode = $n2cref->{$node}->[$i-1];
		if ($n2logprob{$node}->[$i] >= $logprob_th) {
		    print "   ($i) $cnode $n2logprob{$node}->[$i]\n";
		}
	    }
	}
    }
}
close(FD);

