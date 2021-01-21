#!/usr/bin/perl -w
use strict;
use lib '/hwfssz1/ST_DIVERSITY/PUB/USER/zhangpei/perl/pm';
use Pvalue_adj;
use Math::CDF qw(:all);

die "perl $0 <scaf.bed.tab>" unless @ARGV == 1;

my %pvalue;

open(IN, $ARGV[0]) or die $!;
while(<IN>){
	chomp;
	my ($id, $pval) = (split /\t/)[0,11];
	next if($pval eq 'NA');
	$pvalue{$id} = $pval;
}
close IN;

my %qvalue=&BH(%pvalue);

open(IN, $ARGV[0]) or die $!;
while(<IN>){
	chomp;
	my $id = (split /\t/)[0];
	$qvalue{$id} = 'NA' unless(exists $qvalue{$id});
	print "$_\t$qvalue{$id}\n";
}
close IN;

