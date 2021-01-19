#!/usr/bin/perl -w
use strict;


open IN, shift;
my %hash;
my $total = 0;
my $mummer_snp = 0;
my $good = 0;
my $hard2id = 0;
my $hard2id_of_mummer = 0;
while(<IN>){
	next if (/^#/);
    chomp;
    my @a = split /\t/;
	$total ++ ;
	my $mat = $a[2];
	my $pat = $a[4];
    my $tag = $a[50];
	if ($mat ne "-" && $pat ne "-" && $mat ne $pat) {
		$mummer_snp ++;
		if ($tag eq "hard2id"){
			$hard2id_of_mummer ++;
		}
	}
	if ($tag eq "No_doubt"){
		$good ++;

	}elsif ($tag eq "hard2id"){
		$hard2id ++;
	}elsif ($tag =~ /-/){
		my @b = split/-/, $tag;
		my $mat_tag = $b[0];
		my $pat_tag = $b[1];
		if ($mat_tag =~ /sequencing_error/){
			$hash{'mat'}{'sequencing'} ++;
		}elsif ($mat_tag =~ /polishing_error/){
			$hash{'mat'}{'polishing'} ++;

		}elsif ($mat_tag =~ /assigning_error/){
			$hash{'mat'}{'assigning'} ++;
		}
		if ($pat_tag =~ /sequencing_error/){
			$hash{'pat'}{'sequencing'} ++;
		}elsif ($pat_tag =~ /polishing_error/){
			$hash{'pat'}{'polishing'} ++;

		}elsif ($pat_tag =~ /assigning_error/){
			$hash{'pat'}{'assigning'} ++;
		
		}
	}else{
		if ($tag =~ /error/){
			my @c = split/_/,$tag;
			my $type = $c[0];
			$hash{$c[2]}{$type} ++;

		}
	}
}

close IN;
print "total\thard2id\tmummer_snp\tgood\thard2id_of_mummer\t";
print "mat_assigning\tpat_assigning\tmat_polishing\tpat_polishing\tmat_sequencing\tpat_sequencing\n";
print "$total\t$hard2id\t$mummer_snp\t$good\t$hard2id_of_mummer\t";
foreach my $a('assigning', 'polishing', 'sequencing'){
	foreach my $b('mat', 'pat'){
		if (exists $hash{$b}{$a}){
			print "$hash{$b}{$a}\t";
		}else{
			print "0\t";
		}
	}
}
print "\n";

