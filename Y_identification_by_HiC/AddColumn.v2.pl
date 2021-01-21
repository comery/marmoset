#!usr/bin/perl -w
use strict;
use Getopt::Long;
my $header;
GetOptions(
	"h:s"=>\$header,
);

die "perl $0 -h [header] <base table> <add list> <base id col>" unless @ARGV == 3;

my $table = shift;
my $list = shift;

my (%hash, $num);
open(IN, $list) or die $!;
while(<IN>){
	next if(/^#/);
	chomp;
	my @tmp = split /\t/;
	my $id = shift @tmp;
	my $add = join "\t", @tmp;
	$num = @tmp;
	$hash{$id} = $add;
}
close IN;

open(IN, $table) or die $!;
if(defined $header){
	my $line = <IN>;
	chomp $line;
	print "$line\t$header\n";
}
while(<IN>){
	chomp;
	my $id = (split /\t/)[$ARGV[0]-1];
	my $line;
	if(exists $hash{$id}){
		$line = "$_\t$hash{$id}";
	}
	else{
		$line = "$_"."\t-"x$num;
	}
	print "$line\n";
}
