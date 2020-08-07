#!/usr/bin/perl 
#-------------------------------help-info-start--------------------------------#
=head1 Name 

  FILE: Alignmentout_Sum.pl

=head1 Description
 
   Obtain site information from multiple alignment file (in sam, bam, soap or map type).

=head1 Usage

  perl  Alignmentout_Sum.pl [options] -in <input file> -ref <input reference> -type <type>
  
  -in       <STR>   input alignment file
  -ref      <STR>   input reference sequence  (viral genome)
  -type     <STR>   type of input alignment file: sam, bam, soap, map or pileup
  -outpre   [STR]   out prefix for results, default [myout]
  -outdir   [STR]   out directory for results, default [./]
  -filter   [STR]   filter input aligment file [yes/no], default [no]  
  -ff       [NUM]   filter alignment with [int] or more than [int] differences (mismatches+gaps), default [3]
  -rmdup    [STR]   samtools rmdup, [s/S/], defualt as samtools defualt []
  -samtools [STR]   Samtools directory [defualt]
  -iTools   [STR]   iTools directory [defualt]
  -help     [STR]   print this help to screen

=head1 Example

  perl  Alignmentout_Sum.pl -in in.soap -ref ref.fa -type soap -filter yes -outpre myout -outdir ./

=head1 Contact

  Author: Chengran Zhou, zhouchengran@genomics.cn
          Jinmin Ma, majinmin@genomics.cn
		  Modified by yangchentao, 20200324
  Organization: BGI-CNGB; BGI-SHENZHEN; SCU
  Version: 1.1

=cut
#-------------------------------help-info-end--------------------------------#

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use lib "$Bin/";
use File::Basename qw(basename dirname);
use POSIX qw(strftime);

my ($Need_help,$in,$ref,$name_pre,$outdir,$type,$fil,$filf,$rmdup,$samtools,$iTools);
GetOptions(
    "help!"      => \$Need_help,
	"in=s"		 => \$in,
	"ref=s"		 => \$ref,
	"outpre:s"   => \$name_pre,
	"outdir:s"   => \$outdir,
	"type=s"	 => \$type,
	"filter:s"   => \$fil,
	"ff:i"	     => \$filf,
	"rmdup:s"    => \$rmdup,
	"samtools:s"  => \$samtools,
	"iTools:s"   => \$iTools,
);

die `pod2text $0` if ($Need_help || ! $in || ! $type);

if ($type eq 'pileup'){
	die `pod2text $0` unless (defined $in);
}else{
	die `pod2text $0` unless (defined $in && defined $ref && defined $type);
}

#============================================================================#
#                              Global Variable                               #
#============================================================================#
$samtools ||= "samtools";
$iTools   ||= "iTools";

my $currentdir = `pwd`; chomp $currentdir;
$outdir ||= $currentdir;
mkdir $outdir unless (-d $outdir);

$name_pre ||= "myout";
$fil  ||= "no";
$filf ||= "3";
$rmdup||= "";

my $pileup="$outdir/$name_pre.pileup";

#============================================================================#
#                               Main process                                 #
#============================================================================#

##format
if ($type =~/pileup/){
	$pileup=$in;
	
}else {
	if ($type eq 'soap') {
		system "$iTools Formtools Soap2Sam -InSoap $in -OutSam $outdir/$name_pre.sam";
		$type="sam";$in="$outdir/$name_pre.sam";}
	if ($type eq 'map') {
		system "$iTools Formtools Maq2Sam $in > $outdir/$name_pre.sam";
		$type="sam";$in="$outdir/$name_pre.sam";}
	if ($type eq 'sam') { 
		system "$samtools view -b -S -T $ref -o  $outdir/$name_pre.bam $in";
		$in="$outdir/$name_pre.bam";}
	system "cp $in $outdir/$name_pre.bam " unless (-e "$outdir/$name_pre.bam");
	system "$samtools faidx $ref" unless (-e "$ref.fai");

##filter
	if ($fil=~/yes/i){
		system "$samtools view -F 4 -b $outdir/$name_pre.bam | $samtools sort - $outdir/$name_pre.sort ";	
		system "$samtools fillmd $outdir/$name_pre.sort.bam $ref 2>/dev/null|perl -ne 'print if (/^@/||(/NM:i:(\\d+)/&& \$1<=$filf))'|perl -ne 'print if (/^@/||(/NH:i:(\\d+)/&& \$1==1))' > $outdir/$name_pre.filter.sam ";    #NM  The edit distance from read string to reference
		system "$samtools view -b -S -T $ref -o $outdir/$name_pre.filter.bam $outdir/$name_pre.filter.sam";
	
		if ($rmdup){system "$samtools rmdup -$rmdup $outdir/$name_pre.filter.bam $outdir/$name_pre.filter.rmd.bam 1>>$outdir/std.log 2>>$outdir/std.err";}
		else {system "$samtools rmdup $outdir/$name_pre.filter.bam $outdir/$name_pre.filter.rmd.bam 1>>$outdir/std.log 2>>$outdir/std.err";}
	
		system "$samtools view $outdir/$name_pre.filter.rmd.bam > $outdir/$name_pre.filter.rmd.sam";
		system "$samtools mpileup -AB -f $ref -d 1000000 $outdir/$name_pre.filter.rmd.bam > $outdir/$name_pre.pileup  2>>$outdir/std.err";
		system "$iTools Xamtools stat -InFile $outdir/$name_pre.filter.rmd.bam -Bam -Ref $ref -OutStat $outdir/$name_pre.coverage 1> $outdir/std.log";
	
	} else {system "$samtools mpileup -AB -f $ref -d 1000 $outdir/$name_pre.bam > $outdir/$name_pre.pileup  2>>$outdir/std.err";
		  system "$iTools Xamtools stat -InFile $outdir/$name_pre.bam -Bam -Ref $ref -OutStat $outdir/$name_pre.coverage 1> std.log";}
}


##get results
open (IN,$pileup) or die $!;

##print results & entropy
open (OUT,">$outdir/$name_pre.base.xls") or die $!;
print OUT "Scaffold\tSite\tdepth\tdensity\tA\tT\tG\tC\tPi\tInsert\tDel\tReal depth\tEntropy\tmut_rate\tobserved depth\n";
#print OUT "Scaffold\tSite\tdepth\tdensity\tA\tT\tG\tC\tPi\n";
my %entropy;

# modify by yct 20200324
my $currt_chr = "";
while (<IN>) {
	my (%base,$ndepth,$iddepth,$mut,$density, $npi);
	chomp;
	my ($NC,$site,$letter,$depth,$cov_base) = (split /\s+/, $_, 5);
	$letter = uc ($letter);
	if ($currt_chr eq ""){
		$currt_chr = $NC;
	}elsif ($NC ne $currt_chr){
		print "this Script only support one Scaffold!\n";
		exit();
	}
	if ($depth == 0){
		$cov_base="";
	}
	$cov_base =~ s/\^.|\$//g;
	while ($cov_base =~ /[+-]{1}(\d+)\w+/) {
		$cov_base =~ s/([+-]{1})$1(\w{$1})//;
		my $indel = uc($2);
		$base{$1}{$indel}++;
	}
	$base{A} = $cov_base =~ tr/Aa/Aa/;
	$base{T} = $cov_base =~ tr/Tt/Tt/;
	$base{G} = $cov_base =~ tr/Gg/Gg/;
	$base{C} = $cov_base =~ tr/Cc/Cc/;
	$base{$letter} = $cov_base =~ tr/,./,./;
	$base{skip}= $cov_base =~ tr/></></;
	my $nA=$base{A};
	my $nT=$base{T};
	my $nC=$base{C};
	my $nG=$base{G};
	$ndepth = $nA+$nT+$nC+$nG;
	##depth with indels
	$iddepth = $depth-$base{skip};
	my ($max_base) = (sort{$base{$b}<=>$base{$a}}  ('A','T','G','C',$letter));
	##mutation rate
	$mut = $iddepth ? sprintf "%0.4f",($iddepth - $base{$max_base})/($iddepth) : sprintf "%0.4f",0;
	##nucleotide diversity 	
	$npi = ($ndepth>1)? sprintf "%0.4f",(2*($nA*$nT+$nA*$nC+$nA*$nG+$nT*$nC+$nT*$nG+$nC*$nG)/(($ndepth)*($ndepth-1))) :sprintf "%0.4f",0;
	#$ndepth or $iddepth
	##density
	$density = ($ndepth>1)?sprintf "%0.4f",(log($ndepth)/log(2)):sprintf "%0.4f",0;

	my $Entropy = 0;
	print OUT "$NC\t$site\t$ndepth\t$density";
	foreach my $base ('A','T','G','C') {
		print OUT "\t$base{$base}";
		next if ($iddepth==0);
		my $p = $base{$base}/$iddepth;
		next if $p == 0;
		$Entropy += $p * (log($p)/log(2));
	}
	print OUT "\t$npi";
	my%indel = ('+'=>'',
				'-'=>'',);
	foreach my $symbol ('+','-') {
		next unless exists $base{$symbol};
		next if ($iddepth==0);
		foreach my$indel (keys %{ $base{$symbol} }) {
			$indel{$symbol} .= "$indel($base{$symbol}{$indel});";
			my $p = $base{$symbol}{$indel}/$iddepth;
			$Entropy += $p * log($p)/log(2);
		}
		chop $indel{$symbol};
	}
	$Entropy *= -1;
	$Entropy =sprintf "%0.4f",$Entropy;
	print OUT "\t$indel{'+'}\t$indel{'-'}\t$iddepth\t$Entropy\t$mut\t$depth\n";
}

#system "rm $outdir/$name_pre.bam -f ";
if (-e "$outdir/$name_pre.sam") {
	system "rm $outdir/$name_pre.sam";
}
if (-e "$outdir/$name_pre.filter.sam"){
	system "rm $outdir/$name_pre.filter.sam";
}
if (-e "$outdir/$name_pre.filter.rmd.sam"){
	system "rm $outdir/$name_pre.filter.rmd.sam";
}
#system "rm $outdir/$name_pre.bam  $outdir/$name_pre.sam  $outdir/$name_pre.filter.sam  $outdir/$name_pre.filter.bam $outdir/$name_pre.filter.rmd.sam  -f";
close OUT;
