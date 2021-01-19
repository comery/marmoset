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

die `pod2text $0` if ($Need_help);
#die `pod2text $0` unless (defined $in && defined $ref && defined $type);

if ($type eq 'pileup'){
	die `pod2text $0` unless (defined $in);
}else{
	die `pod2text $0` unless (defined $in && defined $ref && defined $type);
}

die "Input file dose not exist" unless (-e $in);
#die "Input reference file dose not exist" unless (-e $ref);
#============================================================================#
#                              Global Variable                               #
#============================================================================#
$samtools ||= "samtools";
$iTools   ||= "/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Bio_packages/iTools_Code/bin/iTools";

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
my (%base,%depth,%ndepth,%iddepth,@scaf,%npi,%mut,%density,%redun);

while (<IN>) {
	chomp;
	my ($NC,$site,$letter,$depth,$cov_base) = (split)[0,1,2,3,4];
	$letter = uc ($letter);
	push @scaf,$NC if ++$redun{$NC} <2;
	if ($depth == 0){
		$cov_base="";
	}
	$cov_base =~ s/\^.|\$//g;
	while ($cov_base =~ /[+-]{1}(\d+)\w+/) {
		$cov_base =~ s/([+-]{1})$1(\w{$1})//;
		my $indel = uc($2);
		$base{$NC}{$site}{$1}{$indel}++;
	}
	$base{$NC}{$site}{A} = $cov_base =~ tr/Aa/Aa/;
	$base{$NC}{$site}{T} = $cov_base =~ tr/Tt/Tt/;
	$base{$NC}{$site}{G} = $cov_base =~ tr/Gg/Gg/;
	$base{$NC}{$site}{C} = $cov_base =~ tr/Cc/Cc/;
	$base{$NC}{$site}{$letter} = $cov_base =~ tr/,./,./;
	$base{$NC}{$site}{skip}= $cov_base =~ tr/></></;
	my $nA=$base{$NC}{$site}{A};
	my $nT=$base{$NC}{$site}{T};
	my $nC=$base{$NC}{$site}{C};
	my $nG=$base{$NC}{$site}{G};
	$ndepth{$NC}{$site} = $nA+$nT+$nC+$nG;
	$depth{$NC}{$site} = $depth;
	##depth with indels
	$iddepth{$NC}{$site} = $depth-$base{$NC}{$site}{skip};
	my ($max_base) = (sort{$base{$NC}{$site}{$b}<=>$base{$NC}{$site}{$a}}  ('A','T','G','C',$letter));
	##mutation rate
	$mut{$NC}{$site} = $iddepth{$NC}{$site} ? sprintf "%0.4f",($iddepth{$NC}{$site} - $base{$NC}{$site}{$max_base})/($iddepth{$NC}{$site}) : sprintf "%0.4f",0;
	##nucleotide diversity 	
	$npi{$NC}{$site} = ($ndepth{$NC}{$site}>1)? sprintf "%0.4f",(2*($nA*$nT+$nA*$nC+$nA*$nG+$nT*$nC+$nT*$nG+$nC*$nG)/(($ndepth{$NC}{$site})*($ndepth{$NC}{$site}-1))) :sprintf "%0.4f",0;
	#$ndepth or $iddepth
	##density
	$density{$NC}{$site} = ($ndepth{$NC}{$site}>1)?sprintf "%0.4f",(log($ndepth{$NC}{$site})/log(2)):sprintf "%0.4f",0;
}
close IN;

##print results & entropy
open (OUT,">$outdir/$name_pre.base.xls") or die $!;
print OUT "Scaffold\tSite\tdepth\tdensity\tA\tT\tG\tC\tPi\tInsert\tDel\tReal depth\tEntropy\tmut_rate\tobserved depth\n";
my %entropy;

foreach my $NC (sort @scaf) {
	foreach my$site (sort{$a<=>$b} keys %{ $base{$NC} }) {
		my $Entropy = 0;
		print OUT "$NC\t$site\t$ndepth{$NC}{$site}\t$density{$NC}{$site}";
		foreach my $base ('A','T','G','C') {
			print OUT "\t$base{$NC}{$site}{$base}";
			next if ($iddepth{$NC}{$site}==0);
			my $p = $base{$NC}{$site}{$base}/$iddepth{$NC}{$site};
			next if $p == 0;
			$Entropy += $p * (log($p)/log(2));
		}
		print OUT "\t$npi{$NC}{$site}";
		my%indel = ('+'=>'',
					'-'=>'',);
		foreach my $symbol ('+','-') {
			next unless exists $base{$NC}{$site}{$symbol};
			next if ($iddepth{$NC}{$site}==0);
			foreach my$indel (keys %{ $base{$NC}{$site}{$symbol} }) {
				$indel{$symbol} .= "$indel($base{$NC}{$site}{$symbol}{$indel});";
				my $p = $base{$NC}{$site}{$symbol}{$indel}/$iddepth{$NC}{$site};
				$Entropy += $p * log($p)/log(2);
			}
			chop $indel{$symbol};
		}
		$Entropy *= -1;
		$Entropy =sprintf "%0.4f",$Entropy;
		print OUT "\t$indel{'+'}\t$indel{'-'}\t$iddepth{$NC}{$site}\t$Entropy\t$mut{$NC}{$site}\t$depth{$NC}{$site}\n";
	}
}
#system "rm $outdir/$name_pre.bam -f ";
system "rm $outdir/$name_pre.sam  $outdir/$name_pre.filter.sam   $outdir/$name_pre.filter.rmd.sam  -f";
#system "rm $outdir/$name_pre.bam  $outdir/$name_pre.sam  $outdir/$name_pre.filter.sam  $outdir/$name_pre.filter.bam $outdir/$name_pre.filter.rmd.sam  -f";
close OUT;
