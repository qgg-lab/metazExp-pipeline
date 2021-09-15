#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--tfGeneFamily \\\n" .
		"--T2H T_DAS_2h.vs.0h.tsv \\\n" .
		"--T10H T_DAS_10h.vs.0h.tsv \\\n" .
		"--V2H V_DAS_2h.vs.0h.tsv \\\n" . 
		"--V10H V_DAS_10h.vs.0h.tsv \\\n" .
		"--outputStat \n";
	exit;
}

my ($tfGeneFamily, $outputStat);
my ($T2H, $T10H, $V2H, $V10H);
GetOptions(
        'tfGeneFamily=s'=>\$tfGeneFamily,
	'T2H=s'=>\$T2H,
	'T10H=s'=>\$T10H,
	'V2H=s'=>\$V2H,
	'V10H=s'=>\$V10H,
	'outputStat=s'=>\$outputStat,
);

my ($line, @tmp);
my (%tfdbTrsptFamily);

# 读取T2H, T10H, V2H, V10H中gene发生各种AS类型的数量
my (%geneWithAsNum, $geneWithAsNumHref, @tmp);
my ($asId, $orthoAsId, $control, $treatment, $dltPsi, $PValue, $QValue, $geneId, $orthoGeneId, $symbol, $GO, $Pathway);
$geneWithAsNumHref = \%geneWithAsNum;
open FF, "<$T2H";
# asId    orthoAsId       control treatment       dltPsi  P-Value Q-Value geneId  orthoGeneId     symbol  GO      Pathway
# OSJGA5SS0000000325      NA      0.3498  0.181666666666667       0.168133333333333       2.69417386568966e-05    0.0241272250252061      Os08g0385900    OG0001626       Os08g0385900    OG0001626       GO:0003676      ko03040
<FF>;
while($line=<FF>){
	chomp($line);
	($asId, $orthoAsId, $control, $treatment, $dltPsi, $PValue, $QValue, $geneId, $orthoGeneId, $symbol, $GO, $Pathway) = split(/\t/, $line);
	if($asId=~/A3SS\d+/){
		$geneWithAsNumHref->{$geneId}->{"T2H"}->{"A3SS"}++;
	}elsif($asId=~/A5SS\d+/){
		$geneWithAsNumHref->{$geneId}->{"T2H"}->{"A5SS"}++;
	}elsif($asId=~/SE\d+/){
		$geneWithAsNumHref->{$geneId}->{"T2H"}->{"SE"}++;
	}elsif($asId=~/RI\d+/){
		$geneWithAsNumHref->{$geneId}->{"T2H"}->{"RI"}++;
	}elsif($asId=~/MXE\d+/){
		$geneWithAsNumHref->{$geneId}->{"T2H"}->{"MXE"}++;
	}
}
close FF;
open FF, "<$T10H";
# asId orthoAsId control treatment dltPsi P-Value Q-Value geneId orthoGeneId symbol GO Pathway
<FF>;
while($line=<FF>){
	chomp($line);
	($asId, $orthoAsId, $control, $treatment, $dltPsi, $PValue, $QValue, $geneId, $orthoGeneId, $symbol, $GO, $Pathway) = split(/\t/, $line);
	if($asId=~/A3SS\d+/){
		$geneWithAsNumHref->{$geneId}->{"T10H"}->{"A3SS"}++;
	}elsif($asId=~/A5SS\d+/){
		$geneWithAsNumHref->{$geneId}->{"T10H"}->{"A5SS"}++;
	}elsif($asId=~/SE\d+/){
		$geneWithAsNumHref->{$geneId}->{"T10H"}->{"SE"}++;
	}elsif($asId=~/RI\d+/){
		$geneWithAsNumHref->{$geneId}->{"T10H"}->{"RI"}++;
	}elsif($asId=~/MXE\d+/){
		$geneWithAsNumHref->{$geneId}->{"T10H"}->{"MXE"}++;
	}
}
close FF;

open FF, "<$V2H";
# asId orthoAsId control treatment dltPsi P-Value Q-Value geneId orthoGeneId symbol GO Pathway
<FF>;
while($line=<FF>){
	chomp($line);
	($asId, $orthoAsId, $control, $treatment, $dltPsi, $PValue, $QValue, $geneId, $orthoGeneId, $symbol, $GO, $Pathway) = split(/\t/, $line);
	if($asId=~/A3SS\d+/){
		$geneWithAsNumHref->{$geneId}->{"V2H"}->{"A3SS"}++;
	}elsif($asId=~/A5SS\d+/){
		$geneWithAsNumHref->{$geneId}->{"V2H"}->{"A5SS"}++;
	}elsif($asId=~/SE\d+/){
		$geneWithAsNumHref->{$geneId}->{"V2H"}->{"SE"}++;
	}elsif($asId=~/RI\d+/){
		$geneWithAsNumHref->{$geneId}->{"T2H"}->{"RI"}++;
	}elsif($asId=~/MXE\d+/){
		$geneWithAsNumHref->{$geneId}->{"V2H"}->{"MXE"}++;
	}
}
close FF;
open FF, "<$V10H";
# asId orthoAsId control treatment dltPsi P-Value Q-Value geneId orthoGeneId symbol GO Pathway
<FF>;
while($line=<FF>){
	chomp($line);
	($asId, $orthoAsId, $control, $treatment, $dltPsi, $PValue, $QValue, $geneId, $orthoGeneId, $symbol, $GO, $Pathway) = split(/\t/, $line);
	if($asId=~/A3SS\d+/){
		$geneWithAsNumHref->{$geneId}->{"V10H"}->{"A3SS"}++;
	}elsif($asId=~/A5SS\d+/){
		$geneWithAsNumHref->{$geneId}->{"V10H"}->{"A5SS"}++;
	}elsif($asId=~/SE\d+/){
		$geneWithAsNumHref->{$geneId}->{"V10H"}->{"SE"}++;
	}elsif($asId=~/RI\d+/){
		$geneWithAsNumHref->{$geneId}->{"V10H"}->{"RI"}++;
	}elsif($asId=~/MXE\d+/){
		$geneWithAsNumHref->{$geneId}->{"V10H"}->{"MXE"}++;
	}
}
close FF;

my $family;
my ($geneId, $trsptId,  $tfdbTrsptId, $family, $type, $T2H, $T10H, $V2H, $V10H);
open WW, ">$outputStat";
open FF, "<$tfGeneFamily";
# Os12g0233800    bZIP
# Os05g0121600    AP2-EREBP
print WW join("\t", "geneId", "family", "Thaibonnet 2h/0h", "Thaibonnet 10h/0h", "Volano 2h/0h", "Volano 10h/0h") . "\n";
while($line=<FF>){
	chomp($line);
	($geneId, $family) = split(/\t/, $line);

	$T2H = "";
	if(exists($geneWithAsNumHref->{$geneId}->{"T2H"}->{"A3SS"})){
		$T2H .= $geneWithAsNumHref->{$geneId}->{"T2H"}->{"A3SS"} . "A3SS,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"T2H"}->{"A5SS"})){
		$T2H .= $geneWithAsNumHref->{$geneId}->{"T2H"}->{"A5SS"} . "A5SS,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"T2H"}->{"MXE"})){
		$T2H .= $geneWithAsNumHref->{$geneId}->{"T2H"}->{"MXE"} . "MXE,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"T2H"}->{"RI"})){
		$T2H .= $geneWithAsNumHref->{$geneId}->{"T2H"}->{"RI"} . "RI,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"T2H"}->{"SE"})){
		$T2H .= $geneWithAsNumHref->{$geneId}->{"T2H"}->{"SE"} . "SE,"
	}
	if($T2H eq ""){
		$T2H = "-";
	}else{
		 $T2H = substr($T2H, 0, length($T2H) - 1);
	}

	$T10H = "";
	if(exists($geneWithAsNumHref->{$geneId}->{"T10H"}->{"A3SS"})){
		$T10H .= $geneWithAsNumHref->{$geneId}->{"T10H"}->{"A3SS"} . "A3SS,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"T10H"}->{"A5SS"})){
		$T10H .= $geneWithAsNumHref->{$geneId}->{"T10H"}->{"A5SS"} . "A5SS,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"T10H"}->{"MXE"})){
		$T10H .= $geneWithAsNumHref->{$geneId}->{"T10H"}->{"MXE"} . "MXE,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"T10H"}->{"RI"})){
		$T10H .= $geneWithAsNumHref->{$geneId}->{"T10H"}->{"RI"} . "RI,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"T10H"}->{"SE"})){
		$T10H .= $geneWithAsNumHref->{$geneId}->{"T10H"}->{"SE"} . "SE,"
	}
	if($T10H eq ""){
		$T10H = "-";
	}else{
		 $T10H = substr($T10H, 0, length($T10H) - 1);
	}

	$V2H = "";
	if(exists($geneWithAsNumHref->{$geneId}->{"V2H"}->{"A3SS"})){
		$V2H .= $geneWithAsNumHref->{$geneId}->{"V2H"}->{"A3SS"} . "A3SS,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"V2H"}->{"A5SS"})){
		$V2H .= $geneWithAsNumHref->{$geneId}->{"V2H"}->{"A5SS"} . "A5SS,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"V2H"}->{"MXE"})){
		$V2H .= $geneWithAsNumHref->{$geneId}->{"V2H"}->{"MXE"} . "MXE,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"V2H"}->{"RI"})){
		$V2H .= $geneWithAsNumHref->{$geneId}->{"V2H"}->{"RI"} . "RI,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"V2H"}->{"SE"})){
		$V2H .= $geneWithAsNumHref->{$geneId}->{"V2H"}->{"SE"} . "SE,"
	}
	if($V2H eq ""){
		$V2H = "-";
	}else{
		 $V2H = substr($V2H, 0, length($V2H) - 1);
	}

	$V10H = "";
	if(exists($geneWithAsNumHref->{$geneId}->{"V10H"}->{"A3SS"})){
		$V10H .= $geneWithAsNumHref->{$geneId}->{"V10H"}->{"A3SS"} . "A3SS,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"V10H"}->{"A5SS"})){
		$V10H .= $geneWithAsNumHref->{$geneId}->{"V10H"}->{"A5SS"} . "A5SS,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"V10H"}->{"MXE"})){
		$V10H .= $geneWithAsNumHref->{$geneId}->{"V10H"}->{"MXE"} . "MXE,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"V10H"}->{"RI"})){
		$V10H .= $geneWithAsNumHref->{$geneId}->{"V10H"}->{"RI"} . "RI,"
	}
	if(exists($geneWithAsNumHref->{$geneId}->{"V10H"}->{"SE"})){
		$V10H .= $geneWithAsNumHref->{$geneId}->{"V10H"}->{"SE"} . "SE,"
	}
	if($V10H eq ""){
		$V10H = "-";
	}else{
		 $V10H = substr($V10H, 0, length($V10H) - 1);
	}

	print WW join("\t", $geneId, $family, $T2H, $T10H, $V2H, $V10H) . "\n";
}
close FF;
close WW;
