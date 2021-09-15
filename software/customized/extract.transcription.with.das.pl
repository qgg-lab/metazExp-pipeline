#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--tfdbFamily tf.family.txt\\\n" .
                "--tfdbTrspt2Ourtrspt tfdbgeneid.to.trsptid.mapping \\\n" .
                "--ourTrspt2OurGene tfdbgeneid.to.trsptid.mapping \\\n" .
		"--T2H T_DAS_2h.vs.0h.tsv \\\n" .
		"--T10H T_DAS_10h.vs.0h.tsv \\\n" .
		"--V2H V_DAS_2h.vs.0h.tsv \\\n" . 
		"--V10H V_DAS_10h.vs.0h.tsv \\\n";
	exit;
}

my ($tfdbFamily, $tfdbTrspt2Ourtrspt, $ourTrspt2OurGene);
my ($T2H, $T10H, $V2H, $V10H);
GetOptions(
        'tfdbFamily=s'=>\$tfdbFamily,
        'tfdbTrspt2Ourtrspt=s'=>\$tfdbTrspt2Ourtrspt,
        'ourTrspt2OurGene=s'=>\$ourTrspt2OurGene,
	'T2H=s'=>\$T2H,
	'T10H=s'=>\$T10H,
	'V2H=s'=>\$V2H,
	'V10H=s'=>\$V10H,
);

my ($line, @tmp);
my (%tfdbTrsptFamily);

my ($href);
$href = \%tfdbTrsptFamily;
open FF, "<$tfdbFamily";
# Species Protein ID      Family  Genome/EST      Category
# Oryza sativa subsp. japonica    LOC_Os02g27060.1        ARID    Genome  Other
# Oryza sativa subsp. japonica    LOC_Os02g27060.2        ARID    Genome  Other
<FF>;
while($line=<FF>){
	chomp($line);
	@tmp = ();
	@tmp = split(/\t/, $line);
	$href->{$tmp[1]}->{"family"} = $tmp[2];
	$href->{$tmp[1]}->{"type"} = $tmp[4];
}
close FF;


my (%trspt2TfDBTrspt, $tfdbTrsptId, @tt);
open FF, "<$tfdbTrspt2Ourtrspt";
# Oryza_sativa_subsp._japonica_LOC_Os01g01290.1   Os01t0102400-01 100.000 421     0       0       1       421     1       421     0.0     874
while($line=<FF>){
	chomp($line);
	@tmp = ();
	@tmp = split(/\t/, $line);
	@tt = ();
	@tt = split(/_japonica_/, $tmp[0]);
	$tfdbTrsptId = $tt[1];
	$trspt2TfDBTrspt{$tmp[1]} = $tfdbTrsptId;
}
close FF;


# 读取T2H, T10H, V2H, V10H中gene发生各种AS类型的数量
my (%geneWithAsNum, $geneWithAsNumHref, @tmp);
my ($asId, $orthoAsId, $control, $treatment, $dltPsi, $PValue, $QValue, $geneId, $orthoGeneId, $symbol, $GO, $Pathway);
$geneWithAsNumHref = \%geneWithAsNum;
open FF, "<$T2H";
# asId orthoAsId control treatment dltPsi P-Value Q-Value geneId orthoGeneId symbol GO Pathway
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

my ($geneId, $trsptId,  $tfdbTrsptId, $family, $type, $T2H, $T10H, $V2H, $V10H);
open FF, "<$ourTrspt2OurGene";
# Os02g0194200    Os02t0194200-01
print join("\t", "geneId", "family", "type", "T2H", "T10H", "V2H", "V10H") . "\n";
while($line=<FF>){
	chomp($line);
	($geneId, $trsptId) = split(/\t/, $line);
	$tfdbTrsptId = $trspt2TfDBTrspt{$trsptId};
	$family = "NA";
	$family = $href->{$tfdbTrsptId}->{"family"} if(exists($href->{$tfdbTrsptId}->{"family"}));
	$type = "NA";
	$type = $href->{$tfdbTrsptId}->{"type"} if(exists($href->{$tfdbTrsptId}->{"type"}));

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

	print join("\t", $geneId, $family, $type, $T2H, $T10H, $V2H, $V10H) . "\n";
}
close FF;
