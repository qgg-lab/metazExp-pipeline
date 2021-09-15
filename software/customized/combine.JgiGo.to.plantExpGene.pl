#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--jgiGeneToGoList \\\n" .
                "--jgiGeneToProtIdList \\\n" .
                "--plantExpGeneBlastToJgiProtId \\\n" .
		"--outputPlantGeneToGo \n";
	exit;
}

my ($jgiGeneToGoList, $jgiGeneToProtIdList, $plantExpGeneBlastToJgiProtId, $outputPlantGeneToGo);

GetOptions(
        'jgiGeneToGoList=s'=>\$jgiGeneToGoList,
        'jgiGeneToProtIdList=s'=>\$jgiGeneToProtIdList,
        'plantExpGeneBlastToJgiProtId=s'=>\$plantExpGeneBlastToJgiProtId,
        'outputPlantGeneToGo=s'=>\$outputPlantGeneToGo,
);

my (@field, $line, $i, $field);
my (%jgiGeneToGoList, %jgiProtToJgiGene);

# 将jgiGene的GoList读入hash
open FF, "<$jgiGeneToGoList";
# geneId  GO
# Gohir.A08G076350        GO:0006412,GO:0005840,GO:0005622,GO:0003735
<FF>;
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	$jgiGeneToGoList{$field[0]} = $field[1];
}
close FF;

# jgi的prot到gene之间映射读入
my ($jgiGeneId, $protIdList, @protId, %jgiProtIdToJgiGeneId);
open FF, "<$jgiGeneToProtIdList";
<FF>;
# Gohir.A07G225800        Gohir.A07G225800.4.p,Gohir.A07G225800.1.p
while($line=<FF>){
	chomp($line);
	($jgiGeneId, $protIdList) = split(/\t/, $line);
	@protId = ();
	@protId = split(/,/, $protIdList);
	foreach my $protId(@protId){
		$jgiProtIdToJgiGeneId{$protId} = $jgiGeneId;
	} 
}
close FF;

# 读取从plantExp geneId到jgiProtId的映射关系，然后获得jgiGeneId，接着获得GOList输出
my (%plantExpGeneIdToJgiProt, $hrf, @geneId, $geneId);
$hrf = \%plantExpGeneIdToJgiProt;
open FF, "<$plantExpGeneBlastToJgiProtId";
# LOC107939291    Gohir.D13G232600.2.p    100.00  36      0       0       34      69      1       36      7e-19   83.2
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$geneId = $field[0];
	if(not exists($plantExpGeneIdToJgiProt{$geneId})){
		$hrf->{$geneId}->{"jgiProt"} = $field[1];
		$hrf->{$geneId}->{"pvalue"} = $field[10];
	}else{
		if($field[10] < $hrf->{$geneId}->{"pvalue"}){
			$hrf->{$geneId}->{"jgiProt"} = $field[1];
			$hrf->{$geneId}->{"pvalue"} = $field[10];
		}
	}
}
close FF;

# 输出最后plantExp gene的GO注释
my ($jgiGeneId);
open WW, ">$outputPlantGeneToGo";
@geneId = ();
@geneId = keys(%plantExpGeneIdToJgiProt);
foreach $geneId(@geneId){
	$jgiGeneId = $jgiProtIdToJgiGeneId{$hrf->{$geneId}->{"jgiProt"}};
	if(exists($jgiGeneToGoList{$jgiGeneId})){
		print WW join("\t", $geneId, $jgiGeneToGoList{$jgiGeneId}) . "\n";
	}
}
close WW;
