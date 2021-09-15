#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--jgiAnno \\\n" .
		"--protFasta \\\n" .
		"--outputGeneToProtIdList \\\n" .
		"--outputGeneGOanno \n";
	exit;
}

my ($jgiAnno, $protFasta, $outputGeneGOanno, $outputGeneToProtIdList);

GetOptions(
        'jgiAnno=s'=>\$jgiAnno,
	'protFasta=s'=>\$protFasta,
	'outputGeneToProtIdList=s'=>\$outputGeneToProtIdList,
        'outputGeneGOanno=s'=>\$outputGeneGOanno,
);

my ($line, @fieldName, @fieldValue, $i, %hs);
my (%geneToGOList, $geneId, $trsptId, $goList, @go);
open FF, "<$jgiAnno";
# #pacId  locusName       transcriptName  peptideName     Pfam    Panther KOG     ec      KO      GO      Best-hit-arabi-name     arabi-symbol    arabi-defline
$line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
while($line = <FF>){
# 42423414 Gohir.1Z047500 Gohir.1Z047500.1 Gohir.1Z047500.1.p PF00582,PF00069 PTHR27001,PTHR27001:SF116       KOG1187 2.7.11.1	GO:0006950,GO:0006468,GO:0005524,GO:0004672     AT5G57670.2             Protein kinase superfamily protein
	chomp($line);
	@fieldValue = split(/\t/, $line);
	for($i=0; $i<=$#fieldValue; $i++){
		$hs{$fieldName[$i]} = $fieldValue[$i];
	}
	$geneId = $hs{"locusName"};
	$trsptId = $hs{"transcriptName"};
	$goList = $hs{"GO"};

	next if($goList eq "");

	if(not exists($geneToGOList{$geneId})){
		$geneToGOList{$geneId} = $goList;
	}else{
		@go = ();
		@go = split(/\t/, $goList);
		foreach my $go(@go){
			if(index($go, $geneToGOList{$geneId})<0){
				$geneToGOList{$geneId} .= "," . $go;
			}
		}
	}
}
close FF;

my @geneId = keys(%geneToGOList);
open WW, ">$outputGeneGOanno";
print WW "geneId\tGO\n";
foreach my $geneId(@geneId){
	print WW join("\t", $geneId, $geneToGOList{$geneId}) . "\n";
}
close WW;

my (%geneToProtIdList, $geneId, $pacid, $trsptId, $protId);
open FF, "<$protFasta";
# >Gohir.1Z049317.1.p pacid=42364406 transcript=Gohir.1Z049317.1 locus=Gohir.1Z049317 ID=Gohir.1Z049317.1.v2.1 annot-version=v2.1
while($line=<FF>){
	next if(not($line=~/>/));
	($geneId, $protId) = ("", "");
	@fieldValue = split(/ /, $line);
	foreach my $field(@fieldValue){
		if($field=~/>(.*)/){
			$protId = $1;
		}
		if($field=~/locus=(.*)/){
			$geneId = $1;
		}
	}
	if($geneId ne "" and $protId ne ""){
		if(not exists($geneToProtIdList{$geneId})){
			$geneToProtIdList{$geneId} = $protId;
		}else{
			$geneToProtIdList{$geneId} .= "," . $protId;
		}
	}
}
close FF;

open WW, ">$outputGeneToProtIdList";
print WW join("\t", "geneId", "protIdList");
@geneId = ();
@geneId = keys(%geneToProtIdList);
foreach $geneId(@geneId){
	print WW join("\t", $geneId, $geneToProtIdList{$geneId}) . "\n";
}
close WW;
