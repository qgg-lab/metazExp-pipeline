#!/usr/bin/perl
use strict;
my ($gtfFileList, $gtfTypeList, $combinedGtfFile) = @ARGV;

my (@gtfFile, $gtfFile, @gtfType, $gtfType, $i);

@gtfFile = split(/,/, $gtfFileList);
@gtfType = split(/,/, $gtfTypeList);

my (%gtf, $gtfHref);
$gtfHref=\%gtf;
my (@geneId, $geneId, @trsptId, $trsptId, $line, @field);
# 将geneGtf读入
for($i=0; $i<=$#gtfFile; $i++){
	$gtfFile = $gtfFile[$i];
	$gtfType = $gtfType[$i];

	open FF, "<$gtfFile";
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		($geneId, $trsptId) = ("", "");
		&getGeneIdAndTrsptId($field[8], \$geneId, \$trsptId);
		if($gtfType eq "gene"){
			$gtfHref->{$geneId}->{"geneFeature"} = $line;
			$geneId[$#geneId+1] = $geneId;
		}else{
			$gtfHref->{$geneId}->{"trspt"}->{$trsptId} .= $line;
		}
	}
	close FF;
}

open WW, ">$combinedGtfFile";
foreach $geneId(@geneId){
	print WW $gtfHref->{$geneId}->{"geneFeature"};
	@trsptId = ();
	@trsptId = keys(%{$gtfHref->{$geneId}->{"trspt"}});
	foreach $trsptId(@trsptId){
		print WW $gtfHref->{$geneId}->{"trspt"}->{$trsptId};
	}
}
close WW;


sub getGeneIdAndTrsptId{
	my ($attrString, $geneId, $trsptId) = @_;
	my (@attr, $attr);
	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$geneId = $1;
		}	
		if($attr=~/transcript_id "(.*)"/){
			$$trsptId = $1;
		}
	}
}
