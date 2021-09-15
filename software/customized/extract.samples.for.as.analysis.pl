#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
        "--inputTaxonInList 9031,9796,9823,9913,9940\\\n" .
        "--inputTotalSampleFile livestock.samples.tsv\\\n" .
        "--tissue blood\\\n" .
        "--layout RF,FR,UN,R,F,U \\\n" .
        "--phredScore 33,64\\\n" .
        "--minAlignPercent 50\\\n" .
        "--maxAlignPercent 100\\\n" .
        "--minReadLength 70\\\n" .
        "--maxReadLength 150\\\n" .
        "--minMappedMillionSpots 30\\\n" .
        "--maxMappedMillionSpots 100\\\n" .
        "--outputSelectedSampleList selected.experiment.list.tsv \\\n" .
        "--outputSelectedSampleSummary selected.experiment.summary.tsv \n";
	
	exit;
}

my ($inputTaxonInList, $inputTotalSampleFile, $tissue, $layout, $phredScore, $minAlignPercent, $maxAlignPercent, $library);
my ($minReadLength, $maxReadLength, $minMappedMillionSpots, $maxMappedMillionSpots, $outputSelectedSampleList, $outputSelectedSampleSummary);

GetOptions(
	'inputTaxonInList=s'=>\$inputTaxonInList, 
	'inputTotalSampleFile=s'=>\$inputTotalSampleFile, 
	'tissue=s'=>\$tissue, 
	'layout=s'=>\$layout, 
	'library=s'=>\$library,
	'phredScore=s'=>\$phredScore, 
	'minAlignPercent=s'=>\$minAlignPercent, 
	'maxAlignPercent=s'=>\$maxAlignPercent,
	'minReadLength=s'=>\$minReadLength, 
	'maxReadLength=s'=>\$maxReadLength, 
	'minMappedMillionSpots=s'=>\$minMappedMillionSpots, 
	'maxMappedMillionSpots=s'=>\$maxMappedMillionSpots, 
	'outputSelectedSampleList=s'=>\$outputSelectedSampleList, 
	'outputSelectedSampleSummary=s'=>\$outputSelectedSampleSummary
);

my (%summary, %sample);
my (@fieldTitle, @field, $fieldTitle, $field, $percentage);
my (@readLength, $readLength, $mappedSpot);
my ($line, $i);

open WW, ">$outputSelectedSampleList";
open FF, "<$inputTotalSampleFile";
$line = <FF>;
chomp($line);
@fieldTitle = split(/\t/, $line);
while($line = <FF>){
	chomp($line);
	@field = split(/\t/, $line);

	%sample = ();
	for($i=0; $i<=$#field; $i++){
		$sample{$fieldTitle[$i]} = $field[$i];
	}
	
	$mappedSpot = 0;
	if($sample{"alignPer(%)"}=~/(.*)%/){
		$percentage=$1;
	}
	$mappedSpot = $sample{"spotNum(M)"} * $percentage / 100;
	
	if($sample{"readLength"}=~/,/){
		my @readLength = ();
		@readLength = split(/,/, $sample{"readLength"});
		$readLength = $readLength[0];
	}else{
		$readLength = $sample{"readLength"};
	}

	if(	
		uc($tissue) eq uc($sample{"Tissue"}) and 
		&checkTaxonId($inputTaxonInList, $sample{"Taxon"}) == 1 and 
		&checkLayout(uc($layout), uc($sample{"layout"})) == 1 and 
		&checkLibrary(uc($library), uc($sample{"library"})) == 1 and 
		&checkPhredScore(uc($phredScore), uc($sample{"phredScore"})) == 1 and 
		$readLength >= $minReadLength and $readLength <= $maxReadLength and 
		$mappedSpot <= $maxMappedMillionSpots and $mappedSpot >= $minMappedMillionSpots
){
		print WW join("\t", $sample{"Taxon"}, $sample{"PrjId"}, $sample{"expId"}, $sample{"phredScore"}, $sample{"library"}, $sample{"layout"}, $sample{"readLength"}, $sample{"spotNum(M)"}, $sample{"alignPer(%)"}, $mappedSpot) . "\n";
		${$summary{$sample{"Taxon"}}}{"totalNum"}++;
		${$summary{$sample{"Taxon"}}}{$sample{"PrjId"}}++;
	}	

}
close FF;
close WW;
my(@taxonId, @prjId, $prjId);
@taxonId = split(/,/, $inputTaxonInList);
open WW, ">$outputSelectedSampleSummary";
foreach my $taxonId(@taxonId){
	print WW $taxonId;
	if(exists(${$summary{$taxonId}}{"totalNum"})){
		print WW "\t" . "Total:" . ${$summary{$taxonId}}{"totalNum"};
		delete(${$summary{$taxonId}}{"totalNum"});
		@prjId = ();
		@prjId = keys(%{$summary{$taxonId}});
		foreach $prjId(@prjId){
			print WW "\t" . $prjId . ":" . ${$summary{$taxonId}}{$prjId};
		}
		print WW "\n";
	}else{
		print WW "\tTotal:0\n";
	}

}
close WW;
sub checkLayout{
	my ($specifiedLayoutList, $sampleLayoutList) = @_;
	my (@specifiedLayout, @sampleLayout, $specifiedLayout, $sampleLayout);
	@specifiedLayout = split(/,/, $specifiedLayoutList);
	@sampleLayout = split(/,/, $sampleLayoutList);
	foreach $specifiedLayout(@specifiedLayout){
		foreach $sampleLayout(@sampleLayout){
			return 1 if($specifiedLayout eq $sampleLayout);
		}
	}
	return 0;
}

sub checkLibrary{
	my ($specifiedLibraryList, $sampleLibraryList) = @_;
	my (@specifiedLibrary, @sampleLibrary, $specifiedLibrary, $sampleLibrary);
	@specifiedLibrary = split(/,/, $specifiedLibraryList);
	@sampleLibrary = split(/,/, $sampleLibraryList);
	foreach $specifiedLibrary(@specifiedLibrary){
		foreach $sampleLibrary(@sampleLibrary){
			return 1 if($specifiedLibrary eq $sampleLibrary);
		}
	}
	return 0;
}

sub checkTaxonId{
	my ($specifiedTaxonIdList, $sampleTaxonId) = @_;
	my (@specifiedTaxonId, $specifiedTaxonId);
	@specifiedTaxonId = split(/,/, $specifiedTaxonIdList);

	foreach $specifiedTaxonId(@specifiedTaxonId){
		return 1 if($specifiedTaxonId eq $sampleTaxonId);
	}
	return 0;
}

sub checkPhredScore{
	my ($specifiedPhredScoreList, $samplePhredScoreList) = @_;
	my (@specifiedPhredScore, @samplePhredScore, $specifiedPhredScore, $samplePhredScore);
	@specifiedPhredScore = split(/,/, $specifiedPhredScoreList);
	@samplePhredScore = split(/,/, $samplePhredScoreList);
	foreach $specifiedPhredScore(@specifiedPhredScore){
		foreach $samplePhredScore(@samplePhredScore){
			return 1 if($specifiedPhredScore eq $samplePhredScore);
		}
	}
	return 0;
}

