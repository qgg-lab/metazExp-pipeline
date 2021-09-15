#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV < 0){
	print "This script is used to select transcriptome from good experiments.\n\n";
	print "perl $0 \\\n" . 
		"\t\t --assemblyDir                /mnt/home/liujind1/workAS/01-cattle/001-submit-assembly-on-goodExps/assemblyDir     \\\n" . 
		"\t\t --specifiedExpIdList         /mnt/home/liujind1/workAS/01-cattle/003-build-complete-annotation/specified.experimentid.list \\\n" . 
		"\t\t --totalSopts        40 \\\n" .
		"\t\t --alignPer          50 \\\n" .
		"\t\t --mappedSpots       40 \\\n" .
		"\t\t --layout            PAIRED,SINGLE \\\n" .
		"\t\t --minReadLength        50 \\\n" .
		"\t\t --libraryType       RF,FR,UN,R,F,U \\\n" . 
		"\t\t --outputSelectedGtfList /mnt/home/liujind1/workAS/01-cattle/002-build-tripletAnnotation/selectedGtfList.txt \n";
	exit(0);
}

my ($assemblyDir, $specifiedExpIdList, $outputExperimentInfo, $outputSelectedGtfList);
my ($totalSpots, $alignPer, $mappedSpots);
my ($layout, $minReadLength, $outputTripletAnno, $libraryType);
my ($expTranscriptomeFileText);
my ($expId, @expTranscriptomeFile);
my ($totalSpotsInExp, $alignPerInExp);

$specifiedExpIdList = "";

GetOptions(
        'assemblyDir=s'=>\$assemblyDir,
        'specifiedExpIdList=s'=>\$specifiedExpIdList,
	'totalSopts=i'=>\$totalSpots,
	'alignPer=i'=>\$alignPer,
	'mappedSpots=i'=>\$mappedSpots,
	'layout=s'=>\$layout,
	'minReadLength=i'=>\$minReadLength,
	'libraryType=s'=>\$libraryType,
        'outputTripletAnno=s'=>\$outputTripletAnno,
	'outputSelectedGtfList=s'=>\$outputSelectedGtfList,
);

#检测assemblyDir目录是否存在
if(not -e $assemblyDir){
	print STDERR "$assemblyDir doesn't exist!\n";
	exit;
}

#extract experiment Id into array
my $transcriptomeNum = 0;
$assemblyDir = $assemblyDir . "/" if(substr($assemblyDir, length($assemblyDir)-1, 1) ne "/");
if(-s $specifiedExpIdList){
	open FF, "<$specifiedExpIdList";
	while(my $expId=<FF>){
		chomp($expId);
		if(not -s $assemblyDir . $expId . "/transcriptomeByStringtie.gtf"){
			print STDERR $assemblyDir . $expId . " hasn't transcriptome assembly.\n";
		}
	}
	close FF;
}else{
	$expTranscriptomeFileText = `find $assemblyDir -name transcriptomeByStringtie.gtf`;
	@expTranscriptomeFile = split(/\n/, $expTranscriptomeFileText);
	$transcriptomeNum = $#expTranscriptomeFile + 1;
	print STDOUT "A total of $transcriptomeNum transcriptome of good experiments were detected.\n";
}


#gather information from experiment dir
my ($experimentDir);
my (@tmpArr);
my ($mergeCmdString);
my (@libraryType, $libraryTypeTmp, $readLengthTmp, @layout, $layoutTmp, $seqInforFile);

my $selectedExpNum = 0;
open SELECTEDEXPERIMENT, ">$outputSelectedGtfList";
for(my $i=0; $i<=$#expTranscriptomeFile; $i++){

	@tmpArr = ();
	@tmpArr = split(/\//, $expTranscriptomeFile[$i]);
	
	#pop transcriptomeByStringtie.gtf from array terminal
	pop(@tmpArr);

	#get experiment Id
	$expId = $tmpArr[$#tmpArr];
	
	#detect sequencing information
	$seqInforFile = join("/", @tmpArr) . "/" . $expId . ".SeqInfo.txt";

	open FF, "<$seqInforFile";
	<FF>;
	while(my $line = <FF>){
		chomp($line);
		@tmpArr = ();
		@tmpArr = split(/\t/, $line);
		if($#tmpArr==7){
			$layoutTmp = $tmpArr[2];
			$readLengthTmp = $tmpArr[3];
			$libraryTypeTmp = $tmpArr[6];
		}elsif($line=~/TotalSpot:(.*)\tAlignPercent:(.*)%/){
			$totalSpotsInExp = int(($1/1000000)*100)/100;
			$alignPerInExp = $2;
		}
	}
	close FF;

	my @libraryArr = split(/,/, $libraryType);
	my @layoutArr = split(/,/, $layout);

	if(not grep {$_ eq $layoutTmp} @layoutArr or 
		not grep {$_ eq $libraryTypeTmp} @libraryArr or 
		$readLengthTmp <$minReadLength or
		$totalSpotsInExp < $totalSpots or
		$alignPerInExp < $alignPer or
		$totalSpotsInExp * $alignPerInExp / 100 < $mappedSpots
	){
		print STDERR "Experiment " . $expId . " can't meet filter demand.\n";
		next;
	}else{
		print SELECTEDEXPERIMENT join("\t", $expTranscriptomeFile[$i], $layoutTmp, $libraryTypeTmp, $readLengthTmp, $totalSpotsInExp, $alignPerInExp, $totalSpotsInExp * $alignPerInExp / 100 ). "\n";
		$selectedExpNum++;		
	}
}
close SELECTEDEXPERIMENT;
print STDOUT "A total of $selectedExpNum transcriptome will be used to build complete annoation.\n";
