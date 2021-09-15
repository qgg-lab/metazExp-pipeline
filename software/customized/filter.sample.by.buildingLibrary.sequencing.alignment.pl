#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV < 0){
	print "This script is used to gather psi of AS in experiment and experiment sequencing information\n\n";
	print "\t perl $0 \\\n" . 
		"\t\t --inputSampleListFile sample.Information.output.by.psiScript.tsv \\\n" .
		"\t\t --filterMinReadLen 30 \\\n" .
		"\t\t --filterLibLayout SINGLE,PAIRED \\\n" .
		"\t\t --filterLibType U,R,F,UN,RF,FR \\\n" .
		"\t\t --filterAlignPer 50 \\\n" .
		"\t\t --filterMappedSpots 20 \\\n" .
		"\t\t --outputExperimentInfo experimentInfo.tsv\n"; 
	exit(0);
}

my ($inputSampleListFile, $psiOutputDir, $outputExperimentInfo);
my ($filterMinReadLen, $filterLibLayout, $filterLibType, $filterAlignPer, $filterTotalSpots, $filterMappedSpots);

GetOptions(
	'inputSampleListFile=s'=>\$inputSampleListFile,
	'filterMinReadLen=s'=>\$filterMinReadLen, 
	'filterLibLayout=s'=>\$filterLibLayout,
	'filterLibType=s'=>\$filterLibType,
	'filterAlignPer=s'=>\$filterAlignPer,
	'filterMappedSpots=s'=>\$filterMappedSpots,
        'outputExperimentInfo=s'=>\$outputExperimentInfo,
);

my ($line, @fields, $totalMappedSpots, $alignPercentage, @tmp);

open FF, "<$inputSampleListFile";
open WW, ">$outputExperimentInfo";
while($line=<FF>){

	print WW $line if($line=~/^expId/);
	next if($line=~/\tERROR\t/);
	# 0     1	2	3	4	5	6	7	   8	 	9	10
	# expId status runNum runId library layout phredScore readLength spotNum(M) alignPer(%) novelSpliceNum
	# SRX73  OK      1    SRR16    UN   PAIRED  33          90         24.86      93.27%    97637  
	@fields = ();
	@fields = split(/\t/, $line);
	$alignPercentage = substr($fields[9], 0, length($fields[9]) -1 );
	next if($alignPercentage * $fields[8] / 100 < $filterMappedSpots);
	next if(&checkMinReadLen($fields[7], $filterMinReadLen) == 0);
	next if(&checkLibLayout($fields[5], $filterLibLayout) == 0);
	next if(&checkLibType($fields[4], $filterLibType)==0);
	print WW $line;
	
	
}
close FF;

#obtain the num of AS with at least one read support in IJC or SJC
sub getNumWithReadSupport{
	my $asFile = $_[0];
	my (@tt, $asNum);
	open ASFF, "<$asFile";	
	<ASFF>;
	while(my $line = <ASFF>){
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       4       21                      65      49
		@tt = ();
		@tt = split(/\t/, $line);
		if($tt[1]!=0 or $tt[2]!=0){
			$asNum++;
		}
	}
	close ASFF;
	return $asNum;
}

sub checkMinReadLen{
	# 77,75		70
	my ($readLengthX, $filterMinReadLenX) = @_;
	my (@tt, $i);
	@tt = split(/,/, $readLengthX);
	foreach my $len(@tt){
		return 0 if($len < $filterMinReadLenX);
	}
	return 1;
}

sub checkLibLayout{
	# single,single  single,paired
	my ($layoutX, $filterLibLayoutX) = @_;
	my (@tt, $i);
	@tt = split(/,/, $layoutX);
	foreach my $tmpLayout(@tt){
		return 0 if(index($filterLibLayoutX, $tmpLayout)<0);
	}
	return 1;
}
sub checkLibType{
	# R,R   R,F,RF,FR
	my ($libraryX, $filterLibTypeX) = @_;
	my (@type, @ftype, $ftype);

	@type = split(/,/, $libraryX);
	@ftype = split(/,/, $filterLibTypeX);
	foreach $ftype(@ftype){
		return 1 if($type[0] eq $ftype);
	}
	return 0;
}
