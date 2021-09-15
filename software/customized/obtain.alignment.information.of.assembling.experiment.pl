#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV < 0){
	print "perl $0 \\\n" . 
		"\t\t --assemblyDir ../002-assemble-trsptome-on-goodExps/psiOutputDir \\\n" . 
		"\t\t --curatedExprimentTsv ../001-prepare-local-datasource/curated.experiment.tsv \\\n" . 
		"\t\t --outputExptAlignInfoTsv  ./alignment.info.of.assmbled.experiment.tsv \n";
	exit(0);
}

my ($assemblyDir, $curatedExprimentTsv, $outputExptAlignInfoTsv);

GetOptions(
        'assemblyDir=s'=>\$assemblyDir,
        'curatedExprimentTsv=s'=>\$curatedExprimentTsv,
	'outputExptAlignInfoTsv=s'=>\$outputExptAlignInfoTsv
);
my (@fieldName, @field, $i, $exptId, $exptIdPos);
my ($experimentLine, %experiment);

# 将手工挑选的experiment信息读入到hash中
open FF, "<$curatedExprimentTsv";
open WW, ">$outputExptAlignInfoTsv";
$experimentLine=<FF>;
chomp($experimentLine);
print WW join("\t", $experimentLine, "alignPercent", "mappedBases", "mappedReadNum", "detectedReadLen", "libraryType", "phredScore") . "\n";

@fieldName = split(/\t/, $experimentLine);
my $exptIdPos = 0;
for($exptIdPos=0; $exptIdPos<=$#fieldName; $exptIdPos++){
        last if($fieldName[$exptIdPos] eq "Experiment");
}

my (@tmpArr, $libraryType, $alignPercent, $mappedBaseNum, $mappedReadNum, $detectedReadLen, $detectedPhredScore);
while($experimentLine=<FF>){
        chomp($experimentLine);
        @field = ();
        @field = split(/\t/, $experimentLine);
	$exptId = $field[$exptIdPos];
        for($i=0; $i<=$#field; $i++){
                ${$experiment{$exptId}}{$fieldName[$i]} = $field[$i];
        }
	next if(${$experiment{$exptId}}{"Assemble"} != 1);
	# 没有检测到组装的transcriptome
	if(not(-e $assemblyDir . "/" . $exptId . "/transcriptomeByStringtie.gtf")){
		# print WW "\t-1\t-1\t-1\t-1\t-\t-1\n";
		next;
	}
	
	print WW $experimentLine;
	# 检测到组装的transcriptome，开始检测相关信息
	open SEQINFO, "<$assemblyDir" . "/" . $exptId . "/" . $exptId . ".SeqInfo.txt";
	<SEQINFO>;
	while(my $line = <SEQINFO>){
		chomp($line);
		@tmpArr = ();
		@tmpArr = split(/\t/, $line);
		if($#tmpArr==7){
			$detectedPhredScore = $tmpArr[4];
			$detectedReadLen = $tmpArr[3];
			$libraryType = $tmpArr[6];
		}elsif($line=~/TotalSpot:(.*)\tAlignPercent:(.*)%/){
			$alignPercent = $2;
			$mappedBaseNum = sprintf("%.2f", ${$experiment{$exptId}}{"Base"} * $alignPercent / 100);
			$mappedReadNum = sprintf("%.2f", ${$experiment{$exptId}}{"ReadNum"} * $alignPercent / 100);
		}
	}
	close SEQINFO;

	print WW "\t$alignPercent\t$mappedBaseNum\t$mappedReadNum\t$detectedReadLen\t$libraryType\t$detectedPhredScore\n";
}
close FF;
