#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--fileListOfExonPosInCodonAlign \\\n" .
                "--taxonIdList \\\n" .
                "--listOfSummarizedRlt \\\n" .
                "--statisticOfSummarizedRlt \\\n" .
                "--outputMappingUniformedIdandOrigId \\\n" .
                "--outputListOfSummarizedRlt \\\n" .
                "--outputFileListOfExonPosInCodonAlign \\\n" .
		"--outputstatisticOfSummarizedRlt \n";
	exit;
}

my ($fileListOfExonPosInCodonAlign, $taxonIdList, $listOfSummarizedRlt, $statisticOfSummarizedRlt, $outputMappingUniformedIdandOrigId, $outputListOfSummarizedRlt, $outputFileListOfExonPosInCodonAlign, $outputstatisticOfSummarizedRlt);

GetOptions(
        'fileListOfExonPosInCodonAlign=s'=>\$fileListOfExonPosInCodonAlign,
        'taxonIdList=s'=>\$taxonIdList,
        'listOfSummarizedRlt=s'=>\$listOfSummarizedRlt,
        'statisticOfSummarizedRlt=s'=>\$statisticOfSummarizedRlt,
	'outputMappingUniformedIdandOrigId=s'=>\$outputMappingUniformedIdandOrigId,
	'outputListOfSummarizedRlt=s'=>\$outputListOfSummarizedRlt,
	'outputFileListOfExonPosInCodonAlign=s'=>\$outputFileListOfExonPosInCodonAlign,
	'outputstatisticOfSummarizedRlt=s'=>\$outputstatisticOfSummarizedRlt,
);

# 读取listOfSummarizedRlt中的orthoCdsAlignPos，对其统一编号
my (%orthoCdsAlignPos, $number);
open FF, "<$listOfSummarizedRlt";
<FF>;
while(my $line=<FF>){
	chomp($line);
	$number++;
	my @field = split(/\t/, $line);
	$orthoCdsAlignPos{$field[0]} = "orthoPos" . sprintf("%06d", $number);
}
close FF;

# 输出新旧id映射文件
open WW, ">$outputMappingUniformedIdandOrigId";
my (@orthoCdsAlignPos, $orthoCdsAlignPos);
@orthoCdsAlignPos = keys(%orthoCdsAlignPos);
foreach $orthoCdsAlignPos(@orthoCdsAlignPos){
	print WW $orthoCdsAlignPos . "\t" . $orthoCdsAlignPos{$orthoCdsAlignPos} . "\n"; 
}
close WW;

# 重新输出list
open FF, "<$listOfSummarizedRlt";
open WW, ">$outputListOfSummarizedRlt";
my $line=<FF>;
print WW $line;
while($line=<FF>){
	my @field = split(/\t/, $line);
	$orthoCdsAlignPos = shift(@field);
	print WW join("\t", $orthoCdsAlignPos{$orthoCdsAlignPos}, @field);
}
close WW;
close FF;

# 重新输出statistic
open FF, "<$statisticOfSummarizedRlt";
open WW, ">$outputstatisticOfSummarizedRlt";
my $line=<FF>;
print WW $line;
while($line=<FF>){
	my @field = split(/\t/, $line);
	$orthoCdsAlignPos = shift(@field);
	print WW join("\t", $orthoCdsAlignPos{$orthoCdsAlignPos}, @field);
}
close WW;
close FF;

# 重新输出各个物种
my (@fileOfExonPosInCodonAlign, $fileOfExonPosInCodonAlign, @outputFileOfExonPosInCodonAlign, $outputFileOfExonPosInCodonAlign);
my (@field, $geneId, $orthoCdsAlignPos);
@fileOfExonPosInCodonAlign=split(/,/, $fileListOfExonPosInCodonAlign);
@outputFileOfExonPosInCodonAlign=split(/,/, $outputFileListOfExonPosInCodonAlign);
for(my $i=0; $i<=$#fileOfExonPosInCodonAlign; $i++){
	$fileOfExonPosInCodonAlign = $fileOfExonPosInCodonAlign[$i];
	$outputFileOfExonPosInCodonAlign = $outputFileOfExonPosInCodonAlign[$i];
	open FF, "<$fileOfExonPosInCodonAlign";
	open WW, ">$outputFileOfExonPosInCodonAlign";
	# AT1G01020       AthaExon000001  1       7157    7232    -       orthoPos108404
	print WW join("\t", "geneId", "exonId", "chr", "exonStart", "exonStop", "exonStrand", "orthoCdsAlignPos") . "\n";
	while($line=<FF>){
		chomp($line);
		@field = split(/\t/, $line);
		my $geneId_orthoCdsAlignPos = pop(@field);
		if($geneId_orthoCdsAlignPos=~/(.*?):(OG\d+:\d+\-\d+)/){
			$geneId = $1;
			$orthoCdsAlignPos = $2;
		}
		print WW join("\t", $geneId, @field, $orthoCdsAlignPos{$orthoCdsAlignPos}) . "\n";
	}
	close FF;
	close WW;
}
