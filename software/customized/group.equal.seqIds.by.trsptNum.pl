#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--trsptNumInSeqList trsptNum.in.seq.list \\\n" .
		"--gtfFile final.complete.trspt.anno.gtf \\\n" .
                "--groupNum 36 \\\n" .
                "--subGtfPrefix subGtf\\\n" .
		"--outputDir ./ \n";
	exit;
}

my ($trsptNumInSeqList, $gtfFile, $groupNum, $subGtfPrefix, $outputDir);

GetOptions(
        'trsptNumInSeqList=s'=>\$trsptNumInSeqList,
        'gtfFile=s'=>\$gtfFile,
        'groupNum=s'=>\$groupNum,
        'subGtfPrefix=s'=>\$subGtfPrefix,
	'outputDir=s'=>\$outputDir,
);

# 将gtf读入hash
my (%gtf, $line, $seqId, @field);
open FF, "<$gtfFile";
while($line=<FF>){
	next if($line=~/^#/);
	@field = split(/\t/, $line);
	$seqId = $field[0];
	$gtf{$seqId}.=$line;
}
close FF;


# 将trsptNumInSeqList读入数组中
my (@trsptNumInSeq, $line, $trsptNum, $seqId, $seqNum);
open FF, "<$trsptNumInSeqList";
# 45      NC_001879.2
# 27      NC_006581.1
$seqNum = 0;
while($line=<FF>){
	chomp($line);
	($trsptNumInSeq[$seqNum][0], $trsptNumInSeq[$seqNum][1]) = split(/\t/, $line);
	$seqNum++;
}
close FF;

my ($i, $j, $k, $groupSize);
$groupSize = int(($#trsptNumInSeq+1)/$groupNum);
for($i=0; $i<$groupSize; $i++){
	# $i为偶数
	if(int($i/2)*2==$i){
		$k=0;
		for($j=$i*$groupNum; $j<($i+1)*$groupNum; $j++){
			open WW, ">>$outputDir" . "/" . $subGtfPrefix . ".seqIdList." . $k;
			open GTF, ">>$outputDir" . "/" . $subGtfPrefix . ".gtf." . $k;
			#print WW join("\t", $trsptNumInSeq[$j][0], $trsptNumInSeq[$j][1]) . "\n";
			print WW join("\t", $trsptNumInSeq[$j][1]) . "\n";
			print GTF $gtf{$trsptNumInSeq[$j][1]};
			close WW;
			close GTF;
			$k++;
		}
	}else{
	# $i为奇数，反过来
		$k=0;
		for($j=($i+1)*$groupNum-1; $j>=($i*$groupNum); $j--){
			open WW, ">>$outputDir" . "/". $subGtfPrefix . ".seqIdList." . $k;
			open GTF, ">>$outputDir" . "/" . $subGtfPrefix . ".gtf." . $k;
			#print WW join("\t", $trsptNumInSeq[$j][0], $trsptNumInSeq[$j][1]) . "\n";
			print WW join("\t", $trsptNumInSeq[$j][1]) . "\n";
			print GTF $gtf{$trsptNumInSeq[$j][1]};
			close WW;
			close GTF;
			$k++;
		}

	}
}

if($groupSize*$groupNum < $#trsptNumInSeq){
	$k=0;
	for($j=$groupSize*$groupNum; $j<=$#trsptNumInSeq; $j++){
		open WW, ">>$outputDir" . "/" . $subGtfPrefix . ".seqIdList." . $k;
		open GTF, ">>$outputDir" . "/" . $subGtfPrefix . ".gtf." . $k;
		print WW join("\t", $trsptNumInSeq[$j][1]) . "\n";
		print GTF $gtf{$trsptNumInSeq[$j][1]};
		close WW;
		close GTF;
		$k++;
	}
}
