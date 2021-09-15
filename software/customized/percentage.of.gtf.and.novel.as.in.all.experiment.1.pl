#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonIdList \\\n" .
                "--baseDir \\\n" .
                "--outputAsTypeOriginPercentTsv \\\n" .
	exit;
}

my ($taxonIdList, $baseDir, $outputAsTypeOriginPercentTsv);

GetOptions(
        'taxonIdList=s'=>\$taxonIdList,
        'baseDir=s'=>\$baseDir,
        'outputAsTypeOriginPercentTsv=s'=>\$outputAsTypeOriginPercentTsv,
);



my (@taxonId, $taxonId, %tmpHash, $line, @field, $exptId);
my ($rltDir, $exptIdListFile);
my ($totalA3ssNum, $totalA5ssNum, $totalMxeNum, $totalRiNum, $totalSeNum);
my ($novelA3ssNum, $novelA5ssNum, $novelMxeNum, $novelRiNum, $novelSeNum);
my ($gtfA3ssNum, $gtfA5ssNum, $gtfMxeNum, $gtfRiNum, $gtfSeNum);
open FF, "<$taxonIdList";
@taxonId = <FF>;
close FF;

#open NUM, ">$outputAsOrigNumTsv";
#print NUM join("\t", "Taxon", "ExptId", "AsNumName", "Num") . "\n";
open PER, ">$outputAsTypeOriginPercentTsv";
print PER join("\t", "Taxon", "ExptId", "total", "A3SS", "A5SS", "MXE", "RI", "SE", "GTF", "Mapping") . "\n";

foreach $taxonId(@taxonId){
	chomp($taxonId);
	my ($rltDir, $exptIdListFile);
	$exptIdListFile = $baseDir . "/" . $taxonId . "/010-gather-alignment-info-of-all-expts/cutoff.info.of.assembled.experiment.tsv";
	open FF, "<$exptIdListFile";
	while($line=<FF>){
		@field = ();
		@field = split(/\//, $line);
		$exptId = $field[$#field-1];
		$rltDir = $baseDir . "/" . $taxonId . "/008-pickup-psi-of-ASs-in-all-expers/psiOutputDir/" . $exptId;

		$novelA3ssNum = `cat $rltDir/fromGTF.novelEvents.A3SS.txt |wc -l`;
		chomp($novelA3ssNum);
		$novelA3ssNum -=1;

		$novelA5ssNum =`cat $rltDir/fromGTF.novelEvents.A5SS.txt |wc -l`;
		chomp($novelA5ssNum);
		$novelA5ssNum -=1;
		
		$novelMxeNum = `cat $rltDir/fromGTF.novelEvents.MXE.txt |wc -l`;
		chomp($novelMxeNum);
		$novelMxeNum -=1;

		$novelRiNum = `cat $rltDir/fromGTF.novelEvents.RI.txt |wc -l`;
		chomp($novelRiNum);
		$novelMxeNum -=1;

		$novelSeNum = `cat $rltDir/fromGTF.novelEvents.SE.txt |wc -l`;
		chomp($novelSeNum);
		$novelSeNum -=1;

		$totalA3ssNum = `cat $rltDir/fromGTF.A3SS.txt |wc -l`;
		chomp($totalA3ssNum);
		$totalA3ssNum -=1;

		$totalA5ssNum =`cat $rltDir/fromGTF.A5SS.txt |wc -l`;
		chomp($totalA5ssNum);
		$totalA5ssNum -=1;
		
		$totalMxeNum = `cat $rltDir/fromGTF.MXE.txt |wc -l`;
		chomp($totalMxeNum);
		$totalMxeNum -=1;

		$totalRiNum = `cat $rltDir/fromGTF.RI.txt |wc -l`;
		chomp($totalRiNum);
		$totalMxeNum -=1;

		$totalSeNum = `cat $rltDir/fromGTF.SE.txt |wc -l`;
		chomp($totalSeNum);
		$totalSeNum -=1;

		$gtfA3ssNum = $totalA3ssNum - $novelA3ssNum;
		$gtfA5ssNum = $totalA5ssNum - $novelA5ssNum;
		$gtfMxeNum = $totalMxeNum - $novelMxeNum;
		$gtfRiNum = $totalRiNum - $novelRiNum;
		$gtfSeNum = $totalSeNum - $novelSeNum;

		my $totalGtfNum = $gtfA3ssNum + $gtfA5ssNum + $gtfMxeNum + $gtfRiNum + $gtfSeNum;
		my $totalNovelNum = $novelA3ssNum + $novelA5ssNum + $novelMxeNum + $novelRiNum + $novelSeNum;
		my $totalNum = $totalA3ssNum + $totalA5ssNum + $totalMxeNum + $totalRiNum + $totalSeNum;

		next if($totalA3ssNum=0 and $totalA5ssNum ==0 and $totalMxeNum==0 and $totalRiNum == 0 and $totalSeNum == 0);

		print PER join("\t", $taxonId, $exptId, $totalNum, $totalA3ssNum/$totalNum, $totalA5ssNum/$totalNum, $totalMxeNum/$totalNum, $totalRiNum/$totalNum, $totalSeNum/$totalNum, $totalGtfNum/$totalNum, $totalNovelNum/$totalNum) . "\n";

	}
	close FF;
}
close PER;
