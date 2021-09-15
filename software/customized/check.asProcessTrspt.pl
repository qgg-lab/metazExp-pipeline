#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--cDNAfasta \\\n" .
                "--asProcessTrsptTsv \\\n" .
		"--outputCheckRlt \n";
	exit;
}

my ($origCdnaFasta, $asProcessTrsptTsv, $outputCheckRlt);

GetOptions(
        'origCdnaFasta=s'=>\$origCdnaFasta,
        'asProcessTrsptTsv=s'=>\$asProcessTrsptTsv,
        'outputCheckRlt=s'=>\$outputCheckRlt,
);

# 首先将cDNAfasta序列读入hash
my ($line, @seqId, $seqId);
my (%cDNAseq);
open FF, "<$origCdnaFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@seqId = split(/ /, $1);
		$seqId = $seqId[0];
	}else{
		$cDNAseq{$seqId}.=$line;
	}
}
close FF;

# 读取asProcess
my (@nameField, @valueField, $nameFieldString, $valueFieldString);
my (%tmpHash, $i, $trsptId, $newCdnaSeq, $preSubCdnaSeq, $tailSubCdnaSeq);
open WW, ">$outputCheckRlt";
open FF, "<$asProcessTrsptTsv";
while($line=<FF>){
	chomp($line);
	($nameFieldString, $valueFieldString) = split(/___/, $line);
	@nameField = split(/, /, $nameFieldString);
	@valueField = split(/, /, $valueFieldString);
	for($i=0; $i<=$#nameField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}
	

	$trsptId = $tmpHash{"trsptId"};
	$preSubCdnaSeq = substr($cDNAseq{$trsptId}, 0, $tmpHash{"editSiteInOrigCdna"}); 
	$tailSubCdnaSeq = substr($cDNAseq{$trsptId}, $tmpHash{"editSiteInOrigCdna"} + $tmpHash{"cutSize"});
	if($tmpHash{"residentType"} eq "inclusion"){
		$newCdnaSeq = $preSubCdnaSeq . $tailSubCdnaSeq;
	}elsif($tmpHash{"residentType"} eq "exclusion"){
		$newCdnaSeq = $preSubCdnaSeq . $tmpHash{"insertSeq"} . $tailSubCdnaSeq;
	}else{
		$newCdnaSeq = $preSubCdnaSeq . $tmpHash{"insertSeq"} . $tailSubCdnaSeq;
	}

	if($newCdnaSeq ne $tmpHash{"newCdnaSeq"}){
		print WW join("\t", "asId:" . $tmpHash{"asId"}, "asType:" . $tmpHash{"asType"}, "resideType:" . $tmpHash{"residentType"}, "trsptId:" . $tmpHash{"trsptId"}, "editSite:" . $tmpHash{"editSiteInOrigCdna"}, "cutSize:" . $tmpHash{"cutSize"}, "insertSeq:" . $tmpHash{"insertSeq"}) ."\n";
		print WW "orignal trspt seq:\n";
		print WW $cDNAseq{$trsptId} . "\n";
		print WW "calculated new seq:\n";
		print WW $newCdnaSeq . "\n";
		print WW "registed new seq:\n";
		print WW $tmpHash{"newCdnaSeq"} . "\n";
	}
}
close FF;
close WW;
