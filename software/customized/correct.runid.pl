#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--runInfoList \\\n" .
                "--exptInfoTsv \\\n" .
		"--outputUpdatedExpts \n";
	exit;
}

my ($runInfoList, $exptInfoTsv, $outputUpdatedExpts);

GetOptions(
        'runInfoList=s'=>\$runInfoList,
        'exptInfoTsv=s'=>\$exptInfoTsv,
        'outputUpdatedExpts=s'=>\$outputUpdatedExpts,
);

# 读取更新的expt->runling
my (%exptToRunList, $line, @nameField, @valueField, $field, %tmpHash, $i, $runId, $srxId);
open FF, "<$runInfoList";
$line = <FF>;
chomp($line);
@nameField = split(/,/, $line);
while($line = <FF>){
	chomp($line);
	@valueField = split(/,/, $line);
	$runId = $valueField[0];
	$srxId = "";
	for($i=0; $i<=$#valueField; $i++){
		if($valueField[$i]=~/^SRX\d+/){
			$srxId = $valueField[$i];
		}
	}
	$exptToRunList{$srxId} .= $runId . ",";
#	print $srxId . ":" . $exptToRunList{$srxId};
#	<STDIN>;
}
close FF;

# 将新的runlist更新到输出文件中
my ($outputLine);
open FF, "<$exptInfoTsv";
open WW, ">$outputUpdatedExpts";
$line = <FF>;
print WW $line;
chomp($line);
@nameField = split(/\t/, $line);
while($line = <FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	$srxId = $tmpHash{"Experiment"};

	if(exists($exptToRunList{$srxId})){
		$tmpHash{"RunList"} = substr($exptToRunList{$srxId}, 0, length($exptToRunList{$srxId})-1);
	}
	$outputLine = "";
	for($i=0; $i<=$#nameField; $i++){
		$outputLine .= $tmpHash{$nameField[$i]} . "\t";
	}
	# 
	$outputLine = substr($outputLine, 0, length($outputLine) - 1);
	print WW $outputLine . "\n";
}
close FF;
close WW;
