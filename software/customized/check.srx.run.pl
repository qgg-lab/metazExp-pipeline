#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--finishedExptTsv \\\n" .
                "--newExptTsv \\\n" .
		"--reportTsv \n";
	exit;
}

my ($finishedExptTsv, $newExptTsv, $reportTsv);

GetOptions(
        'finishedExptTsv=s'=>\$finishedExptTsv,
        'newExptTsv=s'=>\$newExptTsv,
        'reportTsv=s'=>\$reportTsv,
);
my (%finishedExptToRunList, %newExptToRunList);
my (@fieldName, @fieldValue, %tmpHash, $line, $i, $line);
open FF, "<$finishedExptTsv";
$line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
while($line = <FF>){
	chomp($line);
	@fieldValue = ();
	@fieldValue = split(/\t/, $line);
	for($i=0; $i<=$#fieldName; $i++){
		$tmpHash{$fieldName[$i]} = $fieldValue[$i];
	}
	# 对RunList排序
	my @tmpRun = ();
	@tmpRun = split(/,/, $tmpHash{"RunList"});
	@tmpRun = sort{$a cmp $b}@tmpRun;	
	$finishedExptToRunList{$tmpHash{"Experiment"}} = join(",", @tmpRun);
}

open FF, "<$newExptTsv";
$line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
while($line = <FF>){
	chomp($line);
	@fieldValue = ();
	@fieldValue = split(/\t/, $line);
	for($i=0; $i<=$#fieldName; $i++){
		$fieldValue[$i] = substr($fieldValue[$i], 1, length($fieldValue[$i]) - 2);
		$tmpHash{$fieldName[$i]} = $fieldValue[$i];
	}
	my @tmpRun = ();
	@tmpRun = split(/,/, $tmpHash{"runList"});
	@tmpRun = sort{$a cmp $b}@tmpRun;
	$newExptToRunList{$tmpHash{"experiment_accession"}} = join(",", @tmpRun);
}

# 提取已经完成的expeirment编号
my @exptId = keys(%finishedExptToRunList);
foreach my $exptId(@exptId){
	if(not exists($newExptToRunList{$exptId}) or $finishedExptToRunList{$exptId} ne $newExptToRunList{$exptId}){
		#print join("\t", $exptId, "finished:", $finishedExptToRunList{$exptId}, "new:", $newExptToRunList{$exptId}) . "\n";
		print $exptId . "\n";
	}
}

