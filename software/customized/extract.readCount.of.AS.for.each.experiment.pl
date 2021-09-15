#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--totalJCECFileList \\\n" .
                "--totalJCFileList \\\n" .
                "--exptIdListFile \\\n" .
		"--outputDir \n";
	exit;
}

my ($totalJCECFileList, $totalJCFileList, $exptIdListFile, $outputDir);

GetOptions(
        'totalJCECFileList=s'=>\$totalJCECFileList,
        'totalJCFileList=s'=>\$totalJCFileList,
        'exptIdListFile=s'=>\$exptIdListFile,
        'outputDir=s'=>\$outputDir,
);

#system("mkdir -p $outputDir/AS/js");
system("mkdir -p $outputDir/AS/jcec");
my (@exptId, $exptId, $line, $psi);
my ($ASID, $IJC_SAMPLE_1, $SJC_SAMPLE_1, $IncFormLen, $SkipFormLen, $Experiment);
my (@jcecFile, $jcecFile);
my (%jcec);
@jcecFile = split(/,/, $totalJCECFileList);
foreach $jcecFile(@jcecFile){
	open FF, "<$jcecFile";
	<FF>;
	while($line=<FF>){
		chomp($line);
		($ASID, $IJC_SAMPLE_1, $SJC_SAMPLE_1, $IncFormLen, $SkipFormLen, $Experiment) = split(/\t/, $line);

		next if($IJC_SAMPLE_1 + $SJC_SAMPLE_1 == 0);

		$psi = sprintf("%.4f", ($IJC_SAMPLE_1/$IncFormLen)/(($IJC_SAMPLE_1/$IncFormLen)+($SJC_SAMPLE_1/$SkipFormLen)));
		$jcec{$Experiment} .= join("\t", $ASID, $IJC_SAMPLE_1, $SJC_SAMPLE_1, $IncFormLen, $SkipFormLen, $psi) . "\n";
	}
	close FF;
}

# 以exptId为单位输出
open FF, "<$exptIdListFile";
@exptId = <FF>;
close FF;

foreach $exptId(@exptId){
	chomp($exptId);
	open WW, ">$outputDir/AS/jcec/$exptId";
	print WW join("\t", "ASID", "IJC_SAMPLE", "SJC_SAMPLE", "IncFormLen", "SkipFormLen", "psi") . "\n";
	if(exists($jcec{$exptId})){
		print WW $jcec{$exptId};
	}
	close WW;
}

