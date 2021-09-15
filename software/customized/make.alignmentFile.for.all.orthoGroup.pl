#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--orthoGroupListFile \\\n" .
                "--orthoGroupSeqPath \\\n" .
                "--mafft \\\n" .
		"--alignmentDir \n";
	exit;
}

my ($orthoGroupListFile, $orthoGroupSeqPath, $mafft, $alignmentDir);

GetOptions(
        'orthoGroupListFile=s'=>\$orthoGroupListFile,
        'orthoGroupSeqPath=s'=>\$orthoGroupSeqPath,
        'mafft=s'=>\$mafft,
        'alignmentDir=s'=>\$alignmentDir,
);

my ($cmd);
my ($orthoGroupIdList, @orthoGroupId, $orthoGroupId);
# 依次读取orthoGroup的编号
$cmd = "grep -v \"Orthogroup\"  $orthoGroupListFile | awk -F \'\\t\' \'{print \$1}\'";
$orthoGroupIdList = `$cmd`;
@orthoGroupId = split(/\n/, $orthoGroupIdList);
foreach $orthoGroupId(@orthoGroupId){
	
	$cmd = $mafft . " --thread 12 --quiet " . $orthoGroupSeqPath . "/" . $orthoGroupId . ".fa > " . $alignmentDir . "/" . $orthoGroupId . ".align";
	system($cmd);
}
