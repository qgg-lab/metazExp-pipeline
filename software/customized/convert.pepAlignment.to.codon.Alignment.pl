#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--perl \\\n" .
                "--pal2nal \\\n" .
		"--orthoGroupListFile \\\n" .
                "--cdsSeqDir \\\n" .
                "--pepAlignmentDir \\\n" .
                "--codonAlignmentDir \n";
	exit;
}

my ($perl, $pal2nal, $orthoGroupListFile, $cdsSeqDir, $pepAlignmentDir, $codonAlignmentDir);

GetOptions(
        'perl=s'=>\$perl,
        'pal2nal=s'=>\$pal2nal,
        'orthoGroupListFile=s'=>\$orthoGroupListFile,
        'cdsSeqDir=s'=>\$cdsSeqDir,
        'pepAlignmentDir=s'=>\$pepAlignmentDir,
	'codonAlignmentDir=s'=>\$codonAlignmentDir,
);

# 读取orthoGroup，获得所有的orthoGroup Id
my ($cmd, $ogIdList, @ogId,  $ogId, $seqIdList, @seqId, $seqId, $cdsFile, $pepAlignmentFile, $codonAlignmentFile, $errorFile);
$cmd = "grep -v \"Orthogroup\" " . $orthoGroupListFile . " | awk -F \'\\t\' \'{print \$1}\'";
$ogIdList = `$cmd`;
@ogId = split(/\n/, $ogIdList);
foreach $ogId(@ogId){
        $cdsFile = $cdsSeqDir . "/" . $ogId . ".cds.faa";
	$pepAlignmentFile = $pepAlignmentDir . "/" . $ogId . ".align";
	$codonAlignmentFile = $codonAlignmentDir . "/" . $ogId . ".codon.align";
	$errorFile = $codonAlignmentDir . "/log.e." . $ogId;
	$cmd = join(" ", $perl, $pal2nal, $pepAlignmentFile, $cdsFile, "-output fasta", "1>" . $codonAlignmentFile, "2>" . $errorFile);
	#print $cmd;
	#<STDIN>;
	system($cmd);
	system("rm -rf " . $errorFile);
}

