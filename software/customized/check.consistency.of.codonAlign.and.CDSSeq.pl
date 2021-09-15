#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--cdsSeqDir \\\n" .
                "--codonAlignmentDir \\\n" .
                "--orthoGroupListFile \\\n" .
		"--checkRltTsv \n";
	exit;
}

my ($cdsSeqDir, $codonAlignmentDir, $orthoGroupListFile, $checkRltTsv);

GetOptions(
        'cdsSeqDir=s'=>\$cdsSeqDir,
        'codonAlignmentDir=s'=>\$codonAlignmentDir,
        'orthoGroupListFile=s'=>\$orthoGroupListFile,
        'checkRltTsv=s'=>\$checkRltTsv,
);

my (@ogId, $ogIdList, $cmd, $ogId, $seqId, @seqId);
# 将orthoGroup编号读入的数组
$cmd = "grep -v \"Orthogroup\" " . $orthoGroupListFile . " | awk -F \'\\t\' \'{print \$1}\'";
$ogIdList=`$cmd`;
@ogId = split(/\n/, $ogIdList);
print "orthogroup num:" . $#ogId . "\n";

# 将所有CDS序列读入
my (%cdsSeq, $seqId);
foreach $ogId(@ogId){
	open FF, "<$cdsSeqDir" . "/" . $ogId . ".cds.faa";
	while(my $line=<FF>){
		chomp($line);
		if($line=~/>(.*)/){
			$seqId = $1;
		}else{
			$cdsSeq{$seqId}.=$line;
		}
	}
	close FF;
}
print "finish load CDS into hash.\n";


# 将所有codonAlign中的序列读入hash
my (%codonAlignSeq);
foreach $ogId(@ogId){
	open FF, "<$codonAlignmentDir" . "/" . $ogId . ".codon.align";
	while(my $line=<FF>){
		chomp($line);
		if($line=~/>(.*)/){
			$seqId = $1;
		}else{
			$codonAlignSeq{$seqId}.=$line;
		}
	}
	close FF;
}
print "finish load codon alignment into hash.\n";

# 逐个提取codonAlign中的序列，然后剔除掉-,再和CDS进行比较
my ($seq);
open WW, ">$checkRltTsv";
my @seqId = keys(%codonAlignSeq);
print "total seq num:" . $#seqId . "\n";
foreach $seqId(@seqId){
	$seq = $codonAlignSeq{$seqId};
	$seq=~s/\-//g;
	if($seq ne $cdsSeq{$seqId}){
		print WW ">cds_$seqId\n";
		print WW $cdsSeq{$seqId} . "\n";
		print WW ">codonAlig_$seqId\n";
		print WW $seq . "\n";
	}
}

close WW;
