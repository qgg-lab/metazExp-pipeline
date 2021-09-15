#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--mafListFile maf.list \\\n" .
                "--aSpeciesName \"oryza_sativa\"\\\n" .
                "--bSpeciesName \"arabidopsis_thaliana\" \\\n" .
		"--outputCoordFile coords.of.osat.in.atha.tsv \n";
	print "aSpeciesName为第一个物种名称，必须要和maf文件中一致！比如中间有下划线\n";
	print "只输出a物种所有碱基在B物种基因组中的位置，即获得A物种每个碱基在B物种中的坐标。\n";
	exit;
}

my ($mafListFile, $aSpeciesName, $bSpeciesName, $outputCoordFile);
GetOptions(
        'mafListFile=s'=>\$mafListFile,
        'aSpeciesName=s'=>\$aSpeciesName,
        'bSpeciesName=s'=>\$bSpeciesName,
        'outputCoordFile=s'=>\$outputCoordFile,
);

my ($line);
my (@mafFile, $mafFile, @lline, $lline, $aSpeciesLine, $bSpeciesLine);
my ($aChr, $aChain, $aStart, $aSequence);
my ($bChr, $bChain, $bStart, $bSequence);

open WW, ">$outputCoordFile";
print WW join("\t", $aSpeciesName . "_chr", $aSpeciesName . "_chain", $aSpeciesName . "_Pos", $bSpeciesName . "_chr", $bSpeciesName . "_chain", $bSpeciesName . "_Pos") . "\n";
open MAF, "<$mafListFile";
@mafFile = <MAF>;
foreach $mafFile(@mafFile){
	chomp($mafFile);
}
# ##maf version=1 program=LASTZ_NET
# #
# a# id: 91480001883549
#  score=30192
#  s oryza_sativa.12         3432083        281 +   27531856 CGTCGGCCTCG--TTCGACGCGGCGGGGTTCGAGGCCGAGCGGCTCCGCCTCGACGCAGAggcgcgggccgggatggcgtccgcggcggcggtggcgggggcggaggcggcg---------gACCCCAAGGCGTGGAAGTGGGCGATACGGAAGCGGGTGTGGGACGCGCTGGAGGCGGAGGGCGTCGCGCGCGACCCGCGCCCCGTCCACCACCGCATCCCCAACTTCGATGGCGCCGCAGCCGCCGCGGACTCGGTACGCGCGTGctgtctcggtctctctcctctctct
#  s arabidopsis_thaliana.1    1623108        286 -   30427671 CGACGGCGTCGCTTTCGATGCTGTAGCTTATGAGGCTGATCGGTTGAGCCTCGACGCGGCGGCGATGGAGGATATGGCG---GAGACTGCAAAGAAGGAACTAGAATCGGATCCAGATAGTGACCCGAAAGCGTGGAAATGGGTTATCCGAAAGAAGATGTGGGATCTCATGGAAGCTCGGAACTACGCTATGAGCCCTAGACCTGTTCATCATCGAATCCCTAATTTCGTCGGTGCATCGGCTGCTGCAAGAAAAGTTTGCGAATTTTACCT---TCCTTCTTCACTTTCT
#
# 注意： 3432083和1623108都是based 0的坐标，也就是在基因组中的坐标应该为当前值再加1
my ($aBase, $bBase, $aPos, $bPos);
foreach $mafFile(@mafFile){
	$/="\n\n";
	open FF, "<$mafFile";
	while($line=<FF>){
		# 读取一条记录后初始化
		($aChr, $aStart, $aChain, $aSequence) = ("", "", "", "");
		($bChr, $bStart, $bChain, $bSequence) = ("", "", "", "");

		# 提取两个基因组block的信息
		@lline = ();
		@lline = split(/\n/, $line);
		foreach $lline(@lline){
			if($lline=~/s $aSpeciesName\.(.*?) +(\d+) +(\d+) +(\+|\-) +(\d+) +(.*)/){
				$aChr = $1;
				$aStart = $2 + 1;
				$aChain = $4;
				$aSequence = $6 + 1;			
			}
			if($lline=~/s $bSpeciesName\.(.*?) +(\d+) +(\d+) +(\+|\-) +(\d+) +(.*)/){
				$bChr = $1;
				$bStart = $2 + 1;
				$bChain = $4;
				$bSequence = $6 + 1;			
			}
		}

		# 如果记录无效则放弃
		next if($aChr eq "");

		$aPos = $aStart;
		$bPos = $bStart;
		for(my $i=0; $i<length($aSequence); $i++){
			$aBase = substr($aSequence, $i, 1);
			$bBase = substr($bSequence, $i, 1);
			if($aBase ne "-" and $bBase ne "-"){
				print WW join("\t", $aChr, $aChain , $aPos, $bChr, $bChain , $bPos) . "\n";
			}elsif($aBase ne "-" and $bBase eq "-"){
				# #表示a物种的碱基没有找到在B物种中的坐标。
				print WW join("\t", $aChr , $aChain, $aPos, "#", "#", "#") . "\n";
			}

			$aPos++ if($aBase ne "-");
			$bPos++ if($aBase ne "-");
		}

	}
	close FF;
}
close WW;
