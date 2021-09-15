#!/usr/bin/perl
use strict;
my ($gtfFile, $repFile, $genomeFasta, $coveredGtf) = @ARGV;
if($#ARGV<2){
	print "$0 exon.gtf replace.bed genome.fa exon.covered.tsv\n";
	exit;
}

# 读入genome序列
my (%genomeSeq, $genomeSeqHref);
$genomeSeqHref=\%genomeSeq;
my ($line, @chr, $chr);
my $timeLabel = localtime();
print "begin loading genome sequence into hash:" . $timeLabel . "\n";
open FF, "<$genomeFasta";
while($line=<FF>){
        chomp($line);
        if($line=~/>(.*)/){
                @chr = split(/ /, $1);
                $chr = $chr[0];
        }else{
                $genomeSeqHref->{$chr}.=$line;
        }
}
close FF;
$timeLabel = localtime();
print "finish loading genome sequence into hash:" . $timeLabel . "\n";


# 读入gtfFile
$timeLabel = localtime();
print "begin masking gtf region with M: " . $timeLabel . "\n";
my (@field, $len);
open FF, "<$gtfFile";
while($line=<FF>){
#	print $line;
	@field = ();
	@field = split(/\t/, $line);
	$len=$field[4]-$field[3] + 1;
#	print $len . "\n";
	substr($genomeSeqHref->{$field[0]}, $field[3]-1, $len) = "M"x$len;
#	print substr($genomeSeqHref->{$field[0]}, $field[3]-1, $len);
#	<STDIN>;
}
close FF;
$timeLabel = localtime();
print "finish masking gtf region with M: " . $timeLabel . "\n";


# 读入replace，只关注替换后序列变短的情况
$timeLabel = localtime();
print "begin mask replace region with R: " . $timeLabel . "\n";
my ($cutBegin, $cutEnd, $cutLen);
open FF, "<$repFile";
# 1       788799  788805  AGGTCG      A
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	# 关注替换后序列变短的情况
	#   插入
	next if(length($field[4]) > length($field[3]));
	#   整体替换
	next if(length($field[4]) == length($field[3]) and $field[4] ne "*");
#	print $line . "\n";
#	print "original seq:\n";
#	print substr($genomeSeqHref->{$field[0]}, $field[1], $field[2] - $field[1]);
	# 被切掉的区域
	if($field[4] eq "*"){  # 为整体删除类型
		$cutBegin = $field[1];
	}else{
		$cutBegin = $field[1] +  length($field[4]);
	}

	$cutLen = length($field[3]) - length($field[4]);
	substr($genomeSeqHref->{$field[0]}, $cutBegin, $cutLen) = "R"x$cutLen;
#	print "replaced seq:\n";
#	print substr($genomeSeqHref->{$field[0]}, $field[1], $field[2] - $field[1]);
#	<STDIN>;

}

$timeLabel = localtime();
print "finish masking replace region with M: " . $timeLabel . "\n";

# 重新读取gtf的区域，提取genome sequence，判读：
#   如果只有M，那么表示该区域不受replace影响
#   如果只有R，那么表示该区域将被直接切掉
#   如果有M和R，那么表示该区域部分切掉，可能左边，也可能右边
$timeLabel = localtime();
print "begin output replace region: " . $timeLabel . "\n";
my ($seq, $st);
open WW, ">$coveredGtf";
open FF, "<$gtfFile";
while($line=<FF>){
	@field = ();
	@field = split(/\t/, $line);
	$len=$field[4]-$field[3] + 1;
	$seq = substr($genomeSeqHref->{$field[0]}, $field[3]-1, $len);
	$seq = uc($seq);
	$st = $seq;
	$seq=~s/R//g;
	if($seq eq ""){
		print WW "completeCut\t" . $line;
	}else{
		$st=~s/M//g;
		if($st eq ""){
			print WW "NotCut\t" . $line;
		}else{
			print WW "PartCut\t" . $line;
		}
	}
}
close FF;
close WW;
$timeLabel = localtime();
print "finish output replace region:" . $timeLabel . "\n";
