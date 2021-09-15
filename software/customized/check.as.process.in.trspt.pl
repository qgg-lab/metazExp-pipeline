#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--asProcessTrsptTsv \\\n" .
                "--origCdnaFasta \\\n" .
                "--outputErrProcessTsv \n";
	exit;
}

my ($asProcessTrsptTsv, $origCdnaFasta, $outputErrProcessTsv);

GetOptions(
        'asProcessTrsptTsv=s'=>\$asProcessTrsptTsv,
        'origCdnaFasta=s'=>\$origCdnaFasta,
        'outputErrProcessTsv=s'=>\$outputErrProcessTsv,
);

# 将原始的cDNA序列文件读入到hash中
my (%origTrspt, $origTrsptId, $line, @tt);
open FF, "<$origCdnaFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$origTrsptId = $1;
		@tt = split(/ /, $origTrsptId);
		$origTrsptId = $tt[0];
	}else{
		$origTrspt{$origTrsptId}.=$line;
	}
}
close FF;


# 读取process文件，对原始cDNA序列进行操作生成新的序列和process中记录的新序列进行比较.
