#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--orthoGroupListFile \\\n" .
                "--subGroupNum \\\n" .
		"--prefix \\\n" .
		"--outputDir \n";
	exit;
}

my ($orthoGroupListFile, $subGroupNum, $prefix, $outputDir);

GetOptions(
        'orthoGroupListFile=s'=>\$orthoGroupListFile,
        'subGroupNum=s'=>\$subGroupNum,
	'prefix=s'=>\$prefix,
        'outputDir=s'=>\$outputDir,
);

my (@subgroupText, $line, $i);
open FF, "<$orthoGroupListFile";
$line=<FF>;
# 将第一行分到所有的subgroup中
for($i=0; $i<$subGroupNum; $i++){
	$subgroupText[$i] = $line;
}

# 轮流将orthoGroup行分散到各个subgroup中
$i=0;
while($line=<FF>){
	$subgroupText[$i-int($i/$subGroupNum)*$subGroupNum] .= $line;
	$i++;
}
close FF;

# 将所有subgroup写到文件中
for($i=0; $i<$subGroupNum; $i++){
	open WW, ">$outputDir" . "/" . $prefix . $i;
	print WW $subgroupText[$i];
	close WW;
}

