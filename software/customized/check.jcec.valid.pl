#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--experimentIdList \\\n" .
                "--baseDir \n";
	exit;
}

my ($experimentIdList, $baseDir);

GetOptions(
        'experimentIdList=s'=>\$experimentIdList,
        'baseDir=s'=>\$baseDir,
);

my (@experimentId, $experimentId, $jcecRiTsv,  $line, @field);
open FF, "<$experimentIdList";
@experimentId=<FF>;
close FF;
foreach $experimentId (@experimentId){
	chomp($experimentId);
	$jcecRiTsv = "$baseDir/$experimentId/JCEC.raw.input.RI.txt";
	# ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
	# 0_0     2       133                     231     99
	# 1_0     0       0                       214     99
	# 2_0     1       31                      177     99
	my $checkFlag = 0;
	open FF, "<$jcecRiTsv"; 
	<FF>;
	while($line=<FF>){
		chomp($line);
		@field = split(/\t/, $line);
		if($field[1] !=0 or $field[2]!=0){
			$checkFlag = 1;
			last;
		}
	}
	print $experimentId . "\n" if($checkFlag == 0);
}
