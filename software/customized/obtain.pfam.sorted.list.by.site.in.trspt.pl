#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--pfamInCdnaTsv pfam.in.pep.tsv \\\n" .
                "--pfamSrtListBySiteInTrsptTsv pfam.sorted.list.by.site.in.trspt.tsv \n";
	exit;
}

my ($pfamInCdnaTsv, $pfamSrtListBySiteInTrsptTsv);

GetOptions(
        'pfamInCdnaTsv=s'=>\$pfamInCdnaTsv,
        'pfamSrtListBySiteInTrsptTsv=s'=>\$pfamSrtListBySiteInTrsptTsv,
);

my (%trsptToPfamList, $trsptToPfamListHref);
$trsptToPfamListHref = \%trsptToPfamList;
my ($line, $field, @field, @trsptId, $trsptId, $pfamList, @pfam, @tt, $pfamNum);
open FF, "<$pfamInCdnaTsv";
# trsptId pfamId  pfamName        beginInCdna     endInCdna
# AT5G45110.2     PF00651 BTB/POZ domain  680     1054
# AT5G45110.2     PF12313 NPR1/NIM1 like defence protein C terminal       1610    1999
# AT5G45110.2     PF11900 Domain of unknown function (DUF3420)    1172    1315
<FF>;
$pfamNum = 0;
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$trsptToPfamListHref->{$field[0]} .= join("\t", $field[1], $field[3], $field[4]) . "\n";
}
close FF;

open WW, ">$pfamSrtListBySiteInTrsptTsv";
my ($trsptAndSrtPfam);
@trsptId = keys(%trsptToPfamList);
foreach $trsptId(@trsptId){
	@field = ();
	@field = split(/\n/, $trsptToPfamListHref->{$trsptId});
	$pfamNum = 0;
	@pfam = ();
	foreach $field(@field){
		($pfam[$pfamNum][0], $pfam[$pfamNum][1], $pfam[$pfamNum][2]) = split(/\t/, $field);
		$pfamNum = $pfamNum + 1;
	}
	# 按照第1列进行排序
	@pfam = sort {$a->[1] <=> $b->[1]} @pfam;
	
	# 排序后输出pfam
	$trsptAndSrtPfam = "";
	$trsptAndSrtPfam = $trsptId . "\t";
	for(my $i=0; $i<=$#pfam; $i++){
		$trsptAndSrtPfam .= $pfam[$i][0] . "[" . $pfam[$i][1] . "," . $pfam[$i][2] . "]|";
	}
	print WW substr($trsptAndSrtPfam, 0, length($trsptAndSrtPfam) - 1) . "\n";
#	<STDIN>;
}
close WW;
