#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--goTermInTrspt go.term.in.Trspt.tsv\\\n" .
                "--goUniqTermInTrspt go.uniq.term.in.Trspt.tsv \n";
	exit;
}

my ($goTermInTrspt, $goUniqTermInTrspt);

GetOptions(
        'goTermInTrspt=s'=>\$goTermInTrspt,
        'goUniqTermInTrspt=s'=>\$goUniqTermInTrspt,
);

my ($line, @go, $go, $i, $trsptId, $goIdList);
my (%trsptToGoList);
# 打开
open FF, "<$goTermInTrspt";
# SRX399569.19052.6       GO:0008270
# SRX399569.19052.6       GO:0008270
# DRX092559.14802.4       GO:0003735|GO:0005840|GO:0006412
while($line=<FF>){
	chomp($line);
	($trsptId, $goIdList) = split(/\t/, $line);
	@go = ();
	@go = split(/\|/, $goIdList);
	foreach $go(@go){
		if(not exists($trsptToGoList{$trsptId})){
			$trsptToGoList{$trsptId} = $go;
		}elsif(index($trsptToGoList{$trsptId}, $go)>=0){
			
		}else{
			$trsptToGoList{$trsptId} .= "|" . $go;
		}
	}
}
close FF;

open WW, ">$goUniqTermInTrspt";
my @trsptId = keys(%trsptToGoList);
foreach my $trsptId(@trsptId){
	# 重新整理GO词汇，对Go词汇排序后再输出
	@go = ();
	@go = split(/\|/, $trsptToGoList{$trsptId});
	@go = sort { $a cmp $b } @go;
	$goIdList = join("|", @go);
	print WW join("\t", $trsptId, $goIdList) . "\n";
}
close WW;
