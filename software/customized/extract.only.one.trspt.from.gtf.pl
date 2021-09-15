#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtf orign.gtf\\\n" .
                "--outputSingleTrsptGtf \n";
	exit;
}

my ($gtf, $outputSingleTrsptGtf);

GetOptions(
        'gtf=s'=>\$gtf,
        'outputSingleTrsptGtf=s'=>\$outputSingleTrsptGtf,
);

my (%trspt, $trsptHref, $line, $geneId, $trsptId, @geneId, @trsptId, @field);
$trsptHref=\%trspt;

# 将trspt注释分离开来后分别进入hash
open FF, "<$gtf";
while($line=<FF>){
	@field = ();
	@field = split(/\t/, $line);
	($geneId, $trsptId) = ("", "");
	&getTrsptIdAndGeneId($line, \$geneId, \$trsptId);
	next if($trsptId eq "");
	$trsptHref->{$geneId}->{$trsptId}.=$line;
	
}
close FF;

# 提取每个gene，然后只输出一个trspt的注释
open WW, ">$outputSingleTrsptGtf";
@geneId = keys(%trspt);
foreach $geneId(@geneId){
	@trsptId = ();
	@trsptId = keys($trsptHref->{$geneId});
	if($#trsptId>=0){
		$trsptId = $trsptId[0];
		print WW $trsptHref->{$geneId}->{$trsptId};
	}
}
close WW;


sub getTrsptIdAndGeneId{
	my ($line, $geneId, $trsptId)=@_;
	my (@tt);
	@tt = split(/\t/, $line);
	if($tt[8]=~/gene_id "(.*?)"; transcript_id "(.*?)";/){
		$$geneId = $1;
		$$trsptId = $2;
	}
}
