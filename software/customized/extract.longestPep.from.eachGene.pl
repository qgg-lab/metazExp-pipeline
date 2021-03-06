#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--pepFasta \\\n" .
                "--gtf \\\n" .
                "--longestPepFasta \n";
	exit;
}

my ($pepFasta, $gtf, $longestPepFasta);

GetOptions(
        'pepFasta=s'=>\$pepFasta,
        'gtf=s'=>\$gtf,
        'longestPepFasta=s'=>\$longestPepFasta,
);

my (%gene, $geneHref, $trsptText, @trsptLine, $trsptLine, $geneId, $trsptId, %pep,  $line, @tt, %regGene, @geneId);
$geneHref = \%gene;

# 读取gtf，将trspt聚类到gene中
# 将geneId依次登记到geneId
my $cmd = "grep -P \"\\ttranscript\\t\" " . $gtf;
$trsptText = `$cmd`;
@trsptLine = split(/\n/, $trsptText);
foreach $trsptLine(@trsptLine){
	($geneId, $trsptId) = ("", "");
	&getGeneIdAndTrsptId($trsptLine, \$geneId, \$trsptId);
	$geneHref->{$geneId}->{$trsptId} = 1;
}


# 将所有pep读入hash
open FF, "<$pepFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$trsptId = $1;
		@tt = split(/ /, $trsptId);
		$trsptId = $tt[0];
	}else{
		$pep{$trsptId}.=$line;
	}
}
close FF;


# 重新扫描每个基因内部的trspt长度，挑选最长的pep代表这个基因输出
open WW, ">$longestPepFasta";
my (@trsptId, $maxPepSeq, $maxTrsptId);
@geneId = keys(%gene);
foreach $geneId(@geneId){
	@trsptId = ();
	@trsptId = keys(%{$geneHref->{$geneId}});
	$maxPepSeq = "";
	$maxTrsptId = "";
	foreach $trsptId(@trsptId){
		if(length($maxPepSeq) < length($pep{$trsptId})){
			$maxPepSeq = $pep{$trsptId};
			$maxTrsptId = $trsptId;
		}
	}

	if($maxPepSeq ne ""){
		print WW ">$maxTrsptId\n";
		print WW $maxPepSeq . "\n";
	}
}
close WW;

# 获得基因编号和trspt编号
sub getGeneIdAndTrsptId{
	my ($line, $geneId, $trsptId) = @_;
	my (@tt);
	@tt = split(/\t/, $line);
	if($tt[8]=~/gene_id "(.*?)"; transcript_id "(.*?)"/){
		$$geneId = $1;
		$$trsptId = $2;
	}
}
