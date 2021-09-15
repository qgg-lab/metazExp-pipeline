#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--pepFasta \\\n" .
                "--ensemblgtf \\\n" .
		"--longestPepFasta \n";
	exit;
}

my ($pepFasta, $ensemblgtf, $longestPepFasta);

GetOptions(
        'pepFasta=s'=>\$pepFasta,
        'ensemblgtf=s'=>\$ensemblgtf,
        'longestPepFasta=s'=>\$longestPepFasta,
);

my ($line, @field, $field, %pepIdToGeneId, $geneId, $pepId, %geneToPepIdList);

# 读取gtf，获得pep对应geneId
open FF, "<$ensemblgtf";
# scaf1   SilkDB  CDS     17977   18144   .       -       0       gene_id "BGIBMGA001083"; transcript_id "BGIBMGA001083-RA"; exon_number "2"; gene_source "SilkDB"; gene_biotype "protein_coding"; transcript_source "SilkDB"; transcript_biotype "protein_coding"; protein_id "BGIBMGA001083-TA";
while($line=<FF>){
	next if(not($line=~/\tCDS\t/));
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	($geneId, $pepId) = ("", "");
	&getGeneIdAndPepId($field[8], \$geneId, \$pepId);
	if(not exists($geneToPepIdList{$geneId})){
		$geneToPepIdList{$geneId} = $pepId;
	}elsif(index($geneToPepIdList{$geneId}, $pepId)<0){
		$geneToPepIdList{$geneId} .=  "#" . $pepId;
	}
}
close FF;

# 读取蛋白序列进入哈希
my (%pepSeq);
open FF, "<$pepFasta";
# >BGIBMGA010799-TA pep supercontig:ASM15162v1:scaf25:1994697:2000698:1
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@field = ();
		@field = split(/ /, $1);
		$pepId = $field[0];
	}else{
		$pepSeq{$pepId}.=$line;
	}
}
close FF;

# 扫描每一个gene中蛋白序列，输出最长的pep序列
open WW, ">$longestPepFasta";
my (@geneId, @pepId, $longestPepSeq, $longestPepLen);
@geneId = keys(%geneToPepIdList);
foreach $geneId(@geneId){
	@pepId= ();
	@pepId = split(/#/, $geneToPepIdList{$geneId});
	$longestPepLen = 0;
	$longestPepSeq = "";
	foreach $pepId(@pepId){
		if(length($pepSeq{$pepId})>$longestPepLen){
			$longestPepLen = length($pepSeq{$pepId});
			$longestPepSeq = $pepSeq{$pepId};
		}
	}
	if($longestPepSeq ne ""){
		print WW ">$geneId\n";
		print WW $longestPepSeq . "\n"; 
	}
}
close WW;

sub getGeneIdAndPepId{
	my ($attrString, $geneId, $pepId) = @_;
	my (@attr, $attr);
	@attr = split(/;/, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$geneId = $1;
		}
		if($attr=~/protein_id "(.*)"/){
			$$pepId = $1;
		}
	}
}
