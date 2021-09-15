#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputAsGtf \\\n" .
                "--outputAsGff \n";
	exit;
}

my ($inputAsGtf, $outputAsGff);

GetOptions(
        'inputAsGtf=s'=>\$inputAsGtf,
        'outputAsGff=s'=>\$outputAsGff,
);



# 按照geneId将属于该geneId的as并入
my (%asGtf, $asGtfHref);
$asGtfHref = \%asGtf;
open FF, "<$inputAsGtf";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
#	print $field[8];
#	<STDIN>;
	$geneId = &getGeneId($field[8]);
	$asGtfHref->{$geneId} .= $line . "\n";
}
close FF;


# 提取所有geneId输出
my (@geneId);
@geneId = keys(%asGtf);
foreach $geneId(@geneId){
	open WW, ">$outputDir/$geneId.gtf";
	print WW "track name=transcript priority=1\n";
	print WW $trsptGtfHref->{$geneId};
	print WW "track name=AsEvent priority=1\n";
	print WW $asGtfHref->{$geneId};
	close WW;
}

sub getGeneId{
	my ($attrString) = @_;
	my (@attr,  $attr);
#	print $attrString;
	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*?)"/){
			return $1;
		}
	}
}
