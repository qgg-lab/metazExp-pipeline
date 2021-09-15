#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonIdFileList pure.id.iist\\\n" .
                "--dir /mnt/research/qgg/liujind1/workAS1/fungi.workspace/upload_20210205 \n";
	exit;
}

my ($taxonIdFileList, $dir);

GetOptions(
        'taxonIdFileList=s'=>\$taxonIdFileList,
        'dir=s'=>\$dir,
);

# 237631-uma.Gene2KEGG.xls  237631-uma.KEGG2Gene.xls  237631-uma.KEGG.png  237631-uma.matrix
# 237631-uma.htm            237631-uma.KEGG.pdf       237631-uma_map       237631-uma.path 

my (@taxonId, $taxonId, $cmd);
open FF, "<$taxonIdFileList";
@taxonId = <FF>;
close FF;

foreach $taxonId(@taxonId){
	chomp($taxonId);

	$cmd = "mv $dir/$taxonId/result/$taxonId-*.Gene2KEGG.xls $dir/$taxonId/result/$taxonId.Gene2KEGG.xls";
#	system($cmd);

	$cmd = "mv $dir/$taxonId/result/$taxonId-*.htm $dir/$taxonId/result/$taxonId.htm";
#	system($cmd);

	$cmd = "mv $dir/$taxonId/result/$taxonId-*.KEGG2Gene.xls $dir/$taxonId/result/$taxonId.KEGG2Gene.xls";
#	system($cmd);

	$cmd = "mv $dir/$taxonId/result/$taxonId-*.KEGG.pdf $dir/$taxonId/result/$taxonId.KEGG.pdf";
#	system($cmd);

	$cmd = "mv $dir/$taxonId/result/$taxonId-*.KEGG.png $dir/$taxonId/result/$taxonId.KEGG.png";
#	system($cmd);

	$cmd = "mv $dir/$taxonId/result/$taxonId-*_map $dir/$taxonId/result/$taxonId" . "_map";
	system($cmd);

	$cmd = "mv $dir/$taxonId/result/$taxonId-*.matrix $dir/$taxonId/result/$taxonId.matrix";
#	system($cmd);

	$cmd = "mv $dir/$taxonId/result/$taxonId-*.path $dir/$taxonId/result/$taxonId.path";
#	system($cmd);
}
