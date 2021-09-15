#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonIdListFile \\\n" .
                "--desDir /mnt/research/qgg/liujind1/workAS1/39947/031-annotate-symbol-by-swissprot/\\\n" .
                "--srcDir /mnt/research/qgg/liujind1/workAS1/fungi.workspace/swissprot\n";
	exit;
}

my ($taxonIdListFile, $desDir, $srcDir);

GetOptions(
        'taxonIdListFile=s'=>\$taxonIdListFile,
        'desDir=s'=>\$desDir,
        'srcDir=s'=>\$srcDir,
);

my (@taxonId, $taxonId, $srcGzFile, $desGzFile);
open FF, "<$taxonIdListFile";
@taxonId = <FF>;
foreach $taxonId(@taxonId){
	chomp($taxonId);
	# sp
	# uniprot-reviewed_yes+AND+organism__871575_.fasta.gz
	# 39947/031-annotate-symbol-by-swissprot/39947.sp.fasta
	$srcGzFile = "$srcDir/uniprot-reviewed_yes+AND+organism__" . $taxonId . "_.fasta.gz";
	$desGzFile = "$desDir/$taxonId/031-annotate-symbol-by-swissprot/" . $taxonId . ".sp.fasta.gz";
	if(-e $srcGzFile){
		system("cp $srcGzFile $desGzFile");
		system("gunzip $desGzFile");
	}else{
		system("touch $desDir/$taxonId/031-annotate-symbol-by-swissprot/" . $taxonId . ".sp.fasta");
	}

	# tr
	# uniprot-reviewed_no+AND+organism__663331_.fasta.gz
	# 39947/031-annotate-symbol-by-swissprot/39947.tr.fasta
	$srcGzFile = "$srcDir/uniprot-reviewed_no+AND+organism__" . $taxonId . "_.fasta.gz";
	$desGzFile = "$desDir/$taxonId/031-annotate-symbol-by-swissprot/" . $taxonId . ".tr.fasta.gz";
	if(-e $srcGzFile){
		system("cp $srcGzFile $desGzFile");
		system("gunzip $desGzFile");
	}else{
		system("touch $desDir/$taxonId/031-annotate-symbol-by-swissprot/" . $taxonId . ".tr.fasta");
	}
}
close FF;
