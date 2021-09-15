#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonIdListFile \\\n" .
                "--workspaceDir \n";
	exit;
}

my ($taxonIdListFile, $workspaceDir);

GetOptions(
        'taxonIdListFile=s'=>\$taxonIdListFile,
        'workspaceDir=s'=>\$workspaceDir,
);

my (@taxonId, $taxonId, %pepSeqId, %cDNAseqId);
my ($seqIdText, $cmd, $pepSeqFile, $cDNAseqFile);
my (@line, $line, @tmp, $flag, @id, $id);

open FF, "<$taxonIdListFile";
@taxonId = <FF>;
close FF;

foreach $taxonId(@taxonId){
	chomp($taxonId);
	$pepSeqFile = $workspaceDir . "/$taxonId/001-prepare-local-datasource/ensembl.pep.fa";
	$cmd = "grep \">\" $pepSeqFile";
	%pepSeqId = ();
	$seqIdText = `$cmd`;
	@line = ();
	@line = split(/\n/, $seqIdText);
	foreach $line(@line){
		@tmp = ();
		@tmp = split(/ /, $line);
		$id = substr($tmp[0], 1);
		$pepSeqId{$id} = 1;
		#print "pep id: " . $id . "\n";
		#<STDIN>;
	}
	
	$cDNAseqFile = $workspaceDir . "/$taxonId/001-prepare-local-datasource/ensembl.cDNA.fa";
	$cmd = "grep \">\" $cDNAseqFile";
	%cDNAseqId = ();
	$seqIdText = `$cmd`;
	@line = ();
	@line = split(/\n/, $seqIdText);
	foreach $line(@line){
		@tmp = ();
		@tmp = split(/ /, $line);
		$id = substr($tmp[0], 1);
		$cDNAseqId{$id} = 1;
		#print "cDNA id: " . $id . "\n";
		#<STDIN>;
	}
	
	# 检测pep的ID和cDNA的ID是否一致
	$flag = 1;
	@id = ();
	@id = keys(%pepSeqId);
	foreach $id(@id){
		if(not exists($cDNAseqId{$id})){
			print $taxonId . " $id\n";
			last;
		}
	}
}
