#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--ensemblGtf \\\n" .
                "--ensemblPepFile \\\n" .
		"--outputPepWithTrsptIdFile \n";
	exit;
}

my ($ensemblGtf, $ensemblPepFile, $outputPepWithTrsptIdFile);

GetOptions(
        'ensemblGtf=s'=>\$ensemblGtf,
        'ensemblPepFile=s'=>\$ensemblPepFile,
        'outputPepWithTrsptIdFile=s'=>\$outputPepWithTrsptIdFile,
);
my (%pepIdToTrsptId, $line, $pepId, $trsptId, @field, $field, $geneId);

open FF, "<$ensemblGtf";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "CDS");

	($geneId, $trsptId, $pepId) = ("", "", "");
	&getId($field[8], \$geneId, \$trsptId, \$pepId);
	#print join("\t", $geneId, $trsptId, $pepId);
	$pepIdToTrsptId{$pepId} = $trsptId if($pepId ne "");
}
close FF;

my (%pep);
open FF, "<$ensemblPepFile";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$pepId = $1;
		@field = ();
		@field = split(/ /, $pepId);
		$pepId = $field[0];
	}else{
		$pep{$pepId} .= $line;
	}
}
close FF;

my @pepId = keys(%pep);
open WW, ">$outputPepWithTrsptIdFile";
foreach $pepId(@pepId){
#	print $pepId;
#	<STDIN>;
#	print $pepIdToTrsptId{$pepId};
#	<STDIN>;
	next if(not exists($pepIdToTrsptId{$pepId}));
	print WW ">" . $pepIdToTrsptId{$pepId} . "\n";
	print WW $pep{$pepId} . "\n";
}
close WW;


sub getId{
	my ($attrString, $geneId, $trsptId, $pepId) =@_;
	my (@attr, $attr);
#	print $attrString;
#	<STDIN>;
	@attr = split(/;/, $attrString);
#	print $#attr;
#	<STDIN>;
	foreach$attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$geneId = $1;
#			print "geneId: " . $$geneId;
#			<STDIN>;
		}elsif($attr=~/transcript_id "(.*)"/){
			$$trsptId = $1;
#			print "trsptId: " . $$trsptId;
#			<STDIN>;
		}elsif($attr=~/protein_id "(.*)"/){
			$$pepId = $1;
#			print "pepId: " . $$pepId;
#			<STDIN>;
		}
	}
}
