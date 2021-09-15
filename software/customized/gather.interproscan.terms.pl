#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputInterproAnnoNotInPepFile final.complete.trspt.interproscan.annotation.not.in.pep.file.tsv\\\n" .
		"--inputInterproAnnoInPepFile final.complete.trspt.interproscan.annotation.in.pep.file.tsv \\\n" .
		"--outputPurePfamAndGoTermFile final.complete.trspt.pfam.go.anno.tsv \n\n\n";
	exit;
}
my ($id);
my ($inputInterproAnnoNotInPepFile, $inputInterproAnnoInPepFile, $outputPurePfamAndGoTermFile);

GetOptions(
	'inputInterproAnnoNotInPepFile=s'=>\$inputInterproAnnoNotInPepFile,
	'inputInterproAnnoInPepFile=s'=>\$inputInterproAnnoInPepFile,
	'outputPurePfamAndGoTermFile=s'=>\$outputPurePfamAndGoTermFile,
);

my (@pfam, @go, $pfam, $go, $goTerm, $line, %isoformAnno, @fields);
# parse pfam and GO annotation from cDNAseqNotInnPepFile.fa
open FF, "<$inputInterproAnnoNotInPepFile";
# ENSBTAT00000037678_62   102ed91c771096fcc2af830fdceaca84        2549    Pfam    PF00878 Cation-independent mannose-6-phosphate receptor repeat  1666     1791    1.3E-18 T       28-05-2019      IPR000479       Cation-independent mannose-6-phosphate receptor GO:0007041|GO:0038023|GO:0051219
while(my $line=<FF>){

	@fields = split(/\t/, $line);
	$id = $fields[0];
	if($id=~/^(.*)_\d+/){
		$id = $1;
	}

	if($line=~/^.*\tPfam\t(PF\d+)\t.*\n/){
		$pfam = $1;
		if(not exists(${$isoformAnno{$id}}{"pfam"})){
			${$isoformAnno{$id}}{"pfam"} = $pfam . ",";
		}elsif(exists(${$isoformAnno{$id}}{"pfam"}) and index(${$isoformAnno{$id}}{"pfam"}, $pfam)<0){
			${$isoformAnno{$id}}{"pfam"}.= $pfam . ",";
		}
	}

	if($#fields == 13 and $fields[13]=~/^(GO:\d+.*)\n/){
		$goTerm = $1;
		@go = ();
		@go = split(/\|/, $goTerm);
		foreach $go(@go){
			if(not exists(${$isoformAnno{$id}}{"GO"})){
				${$isoformAnno{$id}}{"GO"} = $go . ",";
			}elsif(exists(${$isoformAnno{$id}}{"GO"}) and index(${$isoformAnno{$id}}{"GO"}, $go)<0){
				${$isoformAnno{$id}}{"GO"}.=$go . ",";
			}
		}
	}
}
close FF;

# parse pfam and GO annotation from cDNAseqInnPepFile.fa
open FF, "<$inputInterproAnnoInPepFile";
# ENSBTAT00000037678_62   102ed91c771096fcc2af830fdceaca84        2549    Pfam    PF00878 Cation-independent mannose-6-phosphate receptor repeat  1666     1791    1.3E-18 T       28-05-2019      IPR000479       Cation-independent mannose-6-phosphate receptor GO:0007041|GO:0038023|GO:0051219
while(my $line=<FF>){

	@fields = split(/\t/, $line);
	$id = $fields[0];
	if($line=~/^.*\tPfam\t(PF\d+)\t.*\n/){
		$pfam = $1;
		if(not exists(${$isoformAnno{$id}}{"pfam"})){
			${$isoformAnno{$id}}{"pfam"} = $pfam . ",";
		}elsif(exists(${$isoformAnno{$id}}{"pfam"}) and index(${$isoformAnno{$id}}{"pfam"}, $pfam)<0){
			${$isoformAnno{$id}}{"pfam"}.= $pfam . ",";
		}
	}

	if($#fields == 13 and $fields[13]=~/^(GO:\d+.*)\n/){
		$goTerm = $1;
		@go = ();
		@go = split(/\|/, $goTerm);
		foreach $go(@go){
			if(not exists(${$isoformAnno{$id}}{"GO"})){
				${$isoformAnno{$id}}{"GO"} = $go . ",";
			}elsif(exists(${$isoformAnno{$id}}{"GO"}) and index(${$isoformAnno{$id}}{"GO"}, $go)<0){
				${$isoformAnno{$id}}{"GO"}.=$go . ",";
			}
		}
	}
}
close FF;

# output pfam and go annotation
my ($outputLine);
open WW, ">$outputPurePfamAndGoTermFile";
my @id = keys(%isoformAnno);
foreach $id(@id){
	$outputLine = $id . "\t";
	
	if(exists(${$isoformAnno{$id}}{"pfam"})){
		${$isoformAnno{$id}}{"pfam"}=substr(${$isoformAnno{$id}}{"pfam"}, 0, length(${$isoformAnno{$id}}{"pfam"})-1);
		$outputLine.= ${$isoformAnno{$id}}{"pfam"};
	}else{
		$outputLine.= "";
	}

	$outputLine.= "\t";

	if(exists(${$isoformAnno{$id}}{"GO"})){
		${$isoformAnno{$id}}{"GO"} =substr(${$isoformAnno{$id}}{"GO"}, 0, length(${$isoformAnno{$id}}{"GO"})-1);
		$outputLine.= ${$isoformAnno{$id}}{"GO"};
	}else{
		$outputLine.= "";
	}
	if($outputLine ne ""){
		print WW $outputLine . "\n";
	}
}
close WW;
