#!/usr/bin/perl
use strict;
use LWP::Simple;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonIdListFile \\\n" .
                "--inputDir \\\n" .
                "--outputDir \n";
	exit;
}

my ($inputDir, $outputDir, $taxonIdListFile);

GetOptions(
        'taxonIdListFile=s'=>\$taxonIdListFile,
        'inputDir=s'=>\$inputDir,
        'outputDir=s'=>\$outputDir
);

my (@taxonId, $taxonId,$studyId, @line, $line, $title, $titleLine,$doc);
open FF, "<$taxonIdListFile";
@taxonId = <FF>;
close FF;

foreach $taxonId(@taxonId){
	chomp($taxonId);
	open FF, "<$inputDir/$taxonId";	
	open WW, ">$outputDir/$taxonId.title";
	<FF>;
	while($studyId=<FF>){
		chomp($studyId);
		#print $studyId . ": begin obtaion...\n";
		$doc = get("https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=" . $studyId);
		#open TT, ">DOC";
		#print TT $doc;
		#close FF;
		#<STDIN>;
		$titleLine = "";
		$title = "";
		@line = ();
		@line = split(/\n/, $doc);
		for(my $i=0; $i<=$#line; $i++){
			#print $line[$i];
			#<STDIN>;
			if($line[$i]=~/<div class=\"study-container\">/){
				#print $line[$i];
				#<STDIN>;
				$titleLine = $line[$i+1];
				if($titleLine=~/<h1>(.*)<\/h1>/){
					$title = $1;
					last;
				}
			}
		}
		if($title ne ""){
			print WW join("\t", $studyId, $title) . "\n";
		}else{
			print WW join("\t", $studyId, "None") . "\n";
		}
	}
	close FF;
	close WW;
}
