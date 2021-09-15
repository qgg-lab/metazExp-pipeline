#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputAsFileList A3SS.catalog,A5SS.catalog,MXE.catalog,RI.catalog,SE.catalog \\\n" .
                "--inputGtfFileList ensembl.gtf,refseq.gt \\\n" .
                "--outputAsFileList new.A3SS.catalog,new.A5SS.catalog,new.MXE.catalog,new.RI.catalog,new.SE.catalog \n\n";
	exit;
}

my ($inputAsFileList, $inputGtfFileList, $outputAsFileList);
GetOptions(
        'inputAsFileList=s'=>\$inputAsFileList,
        'inputGtfFileList=s'=>\$inputGtfFileList,
        'outputAsFileList=s'=>\$outputAsFileList,
);

my (@tmp, @fields, $line, @asFile, $asFile, @gtfFile, $gtfFile, $i, @attr, $attr, $geneId, $symbol, @newAsFile);
my (%geneIdToSymbol);

@asFile = split(/,/, $inputAsFileList);
@newAsFile = split(/,/, $outputAsFileList);
@gtfFile = split(/,/, $inputGtfFileList);

foreach $gtfFile(@gtfFile){
	open FF, "<$gtfFile";
	while($line=<FF>){
		next if(not($line=~/\ttranscript\t/));
		chomp($line);
		@fields = ();
		@fields = split(/\t/, $line);
		@attr = ();
		@attr = split(/;/, $fields[8]);
		$geneId = "";
		$symbol = "NA";
		foreach $attr(@attr){
			#gene_id "gene-LOC107049475"; gene_name "LOC107049475";
			if($attr=~/gene_id "(.*)"/){
				$geneId = "\"" . $1 . "\"";
			}
			if($attr=~/gene_name "(.*)"/){
				$symbol = $1;
			}			
		}
		if(not(exists($geneIdToSymbol{$geneId}))){
			if($symbol eq "NA"){
				$geneIdToSymbol{$geneId} = $symbol;
			}else{
				$geneIdToSymbol{$geneId} = "\"" . $symbol . "\"";
			}
		}elsif($geneIdToSymbol{$geneId} eq "NA"){
			if($symbol eq "NA"){
				#$geneIdToSymbol{$geneId} = $symbol;
			}else{
				$geneIdToSymbol{$geneId} = "\"" . $symbol . "\"";
			}
		}
	}
	close FF;
}

for($i=0; $i<=$#asFile; $i++){
	open FF, "<" . $asFile[$i];
	open WW, ">" . $newAsFile[$i];
	while($line=<FF>){
		#ASID    GeneID  geneSymbol      chr     strand  ...
		#GGALA5SS0000014947      "ENSGALG00000000003"    "PANX2" chr1
		if($line=~/^ASID/){
			print WW $line;
		}else{
			@fields = ();
			@fields = split(/\t/, $line);
#			if(exists($geneIdToSymbol{$fields[1]})){
			print WW $fields[0] . "\t" . $fields[1] . "\t" . $geneIdToSymbol{$fields[1]} ;
#			}else{
#				print WW $fields[0] . "\t" . $fields[1] . "\tNA";
#			}

			for(my $j=3; $j<=$#fields; $j++){
				print WW "\t" . $fields[$j];
			}
			print WW;
		}
	}
	close FF;
	close WW;
}
