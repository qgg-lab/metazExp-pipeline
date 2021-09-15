#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputSampleInfoFile sample.tsv \\\n" .
                "--inputSampleTissueFile 9031.manual.annot.with.tissue.tsv \\\n" .
		"--inputSamplePrjIdFile prj.infor.tsv \\\n" .
		"--taxonId 9031 \\\n" .
                "--outputSampleFullInfoFile sample.with.full.infomation.tsv \n";
	exit;
}

my ($inputSampleInfoFile, $inputSampleTissueFile, $outputSampleFullInfoFile, $inputSamplePrjIdFile, $taxonId);
GetOptions(
        'inputSampleInfoFile=s'=>\$inputSampleInfoFile,
        'inputSampleTissueFile=s'=>\$inputSampleTissueFile,
	'taxonId=s'=>\$taxonId,
	'inputSamplePrjIdFile=s'=>\$inputSamplePrjIdFile,
        'outputSampleFullInfoFile=s'=>\$outputSampleFullInfoFile,
);

my (%tissue, %sample, $line, @fields, %prj);
open FF, "<$inputSamplePrjIdFile";
while($line=<FF>){
	chomp($line);
	@fields = split(/\t/, $line);
	$prj{$fields[0]}=$fields[1];
}
close FF;

open FF, "<$inputSampleTissueFile"; 
#DRX001562       embryos
#DRX001570       embryos
#DRX001568       embryos
while($line=<FF>){
	chomp($line);
	@fields = split(/\t/, $line);	
	$tissue{$fields[0]}=$fields[1];
}
close FF;

open FF, "<$inputSampleInfoFile";
open WW, ">$outputSampleFullInfoFile";
#expId   status  runNum  runId   library layout  phredScore      readLength      spotNum(M)      alignPer(%)     novelSpliceNum  totalA5SS       totalA3SS       totalSE totalRI totalMXE        novelA5SS       novelA3SS       novelSE novelRI novelMXE        jcecA5SS        jcecA3SS        jcecSE  jcecRI  jcecMXE jcA5SS  jcA3SS  jcSE    jcRI    jcMXE
#SRX1036607      OK      1       SRR2037196      UN      PAIRED  33      101     69.97   94.41%  180187  8214    9989    30316   20246   2916    0       0       11742   2491    1179    7730    9369    29083   19515   2839    7656    9327    28886   19325   2813
while($line=<FF>){
	chomp($line);
	@fields =();
	@fields = split(/\t/, $line);
	if($fields[0] eq "expId"){
		print WW "Taxon\tPrjId\t" . $line . "\tTissue\n";
	}else{
		print WW $taxonId . "\t" . $prj{$fields[0]} . "\t" . $line . "\t" . $tissue{$fields[0]} . "\n";
	}	
}
close FF;
close WW;
