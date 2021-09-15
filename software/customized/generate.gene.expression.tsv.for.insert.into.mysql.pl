#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputAssemblyDir \\\n" . 
		"--filteredExptInfoTsv \\\n" .
		"--geneIdListTsv \\\n" .
		"--outputGeneExpFile geneExp.tsv \n\n";
	exit;
}

my ($inputAssemblyDir, $filteredExptInfoTsv, $outputGeneExpFile, $geneIdListTsv);

GetOptions(
        'inputAssemblyDir=s'=>\$inputAssemblyDir,
        'filteredExptInfoTsv=s'=>\$filteredExptInfoTsv,
	'geneIdListTsv=s'=>\$geneIdListTsv,
        'outputGeneExpFile=s'=>\$outputGeneExpFile,
);

my (@fieldName, @field, $i, $exptId, $exptIdPos, $geneExpFile);
my ($experimentLine, %experiment, $geneLine, $geneId, $cov, $fpkm, $tpm);

# 把基因编号读入hash
my %geneId;
open FF, "<$geneIdListTsv";
while(my $line=<FF>){
	chomp($line);
	$geneId{$line}=1;	
}
close FF;


open FF, "<$filteredExptInfoTsv";
open WW, ">$outputGeneExpFile";
$experimentLine=<FF>;
chomp($experimentLine);

@fieldName = split(/\t/, $experimentLine);
my $exptIdPos = 0;
for($exptIdPos=0; $exptIdPos<=$#fieldName; $exptIdPos++){
        last if($fieldName[$exptIdPos] eq "Experiment");
}

my (%tmpExptHash, $tmpExptHashHref, $tissue, $treatment);
$tmpExptHashHref=\%tmpExptHash;
# 依次读取每个实验的编号和对应的表达文件
while($experimentLine=<FF>){
        chomp($experimentLine);
        @field = ();
        @field = split(/\t/, $experimentLine);
	for(my $j=0; $j<=$#field; $j++){
		$tmpExptHashHref->{$fieldName[$j]} = $field[$j];
	}
	$tissue = $tmpExptHashHref->{"Tissue"};
	$tissue=~tr/ /_/;
	$treatment= $tmpExptHashHref->{"Treatment"};
	$treatment=~tr/ /_/;

        $exptId = $field[$exptIdPos];
	$geneExpFile = $inputAssemblyDir . "/" . $exptId . "/geneAbundanceByStringtie.tab";
        if(not(-e $geneExpFile)){ 
                next;
        }

	# 读取表达水平文件，只保留geneId, exptId, coverage, fpkm, tpm
	open EXP, "<$geneExpFile";
	<EXP>;
	# GeneID GeneName Reference Strand Start End Coverage FPKM TPM
	@field = ();
	while($geneLine=<EXP>){
		chomp($geneLine);
		@field = split(/\t/, $geneLine);
		$geneId = $field[0];
		$cov = $field[6];
		$fpkm = $field[7];
		$tpm = $field[8];
		# 为了避免麻烦，直接将treatment设置为NA
		$treatment = "NA";
		next if(not exists($geneId{$geneId}));
		print WW "GeneId, ExptId, Tissue, Treatment, Cov, FPKM, TPM___" . $field[0] . ", " . $exptId . ", " . $tissue . ", " . $treatment . ", " . sprintf("%.3f", $cov) . ", " . sprintf("%.3f", $fpkm) . ", " . sprintf("%.3f", $tpm) . "\n";
	}
	close EXP;	
}
close WW;
