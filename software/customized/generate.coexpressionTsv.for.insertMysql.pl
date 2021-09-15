#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--coexpressionStudyTsv \\\n" .
		"--curatedStudyTsv \\\n" .
		"--curatedExperimentTsv \\\n" .
		"--outputMysqlTsv\n";
	exit;
}

my ($coexpressionStudyTsv, $curatedStudyTsv, $curatedExperimentTsv, $outputMysqlTsv);

GetOptions(
        'coexpressionStudyTsv=s'=>\$coexpressionStudyTsv,
        'curatedStudyTsv=s'=>\$curatedStudyTsv,
        'curatedExperimentTsv=s'=>\$curatedExperimentTsv,
        'outputMysqlTsv=s'=>\$outputMysqlTsv,
);

# 读取coexpression中study的基本信息
my (%study, $studyHref, $line, $studyId, $taxonId, $exptNum, $exptIdList);
$studyHref=\%study;
open FF, "<$coexpressionStudyTsv";
# 3702    SRP072300       24      SRX1660576,SRX1660577,SRX1660578
while($line=<FF>){
	chomp($line);
	($taxonId, $studyId, $exptNum, $exptIdList) = split(/\t/, $line);
	$studyHref->{$studyId}->{"taxonId"} = $taxonId;
	$studyHref->{$studyId}->{"exptNum"} = $exptNum;
	$studyHref->{$studyId}->{"exptIdList"} = $exptIdList;
}
close FF;


my ($studyTitle, $studyAbstract);
# 读取curatedExperimentTsv中的Title和abstract
open FF, "<$curatedStudyTsv";
# study_accession study_title     study_abstract 
# DRP001756       illumina mRNA-seq analysis for Arabidopsis thaliana max1-1 mutant       Strigolactones are plant hormone that regu  
<FF>;
while($line=<FF>){
	chomp($line);
	($studyId, $studyTitle, $studyAbstract) = split(/\t/, $line);
	if(exists($studyHref->{$studyId}) and not exists($studyHref->{$studyId}->{"Title"})){

		if($studyTitle ne "" and $studyTitle ne "-" and $studyTitle ne "/" and $studyTitle ne "unknown"){
			$studyHref->{$studyId}->{"Title"} = $studyTitle;
		}else{
			$studyHref->{$studyId}->{"Title"} = "NA";
		}
		
		if($studyAbstract ne "" and $studyAbstract ne "-" and $studyAbstract ne "/" and $studyAbstract ne "unknown"){
			$studyHref->{$studyId}->{"Abstract"} = $studyAbstract;
		}else{
			$studyHref->{$studyId}->{"Abstract"} = "NA";
		}

	}else{
		next;
	}
}
close FF;

# 读取curatedExperimentTsv中的publication信息
my (@fieldName, @fieldValue, %tmpHash, $publication);
open FF, "<$curatedExperimentTsv";
# Ecotype Cultivar  Genotype  Tissue SubTissue  Development  Treatment  treatmentGroup Experiment  Study dataSource  Base Layout SpotsNum  ReadNum SpotLen ReadLen Gather AS  Assemble  RunList Phenotype
$line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@fieldValue = split(/\t/, $line);
	for(my $i=0; $i<=$#fieldValue; $i++){
		$tmpHash{$fieldName[$i]} = $fieldValue[$i];
	}

	$studyId = $tmpHash{"Study"};
	$publication = $tmpHash{"dataSource"};
	if(exists($studyHref->{$studyId}) and not exists($studyHref->{$studyId}->{"publication"})){
		if($publication ne "" and $publication ne "-" and $publication ne "/" and $publication ne "unknown"){
			$studyHref->{$studyId}->{"publication"} = $publication;
		}else{
			$studyHref->{$studyId}->{"publication"} = "NA";
		}
	}
}
close FF;

my ($nameString, $valueString);
# studyId, title, abstract, exptNum, exptIdList, publication

my (@studyId);
@studyId = keys(%study);
open WW, ">$outputMysqlTsv";
foreach $studyId(@studyId){

	$nameString = join("###", 
	"taxonId",
	"studyId", 
	"title",
	"abtract", 
	"exptNum", 
	"exptIdList",
	"publication"
	);

	$studyHref->{$studyId}->{"Title"}=~s/'/|_/g;
	$studyHref->{$studyId}->{"Title"}=~s/"/||_/g;

	$studyHref->{$studyId}->{"Abstract"}=~s/'/|_/g;
	$studyHref->{$studyId}->{"Abstract"}=~s/"/||_/g;

	$valueString = join("###",
		$studyHref->{$studyId}->{"taxonId"},
		$studyId,
		$studyHref->{$studyId}->{"Title"},
		$studyHref->{$studyId}->{"Abstract"},
		$studyHref->{$studyId}->{"exptNum"},
		$studyHref->{$studyId}->{"exptIdList"},
		$studyHref->{$studyId}->{"publication"}
	);
	
	# 输出1条experiment
	print WW $nameString . "___" . $valueString . "\n";
}
close WW;


