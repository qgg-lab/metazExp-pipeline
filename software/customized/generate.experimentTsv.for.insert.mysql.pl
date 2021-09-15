#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--psiTsv \\\n" .
                "--asExptInfoList \\\n" .
                "--assemblExptInfoList \\\n" .
		"--outputExpTsv ";
	exit;
}

my ($psiTsv, $asExptInfoList, $assemblyExptInfoList, $outputExpTsv);

GetOptions(
        'psiTsv=s'=>\$psiTsv,
        'asExptInfoList=s'=>\$asExptInfoList,
        'assemblyExptInfoList=s'=>\$assemblyExptInfoList,
        'outputExpTsv=s'=>\$outputExpTsv,
);

my (%expt, $exptHref, @titleField, @valueField, $line, $i, $j, $exptIdPos, $exptId);
$exptHref = \%expt;

# 将参与AS识别的所有experiment读入到hash中
open FF, "<$asExptInfoList";
$line = <FF>;
chomp($line);
@titleField=split(/\t/, $line);
for($i=0; $i<=$#titleField; $i++){
	if($titleField[$i] eq "Experiment"){
		$exptIdPos = $i;
	}
}
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	$exptId = $valueField[$exptIdPos];
	for($i=0; $i<=$#valueField; $i++){
		$exptHref->{$exptId}->{$titleField[$i]} = $valueField[$i];
		$exptHref->{$exptId}->{"usage"} = "mappingAS";
	}
}
close FF;

# 将参与提高转录本注释的所有experiment读入到hash中
open FF, "<$assemblyExptInfoList";
$line = <FF>;
chomp($line);
@titleField=split(/\t/, $line);
for($i=0; $i<=$#titleField; $i++){
	if($titleField[$i] eq "Experiment"){
		$exptIdPos = $i;
	}
}
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	$exptId = $valueField[$exptIdPos];
	for($i=0; $i<=$#valueField; $i++){
		$exptHref->{$exptId}->{$titleField[$i]} = $valueField[$i];
		$exptHref->{$exptId}->{"usage"} = "improvingAnno";
	}
}
close FF;

# 将psi值读入
# 如果psi不等于-1，那么将当前AS分配到对应的expt中
my ($fieldString, $valueString, @titleField, @valueField, $i, $j, %tmpPsi, $tmpPsiHref);
$tmpPsiHref=\%tmpPsi;
open FF, "<$psiTsv";
# asId, experimentId, jcecIncl, jcecSkip, jcecInclFormLen, jcecSkipFormLen, jcecPsi, jcIncl, jcSkip, jcFormLen, jcSkipFormLen, jcPsi___ATHARI0000005836, SRX1708373, 4, 10, 198, 124, 0.20032310177706, 4, 10, 197, 124, 0.201135442011354
while($line=<FF>){
	chomp($line);
	($fieldString, $valueString) = split(/___/, $line);
	@titleField = split(/, /, $fieldString);
	@valueField = split(/, /, $valueString);
	for($i=0; $i<=$#valueField; $i++){
		$tmpPsiHref->{$titleField[$i]} = $valueField[$i];
	}
	if($tmpPsiHref->{"jcecPsi"}!=-1 and $tmpPsiHref->{"jcPsi"}!=-1){
		if($tmpPsiHref->{"asId"}=~/.*?A3SS\d+/){
			$exptHref->{$tmpPsiHref->{"experimentId"}}->{"A3SSnum"}++;
		}elsif($tmpPsiHref->{"asId"}=~/.*?A5SS\d+/){
			$exptHref->{$tmpPsiHref->{"experimentId"}}->{"A5SSnum"}++;
		}elsif($tmpPsiHref->{"asId"}=~/.*?RI\d+/){
			$exptHref->{$tmpPsiHref->{"experimentId"}}->{"RInum"}++;
		}elsif($tmpPsiHref->{"asId"}=~/.*?MXE\d+/){
			$exptHref->{$tmpPsiHref->{"experimentId"}}->{"MXEnum"}++;
		}elsif($tmpPsiHref->{"asId"}=~/.*?SE\d+/){
			$exptHref->{$tmpPsiHref->{"experimentId"}}->{"SEnum"}++;
		}
	}
}
close FF;

open WW, ">$outputExpTsv";
# 输出所有的experiment的统计数据
my @exptId = keys(%expt);
# Ecotype Cultivar        Genotype        Tissue  SubTissue       Development     TissueGroup     Treatment       TreatmentGroup  Experiment
#       Study   DataSource      Base    Layout  SpotsNum        ReadNum SpotLen ReadLen Gather  AS      Assemble        RunList Phenotype       alignPercent    mappedBases     mappedReadNum   detectedReadLen libraryType     phredScore
foreach $exptId(@exptId){
	$fieldString = join(", ", "Ecotype", "Cultivar", "Genotype", "Tissue", "SubTissue", "Development", "TissueGroup", "Treatment", "TreatmentGroup", "Experiment", "Study", "DataSource", "TotalBase", "Layout", "SpotsNum", "ReadNum", "SpotLen", "ReadLen", "RunList", "Phenotype", "alignPercent", "mappedBases", "mappedReadNum", "detectedReadLen", "libraryType", "phredScore", "A3SSnum", "A5SSnum", "MXEnum", "RInum", "SEnum");
	$valueString = join(", ", $exptHref->{$exptId}->{"Ecotype"}, $exptHref->{$exptId}->{"Cultivar"}, $exptHref->{$exptId}->{"Genotype"}, $exptHref->{$exptId}->{"Tissue"}, $exptHref->{$exptId}->{"SubTissue"}, $exptHref->{$exptId}->{"Development"}, $exptHref->{$exptId}->{"TissueGroup"}, $exptHref->{$exptId}->{"Treatment"}, $exptHref->{$exptId}->{"TreatmentGroup"}, $exptHref->{$exptId}->{"Experiment"}, $exptHref->{$exptId}->{"Study"}, $exptHref->{$exptId}->{"DataSource"}, $exptHref->{$exptId}->{"Base"}, $exptHref->{$exptId}->{"Layout"}, $exptHref->{$exptId}->{"SpotsNum"}, $exptHref->{$exptId}->{"ReadNum"}, $exptHref->{$exptId}->{"SpotLen"}, $exptHref->{$exptId}->{"ReadLen"}, $exptHref->{$exptId}->{"RunList"}, $exptHref->{$exptId}->{"Phenotype"}, $exptHref->{$exptId}->{"alignPercent"}, $exptHref->{$exptId}->{"mappedBases"}, $exptHref->{$exptId}->{"mappedReadNum"}, $exptHref->{$exptId}->{"detectedReadLen"}, $exptHref->{$exptId}->{"libraryType"}, $exptHref->{$exptId}->{"phredScore"}, $exptHref->{$exptId}->{"A3SSnum"}, $exptHref->{$exptId}->{"A5SSnum"}, $exptHref->{$exptId}->{"MXEnum"}, $exptHref->{$exptId}->{"RInum"}, $exptHref->{$exptId}->{"SEnum"});
	print WW $fieldString . "___" . $valueString . "\n";
}
close WW;
