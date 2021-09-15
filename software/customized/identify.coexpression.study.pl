#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--workdir /mnt/research/qgg/liujind1/workAS1\\\n" .
                "--taxonIdList /mnt/research/qgg/liujind1/workAS1/workspace-plant/all.taxonId.list \\\n" .
		"--exptInfoTsv filtered.alignment.info.of.assembled.experiment.tsv \\\n" .
		"--minMappedBase 2 \\\n" .
		"--minExptNum 20 \\\n" .
		"--outputTaxonStudyExpIdListTsv \n";
	exit;
}

my ($workdir, $taxonIdList, $exptInfoTsv, $minMappedBase, $minExptNum, $outputTaxonStudyExpIdListTsv);

GetOptions(
        'workdir=s'=>\$workdir,
        'taxonIdList=s'=>\$taxonIdList,
        'exptInfoTsv=s'=>\$exptInfoTsv,
	'minMappedBase=s'=>\$minMappedBase,
	'minExptNum=s'=>\$minExptNum,
	'outputTaxonStudyExpIdListTsv=s'=>\$outputTaxonStudyExpIdListTsv,
);

my (@taxonId, $taxonId, %study, @study, $study, $studyHref, %tmpHash, @fieldName, @fieldValue, $experimentId);
my ($exptInfoFile, $line);
$studyHref=\%study;
#print $outputTaxonStudy;
#<STDIN>;

open WW, ">$outputTaxonStudyExpIdListTsv";
open FF, "<$taxonIdList";
@taxonId = <FF>;
close FF;
foreach $taxonId(@taxonId){
	chomp($taxonId);
	$exptInfoFile = $workdir . "/" . $taxonId . "/010-gather-alignment-info-of-all-expts/" . $exptInfoTsv;
	%study = ();
	open FF, "<$exptInfoFile";
	#Ecotype Cultivar Genotype Tissue SubTissue Development Treatment Experiment Study Base Layout SpotsNum ReadNum SpotLen ReadLen Gather AS Assemble RunList Phenotype alignPercent mappedBases mappedReadNum detectedReadLen libraryType phredScore
	$line = <FF>;
	chomp($line);
	@fieldName = split(/\t/, $line);
	while($line=<FF>){
		chomp($line);
		@fieldValue = ();
		@fieldValue = split(/\t/, $line);
		for(my $i=0; $i<=$#fieldValue; $i++){
			$tmpHash{$fieldName[$i]} = $fieldValue[$i];
		}

		# 检查mappedBases是否达标
		if($tmpHash{"mappedBases"}>=$minMappedBase and -e $workdir . "/" . $taxonId . "/008-pickup-psi-of-ASs-in-all-expers/psiOutputDir/" . $tmpHash{"Experiment"} . "/geneAbundanceByStringtie.tab"){
			$study = $tmpHash{"Study"};
			$experimentId = $tmpHash{"Experiment"};
			$studyHref->{$study}->{"exptNum"}++;
			$studyHref->{$study}->{"exptIdList"} .= $experimentId . ",";
		}
		
	}
	close FF;

	# 将当前物种下study中实验数量大于20的输出
	@study = ();
	@study = keys(%study);
	foreach $study (@study){
		if($studyHref->{$study}->{"exptNum"}>=$minExptNum){
			print WW join("\t", $taxonId, $study, $studyHref->{$study}->{"exptNum"}, substr($studyHref->{$study}->{"exptIdList"}, 0, length($studyHref->{$study}->{"exptIdList"})-1)) . "\n";
		}
	}	
}
close WW;
