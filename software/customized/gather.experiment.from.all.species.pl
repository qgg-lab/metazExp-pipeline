#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--baseDir \\\n" .
                "--taxonIdList \\\n" .
                "--outputTotalExptTsv \n";
	exit;
}

my ($baseDir, $taxonIdList, $outputTotalExptTsv, $discardExptId);

GetOptions(
        'baseDir=s'=>\$baseDir,
        'taxonIdList=s'=>\$taxonIdList,
	'discardExptId=s'=>\$discardExptId,
        'outputTotalExptTsv=s'=>\$outputTotalExptTsv,
);


my (%expt, $exptHref, @taxonId, $taxonId);
my (%discardExpt, @discardId);

open WW, ">$outputTotalExptTsv";
print WW join("\t", "Taxon", "Study", "Experiment", "Datasource", "Ecotype", "Cultivar", "Genotype", "Tissue", "Development", "TreatmentGroup", "TreatmentDetail", "Volume", "Layout", "phredScore", "SpotNum", "SpotLen", "ReadNum", "ReadLen", "Mapping%", "RunList", "Phenotype", "Improvement", "Status") . "\n";


open FF, "<$taxonIdList";
@taxonId = <FF>;
close FF;


my (@fieldName, @fieldValue, $i, %tmpHash, $exptId, $line, @field);
open FF, "<$discardExptId";
while($exptId=<FF>){
	chomp($exptId);
	$discardExpt{$exptId}=1;
}
close FF;

foreach $taxonId(@taxonId){

	%expt = ();
	$exptHref=\%expt;

	chomp($taxonId);

	# 打开最初总表
	my $curatedExptTsv = $baseDir . "/" . $taxonId . "/001-prepare-local-datasource/curated.experiment.tsv";
	open FF, "<$curatedExptTsv";
	# Ecotype Cultivar        Genotype        Tissue  SubTissue       Development     Treatment       treatmentGroup  Experiment      Study   dataSource      Base    Layout  SpotsNum        ReadNum SpotLen ReadLen Gather  PASS    Assemble        RunList Phenotype
	$line=<FF>;
	chomp($line);
	# 获得字段名
	@fieldName = ();
	@fieldName = split(/\t/, $line);
	# 提取信息读入hash
	while($line=<FF>){
		chomp($line);
		@fieldValue = ();
		@fieldValue = split(/\t/, $line);
		%tmpHash = ();
		for($i=0; $i<=$#fieldValue; $i++){
			$tmpHash{$fieldName[$i]} = $fieldValue[$i];
		}
		$exptId = $tmpHash{"Experiment"};
		for($i=0; $i<=$#fieldValue; $i++){
			$exptHref->{$exptId}->{$fieldName[$i]} = $fieldValue[$i];
		}
		$exptHref->{$exptId}->{"Taxon"} = $taxonId;
		$exptHref->{$exptId}->{"ASSEMBLY"} = 0;
		$exptHref->{$exptId}->{"IDENTIFY"} = 0;
	}
	close FF;


	# 打开组装时的比对表，获得比对信息：alignPercent, libraryType, phredScore
	# Ecotype Cultivar        Genotype        Tissue  SubTissue       Development     Treatment       Experiment      Study   Base    Layout  SpotsNum        ReadNum SpotLen ReadLen Gather  PASS    Assemble        RunList Phenotype       alignPercent    mappedBases     mappedReadNum   detectedReadLen libraryType     phredScore
	my $assemblyAlignTsv = $baseDir . "/" . $taxonId . "/004-combine-assemblies-and-annos/alignment.info.of.assembled.experiment.tsv";
	open FF, "<$assemblyAlignTsv";	
	$line=<FF>;
	chomp($line);
	# 获得字段名
	@fieldName = ();
	@fieldName = split(/\t/, $line);
	# 提取信息读入hash
	while($line=<FF>){
		chomp($line);
		@fieldValue = ();
		%tmpHash = ();
		@fieldValue = split(/\t/, $line);
		for($i=0; $i<=$#fieldValue; $i++){
			$tmpHash{$fieldName[$i]} = $fieldValue[$i];
		}
		$exptId = $tmpHash{"Experiment"};
		$exptHref->{$exptId}->{"alignPercent"} = $tmpHash{"alignPercent"};
		$exptHref->{$exptId}->{"libraryType"} = $tmpHash{"libraryType"};
		$exptHref->{$exptId}->{"phredScore"} = $tmpHash{"phredScore"};
	}
	close FF;


	# 打开组装时过滤表，标注哪些experiment用于组装
	my $filterTsv = $baseDir . "/" . $taxonId . "/004-combine-assemblies-and-annos/cutoff.info.of.assembled.experiment.tsv";
	# /mnt/ufs18/rs-015/qgg/liujind1/workAS1/112509/004-combine-assemblies-and-annos/../002-assemble-trsptome-on-goodExps/psiOutputDir/SRX1140356/transcriptomeByStringtie.gtf	11.56   3.85    1.93 
	open FF, "<$filterTsv";
	while($line=<FF>){
		chomp($line);
		my @field = ();
		@field = split(/\//, $line);
		$exptId = $field[$#field-1];
		$exptHref->{$exptId}->{"ASSEMBLY"} = 1;
	}
	close FF;



	# 打开第10步的比对文件，获得比对信息：alignPercent libraryType phredScore
	my $alignTsv = $baseDir . "/" . $taxonId . "/010-gather-alignment-info-of-all-expts/alignment.info.of.assembled.experiment.tsv";
	open FF, "<$alignTsv";	
	# Cultivar        Genotype        Tissue  SubTissue       Development     Treatment       Experiment      Study   Base    Layout  SpotsNum        ReadNum SpotLen ReadLen Gather  AS      Assemble        RunList Phenotype       alignPercent    mappedBases     mappedReadNum   detectedReadLen libraryType     phredScore 
	$line=<FF>;
	chomp($line);
	# 获得字段名
	@fieldName = ();
	@fieldName = split(/\t/, $line);
	# 提取信息读入hash
	while($line=<FF>){
		chomp($line);
		@fieldValue = ();
		%tmpHash = ();
		@fieldValue = split(/\t/, $line);
		for($i=0; $i<=$#fieldValue; $i++){
			$tmpHash{$fieldName[$i]} = $fieldValue[$i];
		}
		$exptId = $tmpHash{"Experiment"};
		$exptHref->{$exptId}->{"alignPercent"} = $tmpHash{"alignPercent"};
		$exptHref->{$exptId}->{"libraryType"} = $tmpHash{"libraryType"};
		$exptHref->{$exptId}->{"phredScore"} = $tmpHash{"phredScore"};
	}
	close FF;

	# 打开过滤文件表，标注那些experiment被用于最后搜集AS
	my $filterTsv = $baseDir . "/" . $taxonId . "/010-gather-alignment-info-of-all-expts/cutoff.info.of.assembled.experiment.tsv";
	# /mnt/ufs18/rs-015/qgg/liujind1/workAS1/112509/004-combine-assemblies-and-annos/../002-assemble-trsptome-on-goodExps/psiOutputDir/SRX1140356/transcriptomeByStringtie.gtf	11.56   3.85    1.93 
	open FF, "<$filterTsv";
	while($line=<FF>){
		my @field = ();
		@field = split(/\//, $line);
		$exptId = $field[$#field-1];
		$exptHref->{$exptId}->{"IDENTIFY"} = 1;
	}
	close FF;

	# 输出总表
	my @exptId = ();
	@exptId = keys(%expt);
	foreach $exptId(@exptId){


		if($exptHref->{$exptId}->{"ASSEMBLY"}==1){
			$exptHref->{$exptId}->{"ASSEMBLY"} = "Yes";
		}else{
			$exptHref->{$exptId}->{"ASSEMBLY"} = "No";
		}

		if($exptHref->{$exptId}->{"IDENTIFY"}==0){
			$exptHref->{$exptId}->{"IDENTIFY"} = "discarded";
		}else{
			$exptHref->{$exptId}->{"IDENTIFY"} = "adopted";
		}
		if(exists($discardExpt{$exptId})){
			$exptHref->{$exptId}->{"IDENTIFY"} = "discarded";
		}

		print WW join("\t",

		$taxonId,		
		$exptHref->{$exptId}->{"Study"}, 
		$exptId,
		$exptHref->{$exptId}->{"dataSource"}, 
		$exptHref->{$exptId}->{"Ecotype"}, 
		$exptHref->{$exptId}->{"Cultivar"}, 
		$exptHref->{$exptId}->{"Genotype"}, 
		$exptHref->{$exptId}->{"Tissue"}, 
		$exptHref->{$exptId}->{"Development"}, 
		$exptHref->{$exptId}->{"Treatment"}, 
		$exptHref->{$exptId}->{"treatmentGroup"}, 
		$exptHref->{$exptId}->{"Base"}, 
		$exptHref->{$exptId}->{"Layout"}, 
		$exptHref->{$exptId}->{"phredScore"}, 
		$exptHref->{$exptId}->{"SpotsNum"}, 
		$exptHref->{$exptId}->{"SpotLen"}, 
		$exptHref->{$exptId}->{"ReadNum"}, 
		$exptHref->{$exptId}->{"ReadLen"}, 
		$exptHref->{$exptId}->{"alignPercent"}, 
		$exptHref->{$exptId}->{"RunList"}, 
		$exptHref->{$exptId}->{"Phenotype"},
		$exptHref->{$exptId}->{"ASSEMBLY"}, 
		$exptHref->{$exptId}->{"IDENTIFY"}


) . "\n";

	}
#	<STDIN>;
#	exit;
}
