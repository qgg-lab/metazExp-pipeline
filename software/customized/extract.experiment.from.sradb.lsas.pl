#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonmyListFile taxonmy.txt\\\n" .
		"--sqlite3 sqlite3 \\\n" .
                "--sradb SRAmetadb.sqlite\\\n" .
		"--minMillionSpot 20 \\\n" .
		"--workDir ./ \\\n" .
                "--outputExpTsvFile experiment.tsv\n";
	exit;
}

my ($taxonmyListFile, $sradb, $workDir, $outputExpTsvFile, $sqlite3, $minMillionSpot);

GetOptions(
        'taxonmyListFile=s'=>\$taxonmyListFile,
	'sqlite3=s'=>\$sqlite3,
        'sradb=s'=>\$sradb,
	'workDir=s'=>\$workDir,
	'minMillionSpot=s'=>\$minMillionSpot,
        'outputExpTsvFile=s'=>\$outputExpTsvFile,
);

# 抽出信息
my ($cmd, @taxonmy, $taxonmy, %exptAttrCatalog, %expt, @exptId, $exptId, $line);
# exptAttrCatalog中登记的是样本的属性名称。
# 主要用于扫描sample_attribute中登记的属性种类，以便于汇总。最后输出每个expt时，即使没有某个属性值，那么也将保留空位，以便插入mysql

open FF, "<$taxonmyListFile";
@taxonmy = <FF>;
foreach $taxonmy(@taxonmy){
	chomp($taxonmy);
}

$cmd = "rm -rf " . $workDir . "/experiment.sra.txt";
system($cmd);

foreach $taxonmy(@taxonmy){
	$cmd = $sqlite3 . " " . $sradb . " \" select taxon_id, \'____\', study_accession, \'____\', study_title, \'____\', study_abstract, \'____\', sample_accession, \'____\', sample_attribute, \'____\', experiment_accession, \'____\', experiment_title, \'____\', spots, \'____\', run_accession, \'____\', sample_alias, \'____\', library_name, \'____\', bases from sra where platform=\'ILLUMINA\' and library_strategy=\'RNA-Seq\' and taxon_id=" . $taxonmy . "\" >> " . $workDir . "/experiment.sra.txt";
	system($cmd);
}

my (%attr, @attr, $attr, $attrName, $attrValue);
my ($taxon_id, $study_accession, $study_title, $study_abstract, $sample_accession, $sample_attribute, $experiment_accession, $experiment_title, $spots, $run_accession, $library_name,  $sample_alias, $bases);
open FF, "<" . $workDir . "/experiment.sra.txt";
while($line=<FF>){
	chomp($line);
	($taxon_id, $study_accession, $study_title, $study_abstract, $sample_accession, $sample_attribute, $experiment_accession, $experiment_title, $spots, $run_accession, $sample_alias, $library_name, $bases) = ("", "", "", "", "", "", "", "", "", "", "", "", "");
	($taxon_id, $study_accession, $study_title, $study_abstract, $sample_accession, $sample_attribute, $experiment_accession, $experiment_title, $spots, $run_accession, $sample_alias, $library_name, $bases) = split(/\|____\|/, $line);
	
	$study_title=~s/"/\\"/g;
	$study_title=~s/'/\\'/g;
	$study_abstract=~s/"/\\"/g;
	$study_abstract=~s/'/\\'/g;
	$experiment_title=~s/"/\\"/g;
	$experiment_title=~s/'/\\'/g;
	$sample_attribute=~s/"/\\"/g;
	$sample_attribute=~s/'/\\'/g;
	$library_name=~s/"/\\"/g;
	$library_name=~s/'/\\'/g;

	${$expt{$experiment_accession}}{"taxon_id"} = $taxon_id;
	${$expt{$experiment_accession}}{"study_accession"} = $study_accession;
	${$expt{$experiment_accession}}{"study_title"} = $study_title;
	${$expt{$experiment_accession}}{"study_abstract"} = $study_abstract;
	${$expt{$experiment_accession}}{"experiment_title"} = $experiment_title;
	${$expt{$experiment_accession}}{"sample_accession"} = $sample_accession;
	${$expt{$experiment_accession}}{"sample_attribute"} = $sample_attribute;
	${$expt{$experiment_accession}}{"experiment_title"} = $experiment_title;
	${$expt{$experiment_accession}}{"library_name"} = $library_name;
	${$expt{$experiment_accession}}{"sample_alias"} = $sample_alias;
	${$expt{$experiment_accession}}{"bases"} += $bases/1000000000;
	${$expt{$experiment_accession}}{"spots"} += $spots/1000000;
	${$expt{$experiment_accession}}{"runList"} .= $run_accession . ",";
	${$expt{$experiment_accession}}{"runNum"}++;

	${$expt{$experiment_accession}}{"breed"} = "";
	${$expt{$experiment_accession}}{"gender"} = "";
	${$expt{$experiment_accession}}{"tissueType"} = "";
	${$expt{$experiment_accession}}{"tissueSubType"} = "";
	${$expt{$experiment_accession}}{"developmentStage"} = "";
	${$expt{$experiment_accession}}{"treatmentType"} = "";
	${$expt{$experiment_accession}}{"treatmentRegulation"} = "";
	${$expt{$experiment_accession}}{"phenotype"} = "";

	# 取样本的 "breed", "gender", "tissue", "development/stage", "treatment"
	@attr = ();
	@attr = split(/ \|\| /, $sample_attribute);
	foreach $attr(@attr){
		($attrName, $attrValue) = ("", "");
		($attrName, $attrValue) = split(/: /, $attr);

		if(index(uc($attrName), "BREED")>=0 or index(uc($attrName), "STRAIN")>=0){
			${$expt{$experiment_accession}}{"breed"} = $attrValue;
		}elsif(index(uc($attrName), "GENDER")>=0){
			${$expt{$experiment_accession}}{"gender"} = $attrValue;
		}elsif(index(uc($attrName), "TISSUE")>=0){
			${$expt{$experiment_accession}}{"tissueType"} = $attrValue;
		}elsif(index(uc($attrName), "DEVELOPMENT")>=0 or index(uc($attrName), "DEV_STAG")>=0 or index(uc($attrName), "AGE")>=0){
			${$expt{$experiment_accession}}{"developmentStage"} = $attrValue;			
		}elsif(index(uc($attrName), "TREATMENT")>=0){
			${$expt{$experiment_accession}}{"treatmentType"} = $attrValue;
		}elsif(index(uc($attrName), "PHENOTYPE")>=0){
			${$expt{$experiment_accession}}{"phenotype"} = $attrValue;
		}
	}

}
close FF;


# 整合exptAttrCatalog中登记的属性和常规属性，将expt的信息按照sql的tsv格式输出
my ($fieldNameString, $valueString);
$fieldNameString = "experiment_accession, taxon_id, study_accession, study_title, study_abstract, sample_attribute, sample_accession, experiment_title, spots, bases, runList, runNum, library_name, sample_alias, breed, gender, tissueType, tissueSubType, development, treatmentType, treatmentRegulation, phenotype";
@attr = ();
@attr = keys(%exptAttrCatalog);
foreach $attr(@attr){
	$fieldNameString .= ", " . $attr;
}

open WW, ">$outputExpTsvFile";
@exptId = keys(%expt);
foreach $exptId(@exptId){
	
	next if(${$expt{$exptId}}{"spots"} < $minMillionSpot);
	print WW $fieldNameString . "_____";

	$valueString = "\"" . $exptId . "\", \"" . ${$expt{$exptId}}{"taxon_id"} . "\", \"" . ${$expt{$exptId}}{"study_accession"} . "\", \"" . ${$expt{$exptId}}{"study_title"} . "\", \"" . ${$expt{$exptId}}{"study_abstract"} . "\", \"" . ${$expt{$exptId}}{"sample_attribute"} . "\", \"" . ${$expt{$exptId}}{"sample_accession"} . "\", \"" . ${$expt{$exptId}}{"experiment_title"} . "\", " . ${$expt{$exptId}}{"spots"} . ", " . ${$expt{$exptId}}{"bases"} . ", \"" . substr(${$expt{$exptId}}{"runList"}, 0, length(${$expt{$exptId}}{"runList"})-1) . "\", " . ${$expt{$exptId}}{"runNum"} . ", \"" . ${$expt{$exptId}}{"library_name"} . "\", \"" . ${$expt{$exptId}}{"sample_alias"} . "\", \"" . ${$expt{$exptId}}{"breed"} . "\", \"" . ${$expt{$exptId}}{"gender"} . "\", \"" . ${$expt{$exptId}}{"tissueType"} . "\", \"" . ${$expt{$exptId}}{"tissueSubType"} . "\", \"" . ${$expt{$exptId}}{"developmentStage"} . "\", \"" .${$expt{$exptId}}{"treatmentType"} . "\", \"" . ${$expt{$exptId}}{"treatmentRegulation"} . "\", \"" . ${$expt{$exptId}}{"phenotype"} . "\"";


	print WW $valueString . "\n";
}
close WW;
