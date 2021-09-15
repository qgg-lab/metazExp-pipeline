#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--genomeFasta genome.fa \\\n" .
                "--asCatalogFileList A3SS.txt,A5SS.txt,MXE.txt,RI.txt,SE.txt \\\n" .
                "--asTypeList A3SS,A5SS,MXE,RI,SE \\\n" .
                "--asAltExonSeriesAndTrsptIdList as.with.long..short.alt.trsptIdList.tsv \\\n" .
		"--asOriginTsv origin.tsv \\\n" .
		"--tissueSpecialAsTsv tissueSpecialHighPsi.tsv \\\n" .
		"--asProcessTrsptTsv as.process.trspt.mysql.tsv \\\n" .
		"--mirnaTargetTsv target.mysql.tsv \\\n" .
		"--orthAsTsv ../../9902/orthAs.tsv \\\n" .
		"--psiFileList jcec.A3SS.tsv,jcec.A5SS.tsv,jcec.MXE.tsv,jcec.RI.tsv,jcec.SE.tsv\\\n" .
		"--exptInfoTsv ../../010-gather-alignment-info-of-all-expts/filtered.alignment.info.of.assembled.experiment.tsv\\\n" .
		"--unitVolume 3\\\n" .
		"--minReadCovPerUnitVolume 5\\\n" .
		"--minValidedExptNum 3\\\n" .
		"--genePfamGoAnno gene.pfam.and.go.tsv \\\n" .
		"--outputAsMysqlTsv as.mysql.insert.tsv\n";
	exit;
}

my ($genomeFasta, $asCatalogFileList, $asTypeList, 
$genePfamGoAnno, $asAltExonSeriesAndTrsptIdList,
$dupConserAsTsv, $tissueSpecialAsTsv, $treatmentSpecialAsTsv, $asProcessTrsptTsv, $mirnaTargetTsv,
$tissueConstantAsTsv, $treatmentConstantAsTsv,
$psiFileList, $exptInfoTsv, $unitVolume, $minReadCovPerUnitVolume, $minValidedExptNum, $taxonId,
$orthAsTsv, $asOriginTsv, $outputAsMysqlTsv);

GetOptions(
        'genomeFasta=s'=>\$genomeFasta,
	'genePfamGoAnno=s'=>\$genePfamGoAnno,
        'asCatalogFileList=s'=>\$asCatalogFileList,
        'asTypeList=s'=>\$asTypeList,
        'asAltExonSeriesAndTrsptIdList=s'=>\$asAltExonSeriesAndTrsptIdList,
	'tissueSpecialAsTsv=s'=>\$tissueSpecialAsTsv,
	'asProcessTrsptTsv=s'=>\$asProcessTrsptTsv,
	'mirnaTargetTsv=s'=>\$mirnaTargetTsv,
	'taxonId=s'=>\$taxonId,
	'orthAsTsv=s'=>\$orthAsTsv,
	'asOriginTsv=s'=>\$asOriginTsv,
	'psiFileList=s'=>\$psiFileList,
	'exptInfoTsv=s'=>\$exptInfoTsv,
	'unitVolume=s'=>\$unitVolume,
	'minReadCovPerUnitVolume=s'=>\$minReadCovPerUnitVolume,
	'minValidedExptNum=s'=>\$minValidedExptNum,
	'outputAsMysqlTsv=s'=>\$outputAsMysqlTsv,
);


my ($asId);
# 将基因组序列读入到hash中
my (%genomeSeq, @tt, $line, $id);
open FF, "<$genomeFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$id = $1;
		@tt = ();
		@tt = split(/ /, $id);
		$id = $tt[0];
	}else{
		$genomeSeq{$id}.=$line;
	}
}
close FF;

# 打开基因功能注释pfam and go annotation
# geneId  pfamList        GOList
# AT5G16970       PF16884,PF00107 GO:0016491,GO:0055114
my (%geneAnno, $geneAnnoHref, $geneId, $pfamList, $goList);
$geneAnnoHref = \%geneAnno;
open FF, "<$genePfamGoAnno";
<FF>;
while($line=<FF>){
	chomp($line);
	($geneId, $pfamList, $goList) = ("", "", "");
	($geneId, $pfamList, $goList) = split(/\t/, $line);
	$geneAnnoHref->{$geneId}->{"pfamList"} = $pfamList;
	$geneAnnoHref->{$geneId}->{"goList"} = $goList;
}
close FF;


# 将AS位置及类型读入到hash中
my (%as, $asHref, @titleField, @valueField, @asFile, $asFile, @asType, $asType, $i, $j); 
$asHref = \%as;
@asFile = split(/,/, $asCatalogFileList);
@asType = split(/,/, $asTypeList);
# 打开每个AS文件
for($i=0; $i<=$#asFile; $i++){
	$asFile = $asFile[$i];
	$asType = $asType[$i];
	# 将该AS文件读入到hash中
	&readAsIntoHash($asFile, $asType, $asHref);
}





# 将实验信息读入
my (@fieldName, @fieldValue, $line, $exptId, $exptPos);
my (%exptInfo, $exptInfoHref);
$exptInfoHref=\%exptInfo;
open FF, "<$exptInfoTsv";
# Ecotype Cultivar        Genotype        Tissue  SubTissue       Development     TissueGroup     Treatment       TreatmentGroup  Experiment      Study
#    DataSource      Base    Layout  SpotsNum        ReadNum SpotLen ReadLen Gather  AS      Assemble        RunList Phenotype       alignPercent    mappedBases     mappedReadNum   detectedReadLen libraryType     phredScore
#    Ler     -       wt      anther  -       8-10WOld        anther_a        -       anther_a        SRX1892482      SRP075666       PRJNA322733     3.13 PAIRED  15.48   30.95   202     101     3       1       0       SRR3737052      -       91.45   2.86    28.30   101     UN      33
$line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
for($j=0; $j<=$#fieldName; $j++){
	if($fieldName[$j] eq "Experiment"){
		$exptPos = $j;
		last;
	}
}
while($line=<FF>){
	chomp($line);
	@fieldValue = split(/\t/, $line);
	# 将实验结果读入hash
	$exptId = $fieldValue[$exptPos];
	for(my $j=0; $j<=$#fieldValue; $j++){
		$exptInfoHref->{$exptId}->{$fieldName[$j]} = $fieldValue[$j];
	}
}



# 将psi值读入统计AS的inclusion和exclusion上分别被多上experiment验证
my (@psiFile, $psiFile, $asId, $mappedBase, $incluReadNum, $excluReadNum);
@psiFile = split(/,/, $psiFileList);
foreach $psiFile(@psiFile){
	open FF, "<$psiFile";
	<FF>;
	while($line=<FF>){
	# ASID    IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
	# ATHARI0000000001        77      89      535     100     SRX485073
	# ATHARI0000000002        225     0       211     100     SRX485073
		chomp($line);
		@fieldValue = split(/\t/, $line);
		$asId = $fieldValue[0];
		$exptId = $fieldValue[5];
		$mappedBase = $exptInfoHref->{$exptId}->{"mappedBases"};
		$incluReadNum = $fieldValue[1];
		$excluReadNum = $fieldValue[2];
		if($incluReadNum >= $mappedBase/$unitVolume*$minReadCovPerUnitVolume){
			$asHref->{$asId}->{"inclValidedExptNum"}++;
		}
		if($excluReadNum >= $mappedBase/$unitVolume*$minReadCovPerUnitVolume){
			$asHref->{$asId}->{"exclValidedExptNum"}++;
		}	
	}
	close FF;
}

# 将发生在重复基因上的保守AS读入
my ($asIdList, @asId);

### 将AS的两个外显子串位置及其对应的转录本编号等信息读入到hash中
my (%tmpAs, $tmpAsHref);
$tmpAsHref = \%tmpAs;
open FF, "<$asAltExonSeriesAndTrsptIdList";
# ASID chr strand longAltExonSeries longAltEnsemblTrsptIdList longAltImprovedTrsptIdList shortAltExonSeries shortAltEnsemblTrsptIdList shortAltImprovedTrsptIdList
$line = <FF>;
chomp($line);
@titleField = ();
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpAsHref->{$titleField[$i]} = $valueField[$i];
	}
	$asId = $tmpAsHref->{"ASID"};
	# 将AS的相关字段读入到总AS中
	$asHref->{$asId}->{"longAltExonSeries"} = $tmpAsHref->{"longAltExonSeries"};
	$asHref->{$asId}->{"longAltEnsemblTrsptIdList"} = $tmpAsHref->{"longAltEnsemblTrsptIdList"};
	$asHref->{$asId}->{"longAltImprovedTrsptIdList"} = $tmpAsHref->{"longAltImprovedTrsptIdList"};
	$asHref->{$asId}->{"shortAltExonSeries"} = $tmpAsHref->{"shortAltExonSeries"};
	$asHref->{$asId}->{"shortAltEnsemblTrsptIdList"} = $tmpAsHref->{"shortAltEnsemblTrsptIdList"};
	$asHref->{$asId}->{"shortAltImprovedTrsptIdList"} = $tmpAsHref->{"shortAltImprovedTrsptIdList"};

	# 生成 longAltTrsptIdList
	if($asHref->{$asId}->{"longAltEnsemblTrsptIdList"} ne "-" and $asHref->{$asId}->{"longAltImprovedTrsptIdList"} ne "-"){

		$asHref->{$asId}->{"longAltTrsptIdList"} = $asHref->{$asId}->{"longAltEnsemblTrsptIdList"} . "," . $asHref->{$asId}->{"longAltImprovedTrsptIdList"};

	}elsif($asHref->{$asId}->{"longAltEnsemblTrsptIdList"} ne "-" and $asHref->{$asId}->{"longAltImprovedTrsptIdList"} eq "-"){

		$asHref->{$asId}->{"longAltTrsptIdList"} = $asHref->{$asId}->{"longAltEnsemblTrsptIdList"};
		$asHref->{$asId}->{"longAltImprovedTrsptIdList"} = "NA";

	}elsif($asHref->{$asId}->{"longAltEnsemblTrsptIdList"} eq "-" and $asHref->{$asId}->{"longAltImprovedTrsptIdList"} ne "-"){

		$asHref->{$asId}->{"longAltTrsptIdList"} = $asHref->{$asId}->{"longAltImprovedTrsptIdList"};
		$asHref->{$asId}->{"longAltEnsemblTrsptIdList"} = "NA";

	}else{
		$asHref->{$asId}->{"longAltTrsptIdList"} = "NA";
		$asHref->{$asId}->{"longAltImprovedTrsptIdList"} = "NA";
		$asHref->{$asId}->{"longAltEnsemblTrsptIdList"} = "NA";
	}

	# 成 shortAltTrsptIdList
	if($asHref->{$asId}->{"shortAltEnsemblTrsptIdList"} ne "-" and $asHref->{$asId}->{"shortAltImprovedTrsptIdList"} ne "-"){

		$asHref->{$asId}->{"shortAltTrsptIdList"} = $asHref->{$asId}->{"shortAltEnsemblTrsptIdList"} . "," . $asHref->{$asId}->{"shortAltImprovedTrsptIdList"};

	}elsif($asHref->{$asId}->{"shortAltEnsemblTrsptIdList"} ne "-" and $asHref->{$asId}->{"shortAltImprovedTrsptIdList"} eq "-"){

		$asHref->{$asId}->{"shortAltImprovedTrsptIdList"} = "NA";
		$asHref->{$asId}->{"shortAltTrsptIdList"} = $asHref->{$asId}->{"shortAltEnsemblTrsptIdList"};

	}elsif($asHref->{$asId}->{"shortAltEnsemblTrsptIdList"} eq "-" and $asHref->{$asId}->{"shortAltImprovedTrsptIdList"} ne "-"){

		$asHref->{$asId}->{"shortAltEnsemblTrsptIdList"} = "NA";
		$asHref->{$asId}->{"shortAltTrsptIdList"} = $asHref->{$asId}->{"shortAltImprovedTrsptIdList"};

	}else{
		$asHref->{$asId}->{"shortAltTrsptIdList"} = "NA";
		$asHref->{$asId}->{"shortAltImprovedTrsptIdList"} = "NA";
		$asHref->{$asId}->{"shortAltEnsemblTrsptIdList"} = "NA";
	}
}
close FF;


# ########################
# 打开orthAs.tsv，将AS对应的orthId及保守水平读入到hash
my ($asIdPos, @valueField, $asIdList, @asId, $asId);
%tmpAs = ();
open FF, "<$orthAsTsv";
# orthAsId        	108875  13443   				225117  2711
# orthA3SS00000020      -       CARAA3SS0000000172,CARAA3SS0000007825	
$line=<FF>;
chomp($line);
@titleField = ();
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	$asId = $valueField[$asIdPos];
	for($i=0; $i<=$#valueField; $i++){
		$tmpAs{$titleField[$i]} = $valueField[$i];
	}
	# ASIdList
	if($tmpAs{$taxonId} ne "-"){
		$asIdList = $tmpAs{$taxonId};
		@asId = ();
		@asId = split(/,/, $asIdList);
		foreach $asId(@asId){
			$asHref->{$asId}->{"orthAsId"} = $tmpAs{"orthAsId"};
		}
	}
}
close FF;


######################################################
#
#  为AS生成longAltSeq和shortAltSeq
#
my @asId = keys(%as);
my ($exon1Seq, $exon2Seq, $exon3Seq, $exon4Seq, $editSiteInAltSeq);
foreach $asId(@asId){
	# 为longAlt和shortAlt分别生成相互转换操作：
	# editSiteInAltSeq、cutSizeInLongAlt、insertSizeInLongAlt、cutSizeInshortAlt、insertSizeInShortAlt
	# 按照一定规则操作，就可以将longAltSeq转换为shortAltSeq，或者从shortAltSeq转换为longAltSeq
	
	# A3SS
	if($asHref->{$asId}->{"asType"} eq "A3SS"){

		if($asHref->{$asId}->{"strand"} eq "+"){

			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - ($asHref->{$asId}->{"flankingES"} + 1) + 1));
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"longExonStart_0base"}, $asHref->{$asId}->{"longExonEnd"} - ($asHref->{$asId}->{"longExonStart_0base"} + 1) + 1));
			$asHref->{$asId}->{"editSiteInAltSeq"} = length($exon1Seq);
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = $asHref->{$asId}->{"shortES"} - ($asHref->{$asId}->{"longExonStart_0base"} + 1) + 1;
			$asHref->{$asId}->{"insertSizeInLongAlt"} = 0;

			# shortAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - ($asHref->{$asId}->{"flankingES"} + 1) + 1));
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"shortES"}, $asHref->{$asId}->{"shortEE"} - ($asHref->{$asId}->{"shortES"} + 1) + 1));			
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = 0;
			$asHref->{$asId}->{"insertSizeInShortAlt"} = $asHref->{$asId}->{"cutSizeInLongAlt"};
			
		}else{
			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - ($asHref->{$asId}->{"flankingES"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq =~tr/ACGT/TGCA/;
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"longExonStart_0base"}, $asHref->{$asId}->{"longExonEnd"} - ($asHref->{$asId}->{"longExonStart_0base"} + 1) + 1));
			$exon2Seq = reverse($exon2Seq);
			$exon2Seq =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq . $exon2Seq; #
			$asHref->{$asId}->{"editSiteInAltSeq"} = length($exon1Seq);
			$asHref->{$asId}->{"cutSizeInLongAlt"} = $asHref->{$asId}->{"longExonEnd"} - ($asHref->{$asId}->{"shortEE"} + 1) + 1;
			$asHref->{$asId}->{"insertSizeInLongAlt"} = 0;

			# shortAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - ($asHref->{$asId}->{"flankingES"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq =~tr/ACGT/TGCA/;
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"shortES"}, $asHref->{$asId}->{"shortEE"} - ($asHref->{$asId}->{"shortES"} + 1) + 1));
			$exon2Seq = reverse($exon2Seq);
			$exon2Seq =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = 0;
			$asHref->{$asId}->{"insertSizeInShortAlt"} = $asHref->{$asId}->{"cutSizeInLongAlt"};
	
		}	
	}

	# A5SS
	if($asHref->{$asId}->{"asType"} eq "A5SS"){

		if($asHref->{$asId}->{"strand"} eq "+"){

			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"longExonStart_0base"}, $asHref->{$asId}->{"longExonEnd"} - ($asHref->{$asId}->{"longExonStart_0base"} + 1) + 1));
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - ($asHref->{$asId}->{"flankingES"} + 1) + 1));
			$asHref->{$asId}->{"editSiteInAltSeq"} = $asHref->{$asId}->{"shortEE"} - ($asHref->{$asId}->{"shortES"} + 1) + 1;
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = $asHref->{$asId}->{"longExonEnd"} - ($asHref->{$asId}->{"shortEE"} + 1) + 1;
			$asHref->{$asId}->{"insertSizeInLongAlt"} = 0;

			# shortAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"shortES"}, $asHref->{$asId}->{"shortEE"} - ($asHref->{$asId}->{"shortES"} + 1) + 1));
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - ($asHref->{$asId}->{"flankingES"} + 1) + 1));
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = 0;
			$asHref->{$asId}->{"insertSizeInShortAlt"} = $asHref->{$asId}->{"cutSizeInLongAlt"};
			
		}else{
			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"longExonStart_0base"}, $asHref->{$asId}->{"longExonEnd"} - ($asHref->{$asId}->{"longExonStart_0base"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq =~tr/ACGT/TGCA/;
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - ($asHref->{$asId}->{"flankingES"} + 1) + 1));
			$exon2Seq = reverse($exon2Seq);
			$exon2Seq =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"editSiteInAltSeq"} = $asHref->{$asId}->{"shortEE"} - ($asHref->{$asId}->{"shortES"} + 1) + 1;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = $asHref->{$asId}->{"shortES"} - ($asHref->{$asId}->{"longExonStart_0base"} + 1) + 1;
			$asHref->{$asId}->{"insertSizeInLongAlt"} = 0;

			# shortAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"shortES"}, $asHref->{$asId}->{"shortEE"} - ($asHref->{$asId}->{"shortES"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq =~tr/ACGT/TGCA/;
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - ($asHref->{$asId}->{"flankingES"} + 1) + 1));
			$exon2Seq = reverse($exon2Seq);
			$exon2Seq =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = 0;
			$asHref->{$asId}->{"insertSizeInShortAlt"} = $asHref->{$asId}->{"cutSizeInLongAlt"};
		}
	}

	# RI
	if($asHref->{$asId}->{"asType"} eq "RI"){

		if($asHref->{$asId}->{"strand"} eq "+"){

			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"riExonStart_0base"}, $asHref->{$asId}->{"riExonEnd"} - ($asHref->{$asId}->{"riExonStart_0base"} + 1) + 1));
			$asHref->{$asId}->{"editSiteInAltSeq"} = $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1;
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = $asHref->{$asId}->{"downstreamES"} - ($asHref->{$asId}->{"upstreamEE"} + 1) + 1;
			$asHref->{$asId}->{"insertSizeInLongAlt"} = 0;

			# shortAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1));
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1));
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = 0;
			$asHref->{$asId}->{"insertSizeInShortAlt"} = $asHref->{$asId}->{"cutSizeInLongAlt"};
			
		}else{
			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"riExonStart_0base"}, $asHref->{$asId}->{"riExonEnd"} - ($asHref->{$asId}->{"riExonStart_0base"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"editSiteInAltSeq"} = $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1; # 
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = $asHref->{$asId}->{"downstreamES"} - ($asHref->{$asId}->{"upstreamEE"} + 1) + 1;
			$asHref->{$asId}->{"insertSizeInLongAlt"} = 0;

			# shortAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq =~tr/ACGT/TGCA/;
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1));
			$exon2Seq = reverse($exon2Seq);
			$exon2Seq =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = 0;
			$asHref->{$asId}->{"insertSizeInShortAlt"} = $asHref->{$asId}->{"cutSizeInLongAlt"};
		}
	}

	# SE
	if($asHref->{$asId}->{"asType"} eq "SE"){
		if($asHref->{$asId}->{"strand"} eq "+"){
			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1));
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"exonStart_0base"}, $asHref->{$asId}->{"exonEnd"} - ($asHref->{$asId}->{"exonStart_0base"} + 1) + 1));
			$exon3Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1));
			$asHref->{$asId}->{"editSiteInAltSeq"} = length($exon1Seq);
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq . $exon2Seq . $exon3Seq;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = length($exon2Seq);
			$asHref->{$asId}->{"insertSizeInLongAlt"} = 0;

			# shortAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1));
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1));
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = 0;
			$asHref->{$asId}->{"insertSizeInShortAlt"} = $asHref->{$asId}->{"cutSizeInLongAlt"};			
		}else{
			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq=~tr/ACGT/TGCA/;

			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"exonStart_0base"}, $asHref->{$asId}->{"exonEnd"} - ($asHref->{$asId}->{"exonStart_0base"} + 1) + 1));
			$exon2Seq = reverse($exon2Seq);
			$exon2Seq =~tr/ACGT/TGCA/;

			$exon3Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1));
			$exon3Seq = reverse($exon3Seq);
			$exon3Seq =~tr/ACGT/TGCA/;

			$asHref->{$asId}->{"editSiteInAltSeq"} = length($exon1Seq);
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq . $exon2Seq . $exon3Seq;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = length($exon2Seq);
			$asHref->{$asId}->{"insertSizeInLongAlt"} = 0;

			# shortAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq=~tr/ACGT/TGCA/;
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1));
			$exon2Seq = reverse($exon2Seq);
			$exon2Seq =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon2Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = 0;
			$asHref->{$asId}->{"insertSizeInShortAlt"} = $asHref->{$asId}->{"cutSizeInLongAlt"};
		}
	}

	# MXE
	if($asHref->{$asId}->{"asType"} eq "MXE"){

		if($asHref->{$asId}->{"strand"} eq "+"){
			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1));
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"1stExonStart_0base"}, $asHref->{$asId}->{"1stExonEnd"} - ($asHref->{$asId}->{"1stExonStart_0base"} + 1) + 1));
			$exon3Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"2ndExonStart_0base"}, $asHref->{$asId}->{"2ndExonEnd"} - ($asHref->{$asId}->{"2ndExonStart_0base"} + 1) + 1));
			$exon4Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1));
			$asHref->{$asId}->{"editSiteInAltSeq"} = length($exon1Seq);
			$asHref->{$asId}->{"longAltSeq"} = $exon1Seq . $exon2Seq . $exon4Seq;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = length($exon2Seq);
			$asHref->{$asId}->{"insertSizeInLongAlt"} = length($exon3Seq);;

			# shortAlt 
			$asHref->{$asId}->{"shortAltSeq"} = $exon1Seq . $exon3Seq . $exon4Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = length($exon3Seq);
			$asHref->{$asId}->{"insertSizeInShortAlt"} = length($exon2Seq);
			
		}else{
			# longAlt
			$exon1Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - ($asHref->{$asId}->{"upstreamES"} + 1) + 1));
			$exon1Seq = reverse($exon1Seq);
			$exon1Seq =~tr/ACGT/TGCA/;
			$exon2Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"1stExonStart_0base"}, $asHref->{$asId}->{"1stExonEnd"} - ($asHref->{$asId}->{"1stExonStart_0base"} + 1) + 1));
			$exon2Seq = reverse($exon2Seq);
			$exon2Seq =~tr/ACGT/TGCA/;
			$exon3Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"2ndExonStart_0base"}, $asHref->{$asId}->{"2ndExonEnd"} - ($asHref->{$asId}->{"2ndExonStart_0base"} + 1) + 1));
			$exon3Seq = reverse($exon3Seq);
			$exon3Seq =~tr/ACGT/TGCA/;
			$exon4Seq = uc(substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - ($asHref->{$asId}->{"downstreamES"} + 1) + 1));
			$exon4Seq = reverse($exon4Seq);
			$exon4Seq =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"editSiteInAltSeq"} = length($exon4Seq);
			$asHref->{$asId}->{"longAltSeq"} = $exon4Seq . $exon2Seq . $exon1Seq;
			$asHref->{$asId}->{"cutSizeInLongAlt"} = length($exon2Seq);
			$asHref->{$asId}->{"insertSizeInLongAlt"} = length($exon3Seq);;

			# shortAlt
			$asHref->{$asId}->{"shortAltSeq"} = $exon4Seq . $exon3Seq . $exon1Seq;
			$asHref->{$asId}->{"cutSizeInShortAlt"} = length($exon3Seq);
			$asHref->{$asId}->{"insertSizeInShortAlt"} = length($exon2Seq);
		}
	}
}

# 将AS注释来源读入到hash中，注释来源包括：ensemblGtf, improvedGtf, improvedGtf+RNAseqMapping
# ASID    ensemblGtf      improvedGtf     RNAseqMapping
# ATHAA3SS0000002088      0       0       1
my %tmpAs = ();
$tmpAsHref = \%tmpAs;
open FF, "<$asOriginTsv";
$line=<FF>;
chomp($line);
@titleField = ();
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	$asId = $valueField[0];
	for($i=0; $i<=$#valueField; $i++){
		$tmpAsHref->{$titleField[$i]} = $valueField[$i];
	}
	if($tmpAsHref->{"ensemblGtf"} == 1){
		$asHref->{$asId}->{"annoOrigin"} = "origAnnoGtf";
	}
	if($tmpAsHref->{"improvedGtf"} == 1){
		$asHref->{$asId}->{"annoOrigin"} = "improvedAnnoGtf";
	}
	if($tmpAsHref->{"RNAseqMapping"} == 1){
		$asHref->{$asId}->{"annoOrigin"} = "RNAseqMapping";
	}

}
close FF;

# 对AS进行归类，将相同gene的AS归置在一起提供给同组所有as
my (%geneToAsIdList, $geneToAsIdListHref);
$geneToAsIdListHref = \%geneToAsIdList;
my @asId = keys(%as);
foreach $asId(@asId){
	$geneToAsIdListHref->{$asHref->{$asId}->{"GeneID"}} .= $asId . ",";
}

# 读取组织特异的AS，将组织特异性赋给AS
# pvalue0.01      pvalue0.03      pvalue0.03      ASID    Regulate        Tissue=PSI      1stTissue|avgPSI|exptNum|pvalue 1stTissue|avgPSI|exptNum|pvalue ……
# N       N       N       ATHARI0000002764        High    endodermis=0.581|8      silique|0.250|21|0.000  stem|0.126|35|0.000     callus|0.227|5|0.000
#     seedling|0.178|267|0.000        inflorescence|0.145|41|0.000    leaf|0.155|268|0.000    ovule|0.035|8|0.000
my (@asTissuPsiField, $tissueName, $psiValue, $tissueSpecialType,  $avgPsiExptNum, $exptNum);
open FF, "<$tissueSpecialAsTsv";
<FF>;
while($line=<FF>){
	chomp($line);
	@asTissuPsiField = ();
	@asTissuPsiField = split(/\t/, $line);
	
	next if($asTissuPsiField[0] eq "N" and $asTissuPsiField[1] eq "N" and $asTissuPsiField[2] eq "N");

	$asId = $asTissuPsiField[3];
	$tissueSpecialType = $asTissuPsiField[4];

	($tissueName, $avgPsiExptNum) = split(/=/, $asTissuPsiField[5]);
	($psiValue, $exptNum) = split(/\|/, $avgPsiExptNum);
	$asHref->{$asId}->{"PsiSpecial" . $tissueSpecialType . "TissueName"} = $tissueName;
	$asHref->{$asId}->{"PsiSpecial" . $tissueSpecialType . "TissueValue"} = $psiValue;
}
close FF;


# 读取miRNA对mRNA的靶标，建立:
# trsptId->miRNA->{Target_start, Target_end}
my ($fieldList, $valueList, @value, %trsptIdToMirnaTargetPos, $trsptIdToMirnaTargetPosHref, $mirnaId, $trsptId);
$trsptIdToMirnaTargetPosHref = \%trsptIdToMirnaTargetPos;
open FF, "<$mirnaTargetTsv";
# miRNA_Acc., Target_Acc., Expectation, UPE, miRNA_start, miRNA_end, Target_start, Target_end, miRNA_aligned_fragment, alignment, Target_aligned_fragment, Inhibition, Target_Desc., Multiplicity___ath-miR156i, SRX399568.3733.4, 0.0, -1.0, 1, 20, 643, 662, UGACAGAAGAGAGAGAGCAG, ::::::::::::::::::::, CUGCUCUCUCUCUUCUGUCA, Cleavage, transcript_name:NA gene_id:AT1G53160 gene_name:NA, 1
while($line=<FF>){
	($fieldList, $valueList) = split(/___/, $line);
	@value = ();
	@value = split(/, /, $valueList);
	$mirnaId = $value[0];
	$trsptId = $value[1];
#	print "mirnaId:" . $mirnaId . "\t" . "\ttrsptId:" . $trsptId . "\n";
#	<STDIN>;
	$trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_start"} = $value[6];
	$trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_stop"} = $value[7];
}
close FF;

# 从asProcessTrspt中获得作用关系，建立：
# as->trsptI->{Target_start, Target_end}
# 检查trsptI上是否有miRNA的靶标，如果有有重叠，那么:
# 	$asHref->{$asId}->{"mrnaTargetList"}="miRNA1[trsptId1]|mirna2[trsptId2]"
my (%tmpAsTargetTrspt, $tmpAsTargetTrsptHref, @fieldName, @fieldValue, @mirnaId, $asTargetStart, $asTargetStop);
$tmpAsTargetTrsptHref = \%tmpAsTargetTrspt;
open FF, "<$asProcessTrsptTsv";
# asId, trsptId, residentType, asType, editSiteInOrigCdna, cutSize, insertSize, origOrfStart, origOrfStop, origStartCodonBegin, origStartCodonEnd, origStopCodonBegin, origStopCodonEnd, newOrfStart, newOrfStop, newStartCodonBegin, newStartCodonEnd, newStopCodonBegin, newStopCodonEnd, origGoTermList, origPfamList, newGoTermList, newPfamList, goImpact, pfamImpact, hitPfamList, frameImpact, stopCodonImpact, origCdnaSeq, origPepSeq, newCdnaSeq, newPepSeq, insertSeq___ATHAA3SS0000002088, SRX853408.3.6, ... 
while($line=<FF>){
	chomp($line);
	($fieldList, $valueList) = split(/___/, $line);
	@fieldName = split(/, /, $fieldList);;
	@fieldValue = split(/, /, $valueList);

	%tmpAsTargetTrspt = ();
	for(my $j=0; $j<=$#fieldName; $j++){
		$tmpAsTargetTrsptHref->{$fieldName[$j]} = $fieldValue[$j];
	}

	$asId = $tmpAsTargetTrsptHref->{"asId"};
	$trsptId = $tmpAsTargetTrsptHref->{"trsptId"};
#	print join("\t", "asId:" . $asId, "trsptId:" . $trsptId, "residentType:" . $tmpAsTargetTrsptHref->{"residentType"});
#	<STDIN>;
	# 根据as和trspt之间的作用关系分类
	if($tmpAsTargetTrsptHref->{"residentType"} eq "inclusion"){
		$asTargetStart = $tmpAsTargetTrsptHref->{"editSiteInOrigCdna"} + 1;
		$asTargetStop = $tmpAsTargetTrsptHref->{"editSiteInOrigCdna"} + $tmpAsTargetTrsptHref->{"cutSize"};
#		print join("\t", "residentType:inclusion", $asTargetStart, $asTargetStop);
#		<STDIN>;
		# 检测该转录本上是否存在mirna的靶标
		if(exists($trsptIdToMirnaTargetPosHref->{$trsptId})){
			# 获得该转录本上有mirna的靶标
			@mirnaId = ();
			@mirnaId = keys(%{$trsptIdToMirnaTargetPosHref->{$trsptId}});
#			print "mirna num:";
#			print $#mirnaId+1;
#			<STDIN>;
			# 获得该转录本上所有mirna
			foreach $mirnaId(@mirnaId){
#				print join("\t", $mirnaId, $trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_start"}, $trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_stop"});
#				<STDIN>;
				# 计算mirna的靶标位置是否和as的靶标位置重叠
				if(not(($asTargetStop < $trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_start"}) or ($asTargetStart > $trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_stop"}))){
					if(not(exists($asHref->{$asId}->{"mirnaTargetList"}))){
						$asHref->{$asId}->{"mirnaTargetList"} = $mirnaId . "[" . $trsptId . "]";
					}else{
						$asHref->{$asId}->{"mirnaTargetList"} .= "|" . $mirnaId . "[" . $trsptId . "]";
					}
				}
			}
		}
	}elsif($tmpAsTargetTrsptHref->{"residentType"} eq "exclusion"){
		$asTargetStart = $tmpAsTargetTrsptHref->{"editSiteInOrigCdna"} + 1;
		# 检测该转录本上是否存在mirna的靶标
		if(exists($trsptIdToMirnaTargetPosHref->{$trsptId})){
			# 获得转录本上有mirna的靶标
			@mirnaId = ();
			@mirnaId = keys(%{$trsptIdToMirnaTargetPosHref->{$trsptId}});
			# 获得该转录本上所有mirna
			foreach $mirnaId(@mirnaId){
				# 计算mirna的靶标位置是否和as的靶标位置重叠
				if($asTargetStart>=$trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_start"} and $asTargetStart<= $trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_stop"}){
					if(not(exists($asHref->{$asId}->{"mirnaTargetList"}))){
						$asHref->{$asId}->{"mirnaTargetList"} = $mirnaId . "[" . $trsptId . "]";
					}else{
						$asHref->{$asId}->{"mirnaTargetList"} .= "|" . $mirnaId . "[" . $trsptId . "]";
					}
				}
			}
		}
	}elsif($tmpAsTargetTrsptHref->{"residentType"} eq "1st" or $tmpAsTargetTrsptHref->{"residentType"} eq "2nd"){
		$asTargetStart = $tmpAsTargetTrsptHref->{"editSiteInOrigCdna"} + 1;
		$asTargetStart = $tmpAsTargetTrsptHref->{"editSiteInOrigCdna"} + $tmpAsTargetTrsptHref->{"cutSize"};

		# 检测该转录本上是否存在mirna的靶标
		if(exists($trsptIdToMirnaTargetPosHref->{$trsptId})){
			# 获得该转录本上有mirna的靶标
			@mirnaId = ();
			@mirnaId = keys(%{$trsptIdToMirnaTargetPosHref->{$trsptId}});
			# 获得该转录本上所有mirna
			foreach $mirnaId(@mirnaId){
				# 计算mirna的靶标位置是否和as的靶标位置重叠
				if(not(($asTargetStop < $trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_start"}) or ($asTargetStart > $trsptIdToMirnaTargetPosHref->{$trsptId}->{$mirnaId}->{"Target_stop"}))){
					if(not(exists($asHref->{$asId}->{"mirnaTargetList"}))){
						$asHref->{$asId}->{"mirnaTargetList"} = $mirnaId . "[" . $trsptId . "]";
					}else{
						$asHref->{$asId}->{"mirnaTargetList"} .= "|" . $mirnaId . "[" . $trsptId . "]";
					}
				}
			}
		}

	}
}
close FF;

# ===== 生成mysql需要的TSV文件 =====
open WW, ">$outputAsMysqlTsv";
my ($fieldString, $valueString, $start, $stop, $tissueSpecialTag, $treatmentSpecialTag);
@asId = keys(%as);
foreach $asId(@asId){

######### fieldString ########
	$fieldString = join(", ",
	"asId",
	"asType",
	"chr",
	"strand",
	"start",
	"stop",
	"orthAsGroupId",
	"asAnnoSource",
	"geneId", 
	"geneSymbol", 
	"genePfamList",
	"geneGoList",
	"otherAsIdListInSameGene",
	"1stExonEnd", "1stExonStart_0base", "2ndExonEnd", "2ndExonStart_0base", "downstreamEE", "downstreamES", "exonEnd", "exonStart_0base", "flankingEE", "flankingES", "longExonEnd", "longExonStart_0base", "riExonEnd", "riExonStart_0base", "shortEE", "shortES", "upstreamEE", "upstreamES", 
	"inclExonSeries",
	"inclTrsptIdList",
	"inclTrsptIdListInOrigAnnoGtf",
	"inclTrsptIdListInImprAnnoGtf",
	"exclExonSeries",
	"exclTrsptIdList",
	"exclTrsptIdListInOrigAnnoGtf",
	"excluTrsptIdListInImprAnnoGtf",
	"editSiteInBothAltSeq",
	"cutSizeInInclAlt",
	"insertSizeInInclAlt",
	"cutSizeInExclAlt",
	"insertSizeInExclAlt",
	"psiSpecialHighTissueName",
	"psiSpecialLowTissueName",
	"mirnaTargetList",
	"validedByEnoughExpts",
	"1stExon",
	"2ndExon",
	"upstreamExon",
	"downstreamExon",
	"flankingExon",
	"exon",
	"longExon",
	"shortExon",
	"riExon",
	"inclusionAltSeq",
	"exclusionAltSeq"
);

	if(not exists($asHref->{$asId}->{"orthAsId"})){
		$asHref->{$asId}->{"orthAsId"} = "NA";
	}
	if(not exists($asHref->{$asId}->{"orthAsConserLevel"})){
		$asHref->{$asId}->{"orthAsConserLevel"} = "NA";
	}

	# AS关联gene的pfam和GO
	if(not exists($geneAnnoHref->{$asHref->{$asId}->{"GeneID"}}->{"pfamList"}) or $geneAnnoHref->{$asHref->{$asId}->{"GeneID"}}->{"pfamList"} eq ""){
		$geneAnnoHref->{$asHref->{$asId}->{"GeneID"}}->{"pfamList"} = "NA";
	}
	if(not exists($geneAnnoHref->{$asHref->{$asId}->{"GeneID"}}->{"goList"}) or $geneAnnoHref->{$asHref->{$asId}->{"GeneID"}}->{"goList"} eq ""){
		$geneAnnoHref->{$asHref->{$asId}->{"GeneID"}}->{"goList"} = "NA";
	}

	# 获得AS在基因组上位置
	($start, $stop) = (0, 0);
	&getStartAndStop(\$start, \$stop, join(",", $asHref->{$asId}->{"1stExonEnd"}, $asHref->{$asId}->{"1stExonStart_0base"} + 1, $asHref->{$asId}->{"2ndExonEnd"}, $asHref->{$asId}->{"2ndExonStart_0base"} + 1, $asHref->{$asId}->{"downstreamEE"}, $asHref->{$asId}->{"downstreamES"} + 1, $asHref->{$asId}->{"exonEnd"}, $asHref->{$asId}->{"exonStart_0base"} + 1, $asHref->{$asId}->{"flankingEE"}, $asHref->{$asId}->{"flankingES"} + 1, $asHref->{$asId}->{"longExonEnd"}, $asHref->{$asId}->{"longExonStart_0base"} + 1, $asHref->{$asId}->{"riExonEnd"}, $asHref->{$asId}->{"riExonStart_0base"} + 1, $asHref->{$asId}->{"shortEE"}, $asHref->{$asId}->{"shortES"} + 1, $asHref->{$asId}->{"upstreamEE"}, $asHref->{$asId}->{"upstreamES"} + 1));

	# 获得相同相同基因下中其他AS列表
	$asHref->{$asId}->{"otherAsIdListInSameGene"} = $geneToAsIdListHref->{$asHref->{$asId}->{"GeneID"}};
	$asHref->{$asId}->{"otherAsIdListInSameGene"} = &removeAsIdFromAsIdList($asHref->{$asId}->{"otherAsIdListInSameGene"}, $asId);


	# 检测该AS是否为组织特异的
	if(not(exists($asHref->{$asId}->{"PsiSpecialHighTissueName"}))){
		$asHref->{$asId}->{"PsiSpecialHighTissueName"} = "NA";
		$asHref->{$asId}->{"PsiSpecialHighTissueValue"} = -1;
	}
	if(not(exists($asHref->{$asId}->{"PsiSpecialLowTissueName"}))){
		$asHref->{$asId}->{"PsiSpecialLowTissueName"} = "NA";
		$asHref->{$asId}->{"PsiSpecialLowTissueValue"} = -1;
	}

	# mirnaTargetList
	if(not(exists($asHref->{$asId}->{"mirnaTargetList"}))){
		$asHref->{$asId}->{"mirnaTargetList"} = "NA";
	}

	# 检测AS的inclusion和exclusion是否同时被足够数量的read覆盖而验证
	my $validedByEnoughExptTag = "No";
	if($asHref->{$asId}->{"inclValidedExptNum"} >= $minValidedExptNum and $asHref->{$asId}->{"exclValidedExptNum"} >= $minValidedExptNum){
		$validedByEnoughExptTag = "Yes";
	}

	# 1stExon
	if($asHref->{$asId}->{"1stExonEnd"} != -1){
		$asHref->{$asId}->{"1stExon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"1stExonStart_0base"}, $asHref->{$asId}->{"1stExonEnd"} - $asHref->{$asId}->{"1stExonStart_0base"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"1stExon"} = uc($asHref->{$asId}->{"1stExon"});
			$asHref->{$asId}->{"1stExon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"1stExon"} = reverse($asHref->{$asId}->{"1stExon"});
		}
	}else{
		$asHref->{$asId}->{"1stExon"} = "NA";
	}

	# 2ndExon
	if($asHref->{$asId}->{"2ndExonEnd"} != -1){
		$asHref->{$asId}->{"2ndExon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"2ndExonStart_0base"}, $asHref->{$asId}->{"2ndExonEnd"} - $asHref->{$asId}->{"2ndExonStart_0base"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"2ndExon"} = uc($asHref->{$asId}->{"2ndExon"});
			$asHref->{$asId}->{"2ndExon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"2ndExon"} = reverse($asHref->{$asId}->{"2ndExon"});
		}
	}else{
		$asHref->{$asId}->{"2ndExon"} = "NA";
	}

	# downstream
	if($asHref->{$asId}->{"downstreamEE"} != -1){
		$asHref->{$asId}->{"downstreamExon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"downstreamES"}, $asHref->{$asId}->{"downstreamEE"} - $asHref->{$asId}->{"downstreamES"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"downstreamExon"} = uc($asHref->{$asId}->{"downstreamExon"});
			$asHref->{$asId}->{"downstreamExon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"downstreamExon"} = reverse($asHref->{$asId}->{"downstreamExon"});
		}
	}else{
		$asHref->{$asId}->{"downstreamExon"} = "NA";
	}

	# upstream
	if($asHref->{$asId}->{"upstreamEE"} != -1){
		$asHref->{$asId}->{"upstreamExon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"upstreamES"}, $asHref->{$asId}->{"upstreamEE"} - $asHref->{$asId}->{"upstreamES"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"upstreamExon"} = uc($asHref->{$asId}->{"upstreamExon"});
			$asHref->{$asId}->{"upstreamExon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"upstreamExon"} = reverse($asHref->{$asId}->{"upstreamExon"});
		}
	}else{
		$asHref->{$asId}->{"upstreamExon"} = "NA";
	}

	# exon
	if($asHref->{$asId}->{"exonEnd"} != -1){
		$asHref->{$asId}->{"exon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"exonStart_0base"}, $asHref->{$asId}->{"exonEnd"} - $asHref->{$asId}->{"exonStart_0base"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"exon"} = uc($asHref->{$asId}->{"exon"});
			$asHref->{$asId}->{"exon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"exon"} = reverse($asHref->{$asId}->{"exon"});
		}
	}else{
		$asHref->{$asId}->{"exon"} = "NA";
	}

	# flanking
	if($asHref->{$asId}->{"flankingEE"} != -1){
		$asHref->{$asId}->{"flankingExon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"flankingES"}, $asHref->{$asId}->{"flankingEE"} - $asHref->{$asId}->{"flankingES"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"flankingExon"} = uc($asHref->{$asId}->{"flankingExon"});
			$asHref->{$asId}->{"flankingExon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"flankingExon"} = reverse($asHref->{$asId}->{"flankingExon"});
		}
	}else{
		$asHref->{$asId}->{"flankingExon"} = "NA";
	}

	# longExon
	if($asHref->{$asId}->{"longExonEnd"} != -1){
		$asHref->{$asId}->{"longExon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"longExonStart_0base"}, $asHref->{$asId}->{"longExonEnd"} - $asHref->{$asId}->{"longExonStart_0base"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"longExon"} = uc($asHref->{$asId}->{"longExon"});
			$asHref->{$asId}->{"longExon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"longExon"} = reverse($asHref->{$asId}->{"longExon"});
		}
	}else{
		$asHref->{$asId}->{"longExon"} = "NA";
	}

	# shortExon
	if($asHref->{$asId}->{"shortEE"} != -1){
		$asHref->{$asId}->{"shortExon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"shortES"}, $asHref->{$asId}->{"shortEE"} - $asHref->{$asId}->{"shortES"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"shortExon"} = uc($asHref->{$asId}->{"shortExon"});
			$asHref->{$asId}->{"shortExon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"shortExon"} = reverse($asHref->{$asId}->{"shortExon"});
		}
	}else{
		$asHref->{$asId}->{"shortExon"} = "NA";
	}

	# riExon
	if($asHref->{$asId}->{"riExonEnd"} != -1){
		$asHref->{$asId}->{"riExon"} = substr($genomeSeq{$asHref->{$asId}->{"chr"}}, $asHref->{$asId}->{"riExonStart_0base"}, $asHref->{$asId}->{"riExonEnd"} - $asHref->{$asId}->{"riExonStart_0base"});
		if($asHref->{$asId}->{"strand"} eq "-"){
			$asHref->{$asId}->{"riExon"} = uc($asHref->{$asId}->{"riExon"});
			$asHref->{$asId}->{"riExon"} =~tr/ACGT/TGCA/;
			$asHref->{$asId}->{"riExon"} = reverse($asHref->{$asId}->{"riExon"});
		}
	}else{
		$asHref->{$asId}->{"riExon"} = "NA";
	}


####  valueString #######
	$valueString = join(", ", 
	$asId, 
	$asHref->{$asId}->{"asType"}, 
	$asHref->{$asId}->{"chr"}, 
	$asHref->{$asId}->{"strand"},
	$start, 
	$stop,
	$asHref->{$asId}->{"orthAsId"},
	$asHref->{$asId}->{"annoOrigin"},
	$asHref->{$asId}->{"GeneID"}, 
	$asHref->{$asId}->{"geneSymbol"},
	$geneAnnoHref->{$asHref->{$asId}->{"GeneID"}}->{"pfamList"},
	$geneAnnoHref->{$asHref->{$asId}->{"GeneID"}}->{"goList"},
	$asHref->{$asId}->{"otherAsIdListInSameGene"},
	$asHref->{$asId}->{"1stExonEnd"}, $asHref->{$asId}->{"1stExonStart_0base"}, 
	$asHref->{$asId}->{"2ndExonEnd"}, $asHref->{$asId}->{"2ndExonStart_0base"},
	$asHref->{$asId}->{"downstreamEE"}, $asHref->{$asId}->{"downstreamES"},
	$asHref->{$asId}->{"exonEnd"}, $asHref->{$asId}->{"exonStart_0base"},
	$asHref->{$asId}->{"flankingEE"}, $asHref->{$asId}->{"flankingES"},
	$asHref->{$asId}->{"longExonEnd"}, $asHref->{$asId}->{"longExonStart_0base"},
	$asHref->{$asId}->{"riExonEnd"}, $asHref->{$asId}->{"riExonStart_0base"},
	$asHref->{$asId}->{"shortEE"}, $asHref->{$asId}->{"shortES"},
	$asHref->{$asId}->{"upstreamEE"}, $asHref->{$asId}->{"upstreamES"},
	$asHref->{$asId}->{"longAltExonSeries"}, 
	$asHref->{$asId}->{"longAltTrsptIdList"}, $asHref->{$asId}->{"longAltEnsemblTrsptIdList"}, $asHref->{$asId}->{"longAltImprovedTrsptIdList"}, 
	$asHref->{$asId}->{"shortAltExonSeries"},
	$asHref->{$asId}->{"shortAltTrsptIdList"}, $asHref->{$asId}->{"shortAltEnsemblTrsptIdList"}, $asHref->{$asId}->{"shortAltImprovedTrsptIdList"},
	$asHref->{$asId}->{"editSiteInAltSeq"},	
	$asHref->{$asId}->{"cutSizeInLongAlt"}, $asHref->{$asId}->{"insertSizeInLongAlt"},
	$asHref->{$asId}->{"cutSizeInShortAlt"}, $asHref->{$asId}->{"insertSizeInShortAlt"},
	$asHref->{$asId}->{"PsiSpecialHighTissueName"},
	$asHref->{$asId}->{"PsiSpecialLowTissueName"},
	$asHref->{$asId}->{"mirnaTargetList"},
	$validedByEnoughExptTag,
	$asHref->{$asId}->{"1stExon"},
	$asHref->{$asId}->{"2ndExon"},
	$asHref->{$asId}->{"upstreamExon"},
	$asHref->{$asId}->{"downstreamExon"},
	$asHref->{$asId}->{"flankingExon"},
	$asHref->{$asId}->{"exon"},
	$asHref->{$asId}->{"longExon"},
	$asHref->{$asId}->{"shortExon"},
	$asHref->{$asId}->{"riExon"},
	$asHref->{$asId}->{"longAltSeq"}, 
	$asHref->{$asId}->{"shortAltSeq"},
	);
	print WW $fieldString . "___" . $valueString . "\n";
}
close WW;

#######################################################
#
# 将一个AS文件中所有的坐标信息读入到hash中
# #####################################################
sub readAsIntoHash{
	my ($asFile, $asType, $href) = @_;
	my (@titleField, @valueField, $j, $line, $asId);
	# 打开AS文件
	open FF, "<$asFile";
	# 读取字段名称
	# ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
	$line = <FF>;
	chomp($line);
	@titleField = ();
	@titleField = split(/\t/, $line);
	while($line=<FF>){
		chomp($line);
		@valueField = ();
		@valueField = split(/\t/, $line);
		$asId = $valueField[0];
		# 将as的每个字段读入hash中
		for($j=0; $j<=$#valueField; $j++){
			if($valueField[$j]=~/"(.*)"/){
				$href->{$asId}->{$titleField[$j]} = $1;
			}else{
				$href->{$asId}->{$titleField[$j]} = $valueField[$j];
			}
		}
		# 设置AS的类型
		$href->{$asId}->{"asType"} = $asType;
		# 将无值的字段值设为-1，表示不存在
		if(not exists($href->{$asId}->{"1stExonEnd"})){
			$href->{$asId}->{"1stExonEnd"} = -1;
		}
		if(not exists($href->{$asId}->{"1stExonStart_0base"})){
			$href->{$asId}->{"1stExonStart_0base"} = -1;
		}
		if(not exists($href->{$asId}->{"2ndExonEnd"})){
			$href->{$asId}->{"2ndExonEnd"} = -1;
		}
		if(not exists($href->{$asId}->{"2ndExonStart_0base"})){
			$href->{$asId}->{"2ndExonStart_0base"} = -1;
		}
		if(not exists($href->{$asId}->{"downstreamEE"})){
			$href->{$asId}->{"downstreamEE"} = -1;
		}
		if(not exists($href->{$asId}->{"downstreamES"})){
			$href->{$asId}->{"downstreamES"} = -1;
		}
		if(not exists($href->{$asId}->{"exonEnd"})){
			$href->{$asId}->{"exonEnd"} = -1;
		}
		if(not exists($href->{$asId}->{"exonStart_0base"})){
			$href->{$asId}->{"exonStart_0base"} = -1;
		}
		if(not exists($href->{$asId}->{"flankingEE"})){
			$href->{$asId}->{"flankingEE"} = -1;
		}
		if(not exists($href->{$asId}->{"flankingES"})){
			$href->{$asId}->{"flankingES"} = -1;
		}
		if(not exists($href->{$asId}->{"longExonEnd"})){
			$href->{$asId}->{"longExonEnd"} = -1;
		}
		if(not exists($href->{$asId}->{"longExonStart_0base"})){
			$href->{$asId}->{"longExonStart_0base"} = -1;
		}
		if(not exists($href->{$asId}->{"riExonEnd"})){
			$href->{$asId}->{"riExonEnd"} = -1;
		}
		if(not exists($href->{$asId}->{"riExonStart_0base"})){
			$href->{$asId}->{"riExonStart_0base"} = -1;
		}
		if(not exists($href->{$asId}->{"shortEE"})){
			$href->{$asId}->{"shortEE"} = -1;
		}
		if(not exists($href->{$asId}->{"shortES"})){
			$href->{$asId}->{"shortES"} = -1;
		}
		if(not exists($href->{$asId}->{"upstreamEE"})){
			$href->{$asId}->{"upstreamEE"} = -1;
		}
		if(not exists($href->{$asId}->{"upstreamES"})){
			$href->{$asId}->{"upstreamES"} = -1;
		}

	}
	close FF;
}

# 从asIdList中剔除掉指定的asId，返回剩余的asId
sub removeAsIdFromAsIdList{
	my ($asIdList, $asId) = @_;
	my (@asId, $returnAsIdList, $tmpAsId);
	@asId = split(/,/, $asIdList);
	foreach $tmpAsId(@asId){
		$returnAsIdList .= $tmpAsId . "," if($tmpAsId ne $asId);
	}
	return substr($returnAsIdList, 0, length($returnAsIdList) - 1);
}

sub getStartAndStop{
	my ($start, $stop, $coordinateList) = @_;
	my (@coordinate, $coordinate);
	$$start = 100000000000;
	$$stop = -1;
	@coordinate = split(/,/, $coordinateList);
	foreach $coordinate(@coordinate){
		next if($coordinate <= 0);
		if($$start > $coordinate){
			$$start = $coordinate;
		}
		if($$stop < $coordinate){
			$$stop = $coordinate;
		}
	}
}


