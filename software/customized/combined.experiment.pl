#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--firstCuratedExptTsv \\\n" .
                "--secondPairedExptTsv \\\n" .
                "--secondSingleExptTsv \\\n" .
		"--outputNewAddedSecondExptTsv \\\n" .
		"--outputFinalAllExptTsv \n";
	exit;
}

my ($firstCuratedExptTsv, $secondPairedExptTsv, $secondSingleExptTsv, $outputFinalAllExptTsv, $outputNewAddedSecondExptTsv);

GetOptions(
        'firstCuratedExptTsv=s'=>\$firstCuratedExptTsv,
        'secondPairedExptTsv=s'=>\$secondPairedExptTsv,
        'secondSingleExptTsv=s'=>\$secondSingleExptTsv,
	'outputNewAddedSecondExptTsv=s'=>\$outputNewAddedSecondExptTsv,
	'outputFinalAllExptTsv=s'=>\$outputFinalAllExptTsv,
);
my (@nameField, @valueField, $line, %tmpHash, $exptId, $i);
my (%firstCuratedExptId, %firstExpt, $firstExptHref, %newAddedRun, $newAddedRunHref, %newAddedExpt, $newAddedExptHref);
$firstExptHref = \%firstExpt;
$newAddedRunHref = \%newAddedRun;
$newAddedExptHref = \%newAddedExpt;
#Ecotype Cultivar Genotype Tissue SubTissue Development Treatment Experiment Study Base Layout SpotsNum ReadNum SpotLen ReadLen Gather AS Assemble RunList Phenotype
#- - - - - - - SRX362644 SRP027561 0.4  SINGLE 7.4  7.4  50  50  3 1 0 SRR941836 -
# 将firstCuratedExptTsv中的ExptId登记下来
open FF, "<$firstCuratedExptTsv";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line = <FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	$firstCuratedExptId{$tmpHash{"Experiment"}} = 1;

	# 登记到hash
	$exptId = $tmpHash{"Experiment"};
	for($i=0; $i<=$#valueField; $i++){
		$firstExptHref->{$exptId}->{$nameField[$i]} = $valueField[$i];		
	}
	$firstExptHref->{$exptId}->{"SubTissue"} = "-";
	$firstExptHref->{$exptId}->{"treatmentGroup"} = "-";
	$firstExptHref->{$exptId}->{"dataSource"} = "-";
}
close FF;

# 扫描second的paired和single的实验，登记其信息到newAddedRun
my ($runId, $exptId);
open FF, "<$secondPairedExptTsv";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line = <FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}
	$runId = $tmpHash{"RunList"};
	$exptId = $tmpHash{"Experiment"};

	# 如果当前run所在的experiment已经在curated中了，那么放弃该run的登记
	next if(exists($firstCuratedExptId{$exptId}));

	# 将当前的run登记到newAddedRun哈希中
	for($i=0; $i<=$#valueField; $i++){
		$newAddedRunHref->{$runId}->{$nameField[$i]} = $valueField[$i];
	}
	# 将当前的runId登记到newAddedExpt中
	$newAddedExptHref->{$exptId}->{"RunList"} .= $runId . ",";
}
close FF;

open FF, "<$secondSingleExptTsv";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line = <FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	$runId = $tmpHash{"RunList"};
	$exptId = $tmpHash{"Experiment"};

	# 如果当前run所在的experiment已经在curated中了，那么放弃该run的登记
	next if(exists($firstCuratedExptId{$exptId}));

	# 将当前的run登记到newAddedRun哈希中
	for($i=0; $i<=$#valueField; $i++){
		$newAddedRunHref->{$runId}->{$nameField[$i]} = $valueField[$i];
	}
	# 将当前的runId登记到newAddedExpt中
	$newAddedExptHref->{$exptId}->{"RunList"} .= $runId . ",";
}
close FF;

# 重新扫描newAddedExpt中的RunList，执行如下操作：
# (1) 判断spotLen, readLen, layout是否一致，如果不一致，那么放弃该experiment
# (2) 如果readLen和layout一致，那么累加base、spotNum和readNum
#     其他信息：Experiment,Study,SpotLen,ReadLen照抄即可
my @exptId = keys(%newAddedExpt);
my (@runId, %SpotLen, %Layout, %ReadLen, @SpotLen, @Layout, @ReadLen, $totalBase, $totalSpotNum, $totalReadNum);
foreach $exptId(@exptId){
	# 去掉RunList的末尾逗号
	$newAddedExptHref->{$exptId}->{"RunList"} = substr($newAddedExptHref->{$exptId}->{"RunList"}, 0, length($newAddedExptHref->{$exptId}->{"RunList"}) - 1);

	# 检测当前experiment下的RUN的readLen、layout、spotLen是否一致
	%SpotLen = ();
	@SpotLen = ();
	%Layout = ();
	@Layout = ();
	%ReadLen = ();
	@ReadLen = ();

	@runId = ();
	@runId = split(/,/, $newAddedExptHref->{$exptId}->{"RunList"});
	($totalBase, $totalSpotNum, $totalReadNum) = (0, 0, 0);
	foreach $runId(@runId){
		$Layout{$newAddedRunHref->{$runId}->{"Layout"}} = 1;
		$SpotLen{$newAddedRunHref->{$runId}->{"SpotLen"}} = 1;
		$ReadLen{$newAddedRunHref->{$runId}->{"ReadLen"}} = 1;
		$totalBase += $newAddedRunHref->{$runId}->{"Base"};
		$totalSpotNum += $newAddedRunHref->{$runId}->{"SpotsNum"};
		$totalReadNum += $newAddedRunHref->{$runId}->{"ReadNum"};
	}
	@SpotLen = keys(%SpotLen);	
	@Layout = keys(%Layout);
	@ReadLen = keys(%ReadLen);

	if($#SpotLen!=0 or $#Layout!=0 or $#ReadLen!=0){
		delete($newAddedExpt{$exptId});
		next;
	}

 # Ecotype Cultivar Genotype Tissue SubTissue Development Treatment Experiment Study Base Layout SpotsNum ReadNum SpotLen ReadLen Gather AS Assemble RunList Phenotype 
	# 生成其他部分信息
	$newAddedExptHref->{$exptId}->{"Ecotype"} = "-";
	$newAddedExptHref->{$exptId}->{"Cultivar"} = "-";
	$newAddedExptHref->{$exptId}->{"Genotype"} = "-";
	$newAddedExptHref->{$exptId}->{"Tissue"} = "-";
	$newAddedExptHref->{$exptId}->{"SubTissue"} = "-";
	$newAddedExptHref->{$exptId}->{"Development"} = "-";
	$newAddedExptHref->{$exptId}->{"Treatment"} = "-";
	$newAddedExptHref->{$exptId}->{"treatmentGroup"} = "-";
	$newAddedExptHref->{$exptId}->{"Experiment"} = $exptId;
	$newAddedExptHref->{$exptId}->{"Study"} = $newAddedRunHref->{$runId[0]}->{"Study"};
	$newAddedExptHref->{$exptId}->{"dataSource"} = "-";
	$newAddedExptHref->{$exptId}->{"Base"} = $totalBase/1000/1000/1000;
	$newAddedExptHref->{$exptId}->{"Layout"} = $newAddedRunHref->{$runId[0]}->{"Layout"};
	$newAddedExptHref->{$exptId}->{"SpotsNum"} = $totalSpotNum;
	$newAddedExptHref->{$exptId}->{"ReadNum"} = $totalReadNum;
	$newAddedExptHref->{$exptId}->{"SpotLen"} = $newAddedRunHref->{$runId[0]}->{"SpotLen"};
	$newAddedExptHref->{$exptId}->{"ReadLen"} = $newAddedRunHref->{$runId[0]}->{"ReadLen"};
	$newAddedExptHref->{$exptId}->{"Gather"} = "3";
	$newAddedExptHref->{$exptId}->{"AS"} = "1";
	$newAddedExptHref->{$exptId}->{"Assemble"} = "0";
	$newAddedExptHref->{$exptId}->{"RunList"} = $newAddedExptHref->{$exptId}->{"RunList"};
	$newAddedExptHref->{$exptId}->{"Phenotype"} = "-";
}

# 输出新增experiment
open WW, ">$outputNewAddedSecondExptTsv";
open ALL, ">$outputFinalAllExptTsv";
print WW "Ecotype\tCultivar\tGenotype\tTissue\tSubTissue\tDevelopment\tTreatment\ttreatmentGroup\tExperiment\tStudy\tdataSource\tBase\tLayout\tSpotsNum\tReadNum\tSpotLen\tReadLen\tGather\tAS\tAssemble\tRunList\tPhenotype\n";
print ALL "Ecotype\tCultivar\tGenotype\tTissue\tSubTissue\tDevelopment\tTreatment\ttreatmentGroup\tExperiment\tStudy\tdataSource\tBase\tLayout\tSpotsNum\tReadNum\tSpotLen\tReadLen\tGather\tAS\tAssemble\tRunList\tPhenotype\n";

# 首先将first校正experiment输出到总文件
@exptId = ();
@exptId = keys(%firstExpt);
foreach $exptId(@exptId){
	print ALL join("\t", $firstExptHref->{$exptId}->{"Ecotype"}, $firstExptHref->{$exptId}->{"Cultivar"}, $firstExptHref->{$exptId}->{"Genotype"}, $firstExptHref->{$exptId}->{"Tissue"}, $firstExptHref->{$exptId}->{"SubTissue"}, $firstExptHref->{$exptId}->{"Development"}, $firstExptHref->{$exptId}->{"Treatment"}, $firstExptHref->{$exptId}->{"treatmentGroup"}, $firstExptHref->{$exptId}->{"Experiment"}, $firstExptHref->{$exptId}->{"Study"},  $firstExptHref->{$exptId}->{"dataSource"}, $firstExptHref->{$exptId}->{"Base"}, $firstExptHref->{$exptId}->{"Layout"}, $firstExptHref->{$exptId}->{"SpotsNum"}, $firstExptHref->{$exptId}->{"ReadNum"}, $firstExptHref->{$exptId}->{"SpotLen"}, $firstExptHref->{$exptId}->{"ReadLen"}, $firstExptHref->{$exptId}->{"Gather"}, $firstExptHref->{$exptId}->{"AS"}, $firstExptHref->{$exptId}->{"Assemble"}, $firstExptHref->{$exptId}->{"RunList"}, $firstExptHref->{$exptId}->{"Phenotype"}) . "\n";
}

# 把新增的experiment输出到新增实验文件中和总文件中
@exptId = ();
@exptId = keys(%newAddedExpt);
foreach $exptId(@exptId){
	print WW join("\t", $newAddedExptHref->{$exptId}->{"Ecotype"}, $newAddedExptHref->{$exptId}->{"Cultivar"}, $newAddedExptHref->{$exptId}->{"Genotype"}, $newAddedExptHref->{$exptId}->{"Tissue"}, $newAddedExptHref->{$exptId}->{"SubTissue"}, $newAddedExptHref->{$exptId}->{"Development"}, $newAddedExptHref->{$exptId}->{"Treatment"}, $newAddedExptHref->{$exptId}->{"treatmentGroup"}, $newAddedExptHref->{$exptId}->{"Experiment"}, $newAddedExptHref->{$exptId}->{"Study"},  $newAddedExptHref->{$exptId}->{"dataSource"}, $newAddedExptHref->{$exptId}->{"Base"}, $newAddedExptHref->{$exptId}->{"Layout"}, $newAddedExptHref->{$exptId}->{"SpotsNum"}, $newAddedExptHref->{$exptId}->{"ReadNum"}, $newAddedExptHref->{$exptId}->{"SpotLen"}, $newAddedExptHref->{$exptId}->{"ReadLen"}, $newAddedExptHref->{$exptId}->{"Gather"}, $newAddedExptHref->{$exptId}->{"AS"}, $newAddedExptHref->{$exptId}->{"Assemble"}, $newAddedExptHref->{$exptId}->{"RunList"}, $newAddedExptHref->{$exptId}->{"Phenotype"}) . "\n";
	print ALL join("\t", $newAddedExptHref->{$exptId}->{"Ecotype"}, $newAddedExptHref->{$exptId}->{"Cultivar"}, $newAddedExptHref->{$exptId}->{"Genotype"}, $newAddedExptHref->{$exptId}->{"Tissue"}, $newAddedExptHref->{$exptId}->{"SubTissue"}, $newAddedExptHref->{$exptId}->{"Development"}, $newAddedExptHref->{$exptId}->{"Treatment"}, $newAddedExptHref->{$exptId}->{"treatmentGroup"}, $newAddedExptHref->{$exptId}->{"Experiment"}, $newAddedExptHref->{$exptId}->{"Study"},  $newAddedExptHref->{$exptId}->{"dataSource"}, $newAddedExptHref->{$exptId}->{"Base"}, $newAddedExptHref->{$exptId}->{"Layout"}, $newAddedExptHref->{$exptId}->{"SpotsNum"}, $newAddedExptHref->{$exptId}->{"ReadNum"}, $newAddedExptHref->{$exptId}->{"SpotLen"}, $newAddedExptHref->{$exptId}->{"ReadLen"}, $newAddedExptHref->{$exptId}->{"Gather"}, $newAddedExptHref->{$exptId}->{"AS"}, $newAddedExptHref->{$exptId}->{"Assemble"}, $newAddedExptHref->{$exptId}->{"RunList"}, $newAddedExptHref->{$exptId}->{"Phenotype"}) . "\n";
}

close ALL;
close WW;
