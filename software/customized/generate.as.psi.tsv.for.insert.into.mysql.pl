#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--jcecFileList jcec.A3SS.txt,jcec.A5SS.txt,jcec.MXE.txt,jcec.RI.txt,jcec.SE.txt\\\n" . 
		"--jcFileList jc.A3SS.txt,jc.A5SS.txt,jc.MXE.txt,jc.RI.txt,jc.SE.txt\\\n" . 
		"--filteredExptInfoTsv \\\n" .
		"--unitVolum	3 \\\n".
		"--minReadCovPerUnitVolum 5 \\\n" .
		"--outputPsiTsv as.psi.mysql.tsv \n\n";
	exit;
}

my ($jcecFileList, $jcFileList, $unitVolum, $minReadCovPerUnitVolum, $filteredExptInfoTsv, $outputPsiTsv);

GetOptions(
        'jcecFileList=s'=>\$jcecFileList,
        'jcFileList=s'=>\$jcFileList,
        'filteredExptInfoTsv=s'=>\$filteredExptInfoTsv,
	'unitVolum=s'=>\$unitVolum,
	'minReadCovPerUnitVolum=s'=>\$minReadCovPerUnitVolum,
        'outputPsiTsv=s'=>\$outputPsiTsv,
);

my (@fieldName, @field, $i, $exptId, $exptIdPos, @jcecFile, @jcFile, $psiFile);
my ($experimentLine, $jcecLine, $jcLine, @titleField, @valueField, $j);

# 将过滤后的experiment编号登记到hash中
open FF, "<$filteredExptInfoTsv";
$experimentLine=<FF>;
chomp($experimentLine);
@fieldName = split(/\t/, $experimentLine);
my $exptIdPos = 0;
for($exptIdPos=0; $exptIdPos<=$#fieldName; $exptIdPos++){
        last if($fieldName[$exptIdPos] eq "Experiment");
}
my (%experiment, $experimentHref, $exptId);
$experimentHref = \%experiment;
while($experimentLine=<FF>){
        chomp($experimentLine);
        @field = ();
        @field = split(/\t/, $experimentLine);
	# 将每个实验的具体信息导入hash
	$exptId = $field[$exptIdPos];
	for($i=0; $i<=$#field; $i++){
		$experimentHref->{$exptId}->{$fieldName[$i]} = $field[$i]
	}

	# 额外添加有效read覆盖的最低阈值
	$experimentHref->{$exptId}->{"minReadCov"} = ($experimentHref->{$exptId}->{"mappedBases"}/$unitVolum)*$minReadCovPerUnitVolum;
}
close WW;
my (%psi, $psiHref, $line);

# 将整个jcec和jc文件读入到以AS+experiment为单位的hash中
my ($fieldString, $valueString);
my (@asIdExptId, $asIdExptId, $asId);
open WW, ">$outputPsiTsv";

@jcecFile = split(/,/, $jcecFileList);
@jcFile = split(/,/, $jcFileList);

for($i=0; $i<=4; $i++){

	# 将5种类型AS的jcec和jc读入到hash中	
	%psi = ();
	$psiHref=\%psi;

	# 将第i种类型AS的jcec读入
	open JCEC, "<$jcecFile[$i]";
	$line = <JCEC>;
	chomp($line);
	@titleField = split(/\t/, $line);
	while($line=<JCEC>){
		chomp($line);
		@valueField = split(/\t/, $line);
		# ASID    IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
		# ATHAA3SS0000000001      549     189     398     100     SRX485073
		$exptId = $valueField[5];
		next if(not exists($experimentHref->{$exptId}));
		$psiHref->{$valueField[0] . "#" . $exptId}->{"tissue"} = $experimentHref->{$exptId}->{"Tissue"};
		$psiHref->{$valueField[0] . "#" . $exptId}->{"treatment"} = $experimentHref->{$exptId}->{"Treatment"};
		$psiHref->{$valueField[0] . "#" . $exptId}->{"jcecIncl"} = $valueField[1];
		$psiHref->{$valueField[0] . "#" . $exptId}->{"jcecSkip"} = $valueField[2];
		$psiHref->{$valueField[0] . "#" . $exptId}->{"jcecIncFormLen"} = $valueField[3];
		$psiHref->{$valueField[0] . "#" . $exptId}->{"jcecSkipFormLen"} = $valueField[4];
		# 原来根据map的总量进行动态推算，现在在程序中固定死了，至少5个reads
		#if($valueField[1]+$valueField[2] < $experimentHref->{$valueField[5]}->{"minReadCov"} ){
		if($valueField[1]+$valueField[2] < 5 ){
			$psiHref->{$valueField[0] . "#" . $valueField[5]}->{"jcecPsivalid"} = "N";
		}else{
			$psiHref->{$valueField[0] . "#" . $valueField[5]}->{"jcecPsiValid"} = "Y";
		}
	}
	close JCEC;

	# 将第i种类型AS的jc读入
	open JC, "<$jcFile[$i]";
	$line = <JC>;
	chomp($line);
	@titleField = split(/\t/, $line);
	while($line=<JC>){
		chomp($line);
		@valueField = split(/\t/, $line);
		# ASID    IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
		# ATHAA3SS0000000001      549     189     398     100     SRX485073
		next if(not exists($experimentHref->{$valueField[5]}));
		$exptId = $valueField[5];
		$psiHref->{$valueField[0] . "#" . $exptId}->{"tissue"} = $experimentHref->{$exptId}->{"Tissue"};
		$psiHref->{$valueField[0] . "#" . $exptId}->{"treatment"} = $experimentHref->{$exptId}->{"Treatment"};
		$psiHref->{$valueField[0] . "#" . $valueField[5]}->{"jcIncl"} = $valueField[1];
		$psiHref->{$valueField[0] . "#" . $valueField[5]}->{"jcSkip"} = $valueField[2];
		$psiHref->{$valueField[0] . "#" . $valueField[5]}->{"jcIncFormLen"} = $valueField[3];
		$psiHref->{$valueField[0] . "#" . $valueField[5]}->{"jcSkipFormLen"} = $valueField[4];
		
		# 这里有效性验证从原来根据测序量动态调整，转为每个AS上固定至少10个reads覆盖
		#if($valueField[1]+$valueField[2] < $experimentHref->{$valueField[5]}->{"minReadCov"}){
		if($valueField[1]+$valueField[2] < 10){
			$psiHref->{$valueField[0] . "#" . $valueField[5]}->{"jcPsiValid"} = "N";
		}else{
			$psiHref->{$valueField[0] . "#" . $valueField[5]}->{"jcPsiValid"} = "Y";
		}
	}
	close JC;



	@asIdExptId = keys(%psi);
	foreach $asIdExptId(@asIdExptId){
		($asId, $exptId) = split(/#/, $asIdExptId);
############################################################
		## fieldString ###
		$fieldString = join(", ", 
			"asId", 
			"experimentId", 
			"tissue",
			"treatment",
			"jcecIncl", 
			"jcecSkip", 
			"jcecInclFormLen", 
			"jcecSkipFormLen", 
			"jcecPsi",
			"jcecPsiValid",
			"jcIncl", 
			"jcSkip", 
			"jcFormLen", 
			"jcSkipFormLen", 
			"jcPsi",
			"jcPsiValid"
		);

		if($psiHref->{$asIdExptId}->{"jcecIncl"} == 0 and $psiHref->{$asIdExptId}->{"jcecSkip"} != 0){
			$psiHref->{$asIdExptId}->{"jcecPsi"} = 0;
		}elsif($psiHref->{$asIdExptId}->{"jcecIncl"} != 0 and $psiHref->{$asIdExptId}->{"jcecSkip"} == 0){
			$psiHref->{$asIdExptId}->{"jcecPsi"} = 1;
		}elsif($psiHref->{$asIdExptId}->{"jcecIncl"} == 0 and $psiHref->{$asIdExptId}->{"jcecSkip"} == 0){
			$psiHref->{$asIdExptId}->{"jcecPsi"} = -1;
		}else{
			$psiHref->{$asIdExptId}->{"jcecPsi"} = sprintf("%.4f", ($psiHref->{$asIdExptId}->{"jcecIncl"}/$psiHref->{$asIdExptId}->{"jcecIncFormLen"}) / ($psiHref->{$asIdExptId}->{"jcecIncl"}/$psiHref->{$asIdExptId}->{"jcecIncFormLen"} + $psiHref->{$asIdExptId}->{"jcecSkip"}/$psiHref->{$asIdExptId}->{"jcecSkipFormLen"}));
		}

		if($psiHref->{$asIdExptId}->{"jcIncl"} == 0 and $psiHref->{$asIdExptId}->{"jcSkip"} != 0){
			$psiHref->{$asIdExptId}->{"jcPsi"} = 0;
		}elsif($psiHref->{$asIdExptId}->{"jcIncl"} != 0 and $psiHref->{$asIdExptId}->{"jcSkip"} == 0){
			$psiHref->{$asIdExptId}->{"jcPsi"} = 1;
		}elsif($psiHref->{$asIdExptId}->{"jcIncl"} == 0 and $psiHref->{$asIdExptId}->{"jcSkip"} == 0){
			$psiHref->{$asIdExptId}->{"jcPsi"} = -1;
		}else{
			$psiHref->{$asIdExptId}->{"jcPsi"} = sprintf("%.4f", ($psiHref->{$asIdExptId}->{"jcIncl"}/$psiHref->{$asIdExptId}->{"jcIncFormLen"}) / ($psiHref->{$asIdExptId}->{"jcIncl"}/$psiHref->{$asIdExptId}->{"jcIncFormLen"} + $psiHref->{$asIdExptId}->{"jcSkip"}/$psiHref->{$asIdExptId}->{"jcSkipFormLen"}));
		}

		# 用readCoverage检测psi值的有效性
		if(not exists($psiHref->{$asIdExptId}->{"jcecPsiValid"})){
			$psiHref->{$asIdExptId}->{"jcecPsiValid"} = "N";
		}
		if(not exists($psiHref->{$asIdExptId}->{"jcPsiValid"})){
			$psiHref->{$asIdExptId}->{"jcPsiValid"} = "N";
		}



############# valueString ############################################
		
		# 注释这里直接把Treatment设置为"NA"
		$psiHref->{$asIdExptId}->{"treatment"} = "NA";
		$valueString = join(", ",
			$asId,
			$exptId,
			$psiHref->{$asIdExptId}->{"tissue"},
			$psiHref->{$asIdExptId}->{"treatment"},
			$psiHref->{$asIdExptId}->{"jcecIncl"}, $psiHref->{$asIdExptId}->{"jcecSkip"},
			$psiHref->{$asIdExptId}->{"jcecIncFormLen"}, $psiHref->{$asIdExptId}->{"jcecSkipFormLen"},
			$psiHref->{$asIdExptId}->{"jcecPsi"},
			$psiHref->{$asIdExptId}->{"jcecPsiValid"},
			$psiHref->{$asIdExptId}->{"jcIncl"}, $psiHref->{$asIdExptId}->{"jcSkip"},
			$psiHref->{$asIdExptId}->{"jcIncFormLen"}, $psiHref->{$asIdExptId}->{"jcSkipFormLen"},
			$psiHref->{$asIdExptId}->{"jcPsi"},
			$psiHref->{$asIdExptId}->{"jcPsiValid"},
		);

		print WW $fieldString . "___" . $valueString . "\n";

	}# foreach
}
close WW;
