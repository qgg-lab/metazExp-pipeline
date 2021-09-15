#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--curatedExptTsv \\\n" .
		"--filterAssembledTrsptParameterFile \\\n" .
		"--jcecPsiFileList \\\n" .
		"--outputAssemblyDir \\\n" .
                "--outputExptTsv\n";
	exit;
}

my ($curatedExptTsv, $filterAssembledTrsptParameterFile, $jcecPsiFileList, $outputAssemblyDir, $outputExptTsv);

GetOptions(
        'curatedExptTsv=s'=>\$curatedExptTsv,
	'jcecPsiFileList=s'=>\$jcecPsiFileList,
	'filterAssembledTrsptParameterFile=s'=>\$filterAssembledTrsptParameterFile,
	'outputAssemblyDir=s'=>\$outputAssemblyDir,
        'outputExptTsv=s'=>\$outputExptTsv,
);

# 读取过滤trspt的参数值
my (%exptToFilterTrsptPara, $exptToFilterTrsptParaHref, $line, @valueField, $assemblyFile, $experimentId, $i);
my ($a3ssNum, $a5ssNum, $mxeNum, $riNum, $seNum);
$exptToFilterTrsptParaHref = \%exptToFilterTrsptPara;
open FF, "<$filterAssembledTrsptParameterFile";
# fileName   mappedBase minTrsptCov   minExonCov
# /mnt/home/liujind1/workAS1/3702/004-combine-assemblies-and-annos/../002-assemble-trsptome-on-goodExps/psiOutputDir/SRX3088664/transcriptomeByStringtie.gtf      7.23    3.62    1.81
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	if($line=~/.*psiOutputDir\/(.*)\/transcriptomeByStringtie\.gtf/){
		$experimentId = $1;
		$exptToFilterTrsptParaHref->{$experimentId}->{"assemblyFile"} = $valueField[0];
		$exptToFilterTrsptParaHref->{$experimentId}->{"minTrsptReadCov"} = $valueField[2];
		$exptToFilterTrsptParaHref->{$experimentId}->{"minExonReadCov"} = $valueField[3];		
	}
}
close FF;


my ($nameString, $valueString);
my ($a3ssNum, $a5ssNum, $mxeNum, $riNum, $seNum, $totalAsNum);
my (@nameField, @valueField, %tmpHash);
# Ecotype Cultivar        Genotype        Tissue  SubTissue       Development     TissueGroup     Treatment       TreatmentGroup  Experiment
#       Study   DataSource      Base    Layout  SpotsNum        ReadNum SpotLen ReadLen Gather  AS      Assemble        RunList Phenotype       alignPercent    mappedBases     mappedReadNum   detectedReadLen libraryType     phredScore

# 从010:filter中读取参与as识别的experiment
open FF, "<$curatedExptTsv";
open WW, ">$outputExptTsv";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$valueField[$i]=~s/ //;
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	$experimentId = $tmpHash{"Experiment"};

	# 在jcecPsiFileList中获得experimentId下5种类型AS的数量
	($totalAsNum, $a3ssNum, $a5ssNum, $mxeNum, $riNum, $seNum) = (0, 0, 0, 0, 0, 0);
	&getAsNum($experimentId, $jcecPsiFileList, \$a3ssNum, \$a5ssNum, \$mxeNum, \$riNum, \$seNum);
	$totalAsNum = $a3ssNum + $a5ssNum + $mxeNum + $riNum + $seNum;

	# 如果某个实验中的可变剪接数为0，那么放弃该实验
	next if($totalAsNum ==0 );

######### 生成experiment字段信息 ##############
	$nameString = join(", ", 
	"ecotype", 
	"cultivar", 
	"genotype", 
	"tissue", 
	"subTissue", 
	"development", 
	"tissueGroup", 
	"treatment", 
	"treatmentGroup", 
	"experimentId", 
	"studyId", 
	"dataSource", 
	"totalBasedVolume", 
	"layout", 
	"spotNum", 
	"readNum", 
	"splotLen", 
	"readLen", 
	"assemblyTag", 
	"runIdList", 
	"phenotype", 
	"alignPercentage", 
	"mappedBaseVolume", 
	"mappedReadNum", 
	"detectedReadLen", 
	"libraryType", 
	"phredScore", 
	"fileterTrsptMinTrsptReadCov", 
	"filterTrsptMinExonReadCov", 
	"asNum",
	"a3ssNum", 
	"a5ssNum", 
	"mxeNum", 
	"riNum", 
	"seNum"
);

	# 生成experiment的值信息
	if($tmpHash{"Assemble"} == 1 and (-e $exptToFilterTrsptParaHref->{$experimentId}->{"assemblyFile"})){
		#system("cp " . $exptToFilterTrsptParaHref->{$experimentId}->{"assemblyFile"} . " " . $outputAssemblyDir . "/" . $experimentId . ".gtf");
		$tmpHash{"Assemble"} = "Y";
	}else{
		$tmpHash{"Assemble"} = "N";
	}
	if(not exists($exptToFilterTrsptParaHref->{$experimentId}->{"minTrsptReadCov"})){
		$exptToFilterTrsptParaHref->{$experimentId}->{"minTrsptReadCov"} = -1;
		$exptToFilterTrsptParaHref->{$experimentId}->{"minExonReadCov"} = -1;
	}

	$tmpHash{"Ecotype"} = "NA" if($tmpHash{"Ecotype"} eq "-" or $tmpHash{"Ecotype"} eq "");
	$tmpHash{"Ecotype"}=~s/, /,_/g;
	$tmpHash{"Cultivar"} = "NA" if($tmpHash{"Cultivar"} eq "-" or $tmpHash{"Cultivar"} eq "");
	$tmpHash{"Cultivar"}=~s/, /,_/g;
	$tmpHash{"Genotype"} = "NA" if($tmpHash{"Genotype"} eq "-" or $tmpHash{"Genotype"} eq "");
	$tmpHash{"Genotype"}=~s/, /,_/g;
	$tmpHash{"Tissue"} = "NA" if($tmpHash{"Tissue"} eq "-" or $tmpHash{"Tissue"} eq "");
	$tmpHash{"Tissue"}=~s/, /,_/g;
	$tmpHash{"SubTissue"} = "NA" if($tmpHash{"SubTissue"} eq "-" or $tmpHash{"SubTissue"} eq "");
	$tmpHash{"SubTissue"}=~s/, /,_/g;
	$tmpHash{"SubTissue"}=~s/, /,_/g;
	$tmpHash{"Development"} = "NA" if($tmpHash{"Development"} eq "-" or $tmpHash{"Development"} eq "");
	$tmpHash{"Development"}=~s/, /,_/g;
	$tmpHash{"TissueGroup"} = "NA" if($tmpHash{"TissueGroup"} eq "-" or $tmpHash{"TissueGroup"} eq "");
	$tmpHash{"TissueGroup"}=~s/, /,_/g;
	$tmpHash{"Treatment"} = "NA" if($tmpHash{"Treatment"} eq "-" or $tmpHash{"Treatment"} eq "");
	$tmpHash{"Treatment"}=~s/, /,_/g;
	$tmpHash{"TreatmentGroup"} = "NA" if($tmpHash{"TreatmentGroup"} eq "-" or $tmpHash{"TreatmentGroup"} eq "");
	$tmpHash{"TreatmentGroup"}=~s/, /,_/g;
	$tmpHash{"DataSource"} = "NA" if($tmpHash{"DataSource"} eq "-" or $tmpHash{"DataSource"} eq "");
	$tmpHash{"DataSource"}=~s/, /,_/g;
	$tmpHash{"Phenotype"} = "NA" if($tmpHash{"Phenotype"} eq "-" or $tmpHash{"Phenotype"} eq "");
	$tmpHash{"Phenotype"}=~s/, /,_/g;
	$tmpHash{"libraryType"} = "NA" if($tmpHash{"libraryType"} eq "-" or $tmpHash{"libraryType"} eq "");
	$tmpHash{"libraryType"}=~s/, /,_/g;
	$valueString = join(", ", 
		$tmpHash{"Ecotype"},  
		$tmpHash{"Cultivar"},  
		$tmpHash{"Genotype"},  
		$tmpHash{"Tissue"},  
		$tmpHash{"SubTissue"},  
		$tmpHash{"Development"},  
		$tmpHash{"TissueGroup"},  
		$tmpHash{"Treatment"},  
		$tmpHash{"TreatmentGroup"},  
		$tmpHash{"Experiment"},  
		$tmpHash{"Study"},  
		$tmpHash{"DataSource"},  
		$tmpHash{"Base"},  
		$tmpHash{"Layout"},  
		$tmpHash{"SpotsNum"},  
		$tmpHash{"ReadNum"},  
		$tmpHash{"SpotLen"}, 
		$tmpHash{"ReadLen"}, 
		$tmpHash{"Assemble"}, 
		$tmpHash{"RunList"}, 
		$tmpHash{"Phenotype"}, 
		$tmpHash{"alignPercent"}, 
		$tmpHash{"mappedBases"}, 
		$tmpHash{"mappedReadNum"}, 
		$tmpHash{"detectedReadLen"}, 
		$tmpHash{"libraryType"}, 
		$tmpHash{"phredScore"}, 
		$exptToFilterTrsptParaHref->{$experimentId}->{"minTrsptReadCov"}, 
		$exptToFilterTrsptParaHref->{$experimentId}->{"minExonReadCov"}, 
		$totalAsNum,
		$a3ssNum,
		$a5ssNum, 
		$mxeNum, 
		$riNum, 
		$seNum
);
	
	# 输出1条experiment
	print WW $nameString . "___" . $valueString . "\n";
}
close FF;
close WW;

# 从总的psiFile中获得指定experiment的各类AS数量
# inclusion和skiping两者至少有1个read覆盖
sub getAsNum{
	my ($exptId, $psiFileList, $a3ssNum, $a5ssNum, $mxeNum, $riNum, $seNum) =@_;
	my (@asNum, $i, @psiFile, $psiFile, @line, $line, @field, $cmd, $mulLineText);
	
	@psiFile = split(/,/, $psiFileList);
	for($i=0; $i<=$#psiFile; $i++){
		$asNum[$i]=0;
		$cmd = "grep \"" . $exptId . "\" " . $psiFile[$i];
		$mulLineText = `$cmd`;
		@line = ();
		@line = split(/\n/, $mulLineText);
		foreach $line(@line){
			@field = split(/\t/, $line);
			next if($field[1] + $field[2] == 0);
			$asNum[$i]++;
		}
	}
	($$a3ssNum, $$a5ssNum, $$mxeNum, $$riNum, $$seNum) = ($asNum[0], $asNum[1], $asNum[2], $asNum[3], $asNum[4]);
}
