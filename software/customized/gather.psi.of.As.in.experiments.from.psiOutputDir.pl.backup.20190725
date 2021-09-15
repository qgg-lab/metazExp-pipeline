#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV < 0){
	print "This script is used to gather psi of AS in experiment and experiment sequencing information\n\n";
	print "\t perl $0 \\\n" . 
		"\t\t --psiOutputDir    ../008-pickup-psi-of-ASs-in-all-expers/psiOutputDir \\\n" . 
		"\t\t --expIdList       experiment.Id.list.txt \\\n" .
		"\t\t --speciesAbbr     ECAB \\\n" .
		"\t\t --filterMinReadLen 30 \\\n" .
		"\t\t --filterLibLayout SINGLE,PAIRED \\\n" .
		"\t\t --filterLibType U,R,F,UN,RF,FR \\\n" .
		"\t\t --filterAlignPer 50 \\\n" .
		"\t\t --filterTotalSpots 20 \\\n" .
		"\t\t --filterMappedSpots 20 \\\n" .

		"\t\t --outputExperimentInfo experimentInfo.tsv \\\n" . 

		"\t\t --outputTotalA5SS  totalA5SS.tsv \\\n" .
		"\t\t --outputTotalA3SS  totalA3SS.tsv \\\n" .
		"\t\t --outputTotalSE    totalSE.tsv \\\n" .
		"\t\t --outputTotalRI    totalRI.tsv \\\n" .
		"\t\t --outputTotalMXE   totalMXE.tsv \\\n" .

		"\t\t --outputJcecA5SS   jcecA5SS.tsv \\\n" .
		"\t\t --outputJcecA3SS   jcecA3SS.tsv \\\n" .
		"\t\t --outputJcecSE     jcecSE.tsv \\\n" .
		"\t\t --outputJcecRI     jcecRI.tsv \\\n" .
		"\t\t --outputJcecMXE    jcecMXE.tsv \\\n" .

		"\t\t --outputJcA5SS     jcA5SS.tsv \\\n" .
		"\t\t --outputJcA3SS     jcA3SS.tsv \\\n" .
		"\t\t --outputJcSE       jcSE.tsv \\\n" .
		"\t\t --outputJcRI       jcRI.tsv \\\n" .
		"\t\t --outputJcMXE      jcMXE.tsv \n";


	exit(0);
}

my ($speciesAbbr, $expId, @expDir, $psiOutputDir, $experimentId, $outputExperimentInfo);
my ($filterMinReadLen, $filterLibLayout, $filterLibType, $filterAlignPer, $filterTotalSpots, $filterMappedSpots);
my ($outputTotalA5SS, $outputTotalA3SS, $outputTotalSE, $outputTotalRI, $outputTotalMXE);
my ($outputNovelA5SS, $outputNovelA3SS, $outputNovelSE, $outputNovelRI, $outputNovelMXE);
my ($outputJcecA5SS, $outputJcecA3SS, $outputJcecSE, $outputJcecRI, $outputJcecMXE);
my ($outputJcA5SS, $outputJcA3SS, $outputJcSE, $outputJcRI, $outputJcMXE);

GetOptions(
	'speciesAbbr=s'=>\$speciesAbbr,
        'psiOutputDir=s'=>\$psiOutputDir,
        'expIdList=s'=>\$experimentId,
	'filterMinReadLen=s'=>\$filterMinReadLen, 
	'filterLibLayout=s'=>\$filterLibLayout,
	'filterLibType=s'=>\$filterLibType,
	'filterAlignPer=s'=>\$filterAlignPer,
	'filterTotalSpots=s'=>\$filterTotalSpots,
	'filterMappedSpots=s'=>\$filterMappedSpots,
        'outputExperimentInfo=s'=>\$outputExperimentInfo,
        'outputTotalA5SS=s'=>\$outputTotalA5SS,
        'outputTotalA3SS=s'=>\$outputTotalA3SS,
        'outputTotalSE=s'=>\$outputTotalSE,
        'outputTotalRI=s'=>\$outputTotalRI,
        'outputTotalMXE=s'=>\$outputTotalMXE,
        'outputJcecA5SS=s'=>\$outputJcecA5SS,
        'outputJcecA3SS=s'=>\$outputJcecA3SS,
        'outputJcecSE=s'=>\$outputJcecSE,
        'outputJcecRI=s'=>\$outputJcecRI,
        'outputJcecMXE=s'=>\$outputJcecMXE,
        'outputJcA5SS=s'=>\$outputJcA5SS,
        'outputJcA3SS=s'=>\$outputJcA3SS,
        'outputJcSE=s'=>\$outputJcSE,
        'outputJcRI=s'=>\$outputJcRI,
        'outputJcMXE=s'=>\$outputJcMXE,
);

# detect psiOutputDir
if(not -e $psiOutputDir){
	print STDERR "$psiOutputDir doesn't exist!\n";
	exit;
}

# process psiOutputDir to ensure this dir with "/" in the end.
my $expDirText = "";
$psiOutputDir = $psiOutputDir . "/" if(substr($psiOutputDir, length($psiOutputDir)-1, 1) ne "/");


# experiment directory will be extracted from experiment list file
if(uc($experimentId) ne "NONE" and -e $experimentId){
	# 抽取 psiOutputDir 目录下指定的experiment(不是全部experiment),指定的experiment的Id存放在$experimentId文件中
	open FF, "<$experimentId";
	while(my $expId=<FF>){
		chomp($expId);
		if(not -e $psiOutputDir . $expId . "/"){

			# 判断出指定的experiment对应的目录不存在
			print STDERR $psiOutputDir . $expId . " doesn't exist!\n";

		}else{
			
			# 如果指定的experiment存在，那么将experiment对应的全路径目录存放在expDir数组中
			$expDir[$#expDir+1] = $psiOutputDir . $expId;
			print STDOUT "Gather experiment information from " . $psiOutputDir . $expId . "/\n";

		}
	}
	close FF;	

}else{
	
	# if without specified expIdList file, experiment id will be detected from specified psiOutputDir
	# 没有通过文件指定抽取的experiment，那么遍历整个psiOutputDir目录,将每个experiment对应的结果目录存放在expDir数组中
	$expDirText = `find $psiOutputDir -type d -maxdepth 1`;
	@expDir = split(/\n/, $expDirText);
	shift(@expDir);
	for(my $i=0; $i<=$#expDir; $i++){
		print STDOUT "Gather experiment information from " . $expDir[$i] . "/\n";
	}
	
}


# gather information from experiment dir
my ($experimentDir);
my ($status, $runNum, $runId, $library, $layout, $phredScore);
my ($spotNum, $readLength, $alignPer);
my ($novelSpliceNum);
my ($annoCovGeneNum, $annoCovTranscriptNum);
my (%totalAS, %novelAS, %JCECAS, %JCAS);
my (@tmpArr, @tmpAsField, $asId);
my (%totalAsA5ssHash, %totalAsA3ssHash, %totalAsSeHash, %totalAsRiHash, %totalAsMxeHash);

#############################
# 收集信息产生16个文件
# (1) experimentInfo : A)测序和比对信息; B)实验中识别的各类总AS数量; C)novel AS的数量; D)jcec和jc登记的各类AS的数量(与AS总数是相等的)
# (2) 5个JCEC: ASID已经重新统一编过号，编号来自于(4)
# (3) 5个JC: ASID已经重新编过号，编号来自于(4)
# (4) 5个总的AS位置信息:
# ###########################
open EXPERIMENT, ">$outputExperimentInfo";
print EXPERIMENT join("\t", "expId", "status", "runNum", "runId", "library", "layout", "phredScore", "readLength",
                "spotNum(M)", "alignPer(%)", "novelSpliceNum",
                "totalA5SS", "totalA3SS", "totalSE", "totalRI", "totalMXE",
                "novelA5SS", "novelA3SS", "novelSE", "novelRI", "novelMXE",
                "jcecA5SS", "jcecA3SS", "jcecSE", "jcecRI", "jcecMXE",
                "jcA5SS", "jcA3SS", "jcSE", "jcRI", "jcMXE") . "\n";

open JCECA5SS, ">$outputJcecA5SS";
print JCECA5SS join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCECA3SS, ">$outputJcecA3SS";
print JCECA3SS join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCECSE, ">$outputJcecSE";
print JCECSE join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCECRI, ">$outputJcecRI";
print JCECRI join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCECMXE, ">$outputJcecMXE";
print JCECMXE join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCA5SS, ">$outputJcA5SS";
print JCA5SS join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCA3SS, ">$outputJcA3SS";
print JCA3SS join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCSE, ">$outputJcSE";
print JCSE join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCRI, ">$outputJcRI";
print JCRI join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

open JCMXE, ">$outputJcMXE";
print JCMXE join("\t", "ASID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncFormLen", "SkipFormLen", "Experiment") . "\n";

my (%catalogA5ss, %catalogA3ss, %catalogRi, %catalogMxe, %catalogSe);
my ($a5ssNum, $a3ssNum, $riNum, $seNum, $mxeNum);
my ($asPosition);


# scan each experiment to gather psi and experiment information
for(my $i=0; $i<=$#expDir; $i++){
	$runNum = 0;
	$runId = "";
	$library = "";
	$layout = "";
	$readLength = "";
	$phredScore = "";
	$spotNum = 0;
	$alignPer = "";
	$annoCovGeneNum = 0;
	$annoCovTranscriptNum = 0;
	
	%totalAS = ();
	%novelAS = ();
	%JCECAS = ();
	%JCAS =();

	@tmpArr = ();
	@tmpArr = split(/\//, $expDir[$i]);
	$expId = $tmpArr[$#tmpArr];

	$experimentDir = $expDir[$i] . "/";

	#detect running status
	if(-e $experimentDir . "/" . "fromGTF.A5SS.txt"){
		$status = "OK";
	}else{
		$status = "ERROR";
	}

	#detect novel splice num
	$novelSpliceNum = "";
	if(-e $experimentDir . "/" . $expId . ".novel.splicesite.from.hisat.txt"){
		my $cmd = "wc -l " . $experimentDir . "/" . $expId . ".novel.splicesite.from.hisat.txt";
		$novelSpliceNum = `$cmd`;
		chomp($novelSpliceNum);
		if($novelSpliceNum=~/(\d+) .*/){
			$novelSpliceNum = $1;
		}
	}

	#extract run information
	if(-e $experimentDir . "/" . $expId . ".SeqInfo.txt"){
		open FF, "<" . $experimentDir . "/" . $expId . ".SeqInfo.txt";		
		<FF>;
		while(my $line = <FF>){
			chomp($line);
			@tmpArr = ();
			@tmpArr = split(/\t/, $line);		
			if($#tmpArr==7){
				$runNum = $runNum + 1;
				$runId .= $tmpArr[1] . ",";
				$layout .= $tmpArr[2] . ",";
				$readLength .= $tmpArr[3] . ",";
				$phredScore .= $tmpArr[4] . ",";
				$library .= $tmpArr[6] . ",";
			}elsif($line=~/TotalSpot:(.*)\tAlignPercent:(.*)/){
				$spotNum = int(($1/1000000)*100)/100;
				$alignPer = $2;
			}
		}
		close FF;
		if($runNum > 0){
			$runId = substr($runId, 0, length($runId) -1);
			$layout =substr($layout, 0, length($layout) -1);
			$readLength = substr($readLength, 0, length($readLength) -1);
			$library = substr($library, 0, length($library) - 1);
			$phredScore = substr($phredScore, 0, length($phredScore) -1 );
		}
	}

# filter experiment with specified conditions
	next if($status eq "ERROR");

	my $tmpAlignPer = substr($alignPer, 0, length($alignPer)-1);
	#print join("\t", $readLength, $layout, $library, $spotNum, $tmpAlignPer);
	#<STDIN>;
	#print join("\t", $filterMinReadLen, $filterLibLayout, $filterLibType, $filterTotalSpots, $filterMappedSpots, $filterAlignPer);
	#<STDIN>;

	next if(&checkMinReadLen($readLength, $filterMinReadLen)==0);
	next if(&checkLibLayout($layout, $filterLibLayout)==0);
	next if(&checkLibType($library, $filterLibType)==0);

#	print join("\t", $spotNum , $filterTotalSpots);
#	<STDIN>;
	next if($spotNum < $filterTotalSpots);

#	print join("\t", $spotNum*$tmpAlignPer/100, $filterMappedSpots);
#	<STDIN>;
	next if($spotNum*$tmpAlignPer/100 < $filterMappedSpots);

#	print join("\t", $tmpAlignPer, $filterAlignPer);
#	<STDIN>;
	next if($tmpAlignPer < $filterAlignPer);

		
	#register total AS and extract total AS num
	%totalAsA5ssHash = ();
	%totalAsA3ssHash = ();
	%totalAsSeHash = ();
	%totalAsRiHash = ();
	%totalAsMxeHash = ();
	if(-e $experimentDir . "/" . "fromGTF.A5SS.txt"){

		############################################################################
		#									   #
		# 1. gather experiment information of sequencing, mapping and As num       #
		#									   #
		############################################################################

		my $cmd = "wc -l " . $experimentDir . "/" . "fromGTF.A5SS.txt";
		$totalAS{"A5SS"} = `$cmd`;
		chomp($totalAS{"A5SS"});
		if($totalAS{"A5SS"}=~/(\d+) .*/){
			$totalAS{"A5SS"} = $1 - 1;
		}
	
		$cmd = "wc -l " . $experimentDir . "/" . "fromGTF.A3SS.txt";
		$totalAS{"A3SS"} = `$cmd`;
		chomp($totalAS{"A3SS"});
		if($totalAS{"A3SS"}=~/(\d+) .*/){
			$totalAS{"A3SS"} = $1 - 1;
		}
	
		$cmd = "wc -l " . $experimentDir . "/" . "fromGTF.SE.txt";
		$totalAS{"SE"} = `$cmd`;
		chomp($totalAS{"SE"});
		if($totalAS{"SE"}=~/(\d+) .*/){
			$totalAS{"SE"} = $1 - 1;
		}

		$cmd = "wc -l " . $experimentDir . "/" . "fromGTF.RI.txt";
		$totalAS{"RI"} = `$cmd`;
		chomp($totalAS{"RI"});
		if($totalAS{"RI"}=~/(\d+) .*/){
			$totalAS{"RI"} = $1 -1;
		}

		$cmd = "wc -l " . $experimentDir . "/" . "fromGTF.MXE.txt";
		$totalAS{"MXE"} = `$cmd`;
		chomp($totalAS{"MXE"});
		if($totalAS{"MXE"}=~/(\d+) .*/){
			$totalAS{"MXE"} = $1 - 1;
		}

		############################################################################
		#									   #
		# 2. register AS into %catalogA5ss, %catalogA3ss ..., %catalogSe, 	   #
		#    as well as into %totalA5ss, totalA3SS ..., %totalSE.		   #
		#    For %catalogXx, position --> uniqId. such as 			   #
		#    "geneId#chr1#345612#23443..." --> A5SS0000001			   #
		#    For %totalA5ss, inner Id --> position
		#    12 --> "geneId#chr1#345612#23443..."
		#									   #
		############################################################################

		# read AS and save it into %catalogA5ss and %totalAsA5ssHash
		open FA5SS, "<" . $experimentDir . "/" . "fromGTF.A5SS.txt";
		<FA5SS>;
		while(my $line = <FA5SS>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			$asPosition = join("#", @tmpAsField);
			# save into %totalA5ss
			$totalAsA5ssHash{$asId} = $asPosition;
			# save into %catalogA5ss		
			if(not exists($catalogA5ss{$asPosition})){
				$a5ssNum++;
				$catalogA5ss{$asPosition}= $speciesAbbr . "A5SS" . sprintf("%010d", $a5ssNum);
			}	
		}
		close FA5SS;


		open FA3SS, "<" . $experimentDir . "/" . "fromGTF.A3SS.txt";
		<FA3SS>;
		while(my $line = <FA3SS>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			$asPosition = join("#", @tmpAsField);
			# save into %totalA3ss
			$totalAsA3ssHash{$asId} = $asPosition;
			# save into %catalogA3ss
			if(not exists($catalogA3ss{$asPosition})){
				$a3ssNum++;
				$catalogA3ss{$asPosition}= $speciesAbbr . "A3SS" . sprintf("%010d", $a3ssNum);
			}
		}
		close FA3SS;

		open FSE, "<" . $experimentDir . "/" . "fromGTF.SE.txt";
		<FSE>;
		while(my $line = <FSE>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			$asPosition = join("#", @tmpAsField);
			# save into %totalSE
			$totalAsSeHash{$asId} = $asPosition;
			# save into %catalogSe
			if(not exists($catalogSe{$asPosition})){
				$seNum++;
				$catalogSe{$asPosition}= $speciesAbbr . "SE" . sprintf("%010d", $seNum);
			}

		}
		close FSE;

		open FRI, "<" . $experimentDir . "/" . "fromGTF.RI.txt";
		<FRI>;
		while(my $line = <FRI>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			$asPosition = join("#", @tmpAsField);
			# save into %totalRi
			$totalAsRiHash{$asId} = $asPosition;
			# save into %catalogRi
			if(not exists($catalogRi{$asPosition})){
				$riNum++;
				$catalogRi{$asPosition}= $speciesAbbr . "RI" . sprintf("%010d", $riNum);
			}

		}
		close FRI;

		open FMXE, "<" . $experimentDir . "/" . "fromGTF.MXE.txt";
		<FMXE>;
		while(my $line = <FMXE>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			$asPosition = join("#", @tmpAsField);
			# save into %totalMxe
			$totalAsMxeHash{$asId} = $asPosition;
			# save into %catalogMxe
			if(not exists($catalogMxe{$asPosition})){
				$mxeNum++;
				$catalogMxe{$asPosition}= $speciesAbbr . "MXE" . sprintf("%010d", $mxeNum);
			}

		}
		close FMXE;
	}

	# obtain novel AS num
	if(-e $experimentDir . "/" . "fromGTF.novelEvents.A5SS.txt"){

		$novelAS{"A5SS"} = `cat $experimentDir/fromGTF.novelEvents.A5SS.txt | wc -l`;
		$novelAS{"A5SS"} = $novelAS{"A5SS"} - 1;

		$novelAS{"A3SS"} = `cat $experimentDir/fromGTF.novelEvents.A3SS.txt | wc -l`;
		$novelAS{"A3SS"} = $novelAS{"A3SS"} - 1;

		$novelAS{"SE"} = `cat $experimentDir/fromGTF.novelEvents.SE.txt | wc -l`;
		$novelAS{"SE"} = $novelAS{"SE"} - 1;

		$novelAS{"RI"} = `cat $experimentDir/fromGTF.novelEvents.RI.txt | wc -l`;
		$novelAS{"RI"} = $novelAS{"RI"} - 1;

		$novelAS{"MXE"} = `cat $experimentDir/fromGTF.novelEvents.MXE.txt | wc -l`;
		$novelAS{"MXE"} = $novelAS{"MXE"} - 1;
	
	}

	############################################################################
	#									   #
	# 3. scan AS to obtain inner Id and psi values then output glob Id and psi #
	# values 	   							   #
	#									   #
	############################################################################

        #extract AS num in JCEC
	if(-e $experimentDir . "/" . "JCEC.raw.input.A5SS.txt"){
		$JCECAS{"A5SS"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.A5SS.txt");
		$JCECAS{"A3SS"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.A3SS.txt");
		$JCECAS{"SE"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.SE.txt");
		$JCECAS{"RI"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.RI.txt");
		$JCECAS{"MXE"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.MXE.txt");

		#read psi of A5SS from JECE and write psi into output file
		open FA5SS, "<" . $experimentDir . "/" . "JCEC.raw.input.A5SS.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		<FA5SS>;
		while(my $line = <FA5SS>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCECA5SS join("\t", $catalogA5ss{$totalAsA5ssHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FA5SS;

		open FA3SS, "<" . $experimentDir . "/" . "JCEC.raw.input.A3SS.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		<FA3SS>;
		while(my $line = <FA3SS>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCECA3SS join("\t", $catalogA3ss{$totalAsA3ssHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FA3SS;

		open FSE, "<" . $experimentDir . "/" . "JCEC.raw.input.SE.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		<FSE>;
		while(my $line = <FSE>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCECSE join("\t", $catalogSe{$totalAsSeHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FSE;

		open FRI, "<" . $experimentDir . "/" . "JCEC.raw.input.RI.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       0       0                       102     99
		<FRI>;
		while(my $line = <FRI>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCECRI join("\t", $catalogRi{$totalAsRiHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FRI;

		open FMXE, "<" . $experimentDir . "/" . "JCEC.raw.input.MXE.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       0       0                       102     99
		<FMXE>;
		while(my $line = <FMXE>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCECMXE join("\t", $catalogMxe{$totalAsMxeHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FMXE;
	}
        
	#extract AS num in JC
	if(-e $experimentDir . "/" . "JC.raw.input.A5SS.txt"){
		$JCAS{"A5SS"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.A5SS.txt");
		$JCAS{"A3SS"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.A3SS.txt");
		$JCAS{"SE"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.SE.txt");
		$JCAS{"RI"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.RI.txt");
		$JCAS{"MXE"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.MXE.txt");

		open FA5SS, "<" . $experimentDir . "/" . "JC.raw.input.A5SS.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       0       0                       102     99
		<FA5SS>;
		while(my $line = <FA5SS>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCA5SS join("\t", $catalogA5ss{$totalAsA5ssHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FA5SS;

		open FA3SS, "<" . $experimentDir . "/" . "JC.raw.input.A3SS.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       0       0                       102     99
		<FA3SS>;
		while(my $line = <FA3SS>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCA3SS join("\t", $catalogA3ss{$totalAsA3ssHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FA3SS;

		open FSE, "<" . $experimentDir . "/" . "JC.raw.input.SE.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       0       0                       102     99
		<FSE>;
		while(my $line = <FSE>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCSE join("\t", $catalogSe{$totalAsSeHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FSE;

		open FRI, "<" . $experimentDir . "/" . "JC.raw.input.RI.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       0       0                       102     99
		<FRI>;
		while(my $line = <FRI>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCRI join("\t", $catalogRi{$totalAsRiHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FRI;

		open FMXE, "<" . $experimentDir . "/" . "JC.raw.input.MXE.txt";
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       0       0                       102     99
		<FMXE>;
		while(my $line = <FMXE>){
			chomp($line);
			@tmpAsField = ();
			@tmpAsField = split(/\t/, $line);		
			$asId = $tmpAsField[0];
			shift(@tmpAsField);
			print JCMXE join("\t", $catalogMxe{$totalAsMxeHash{$asId}}, $tmpAsField[0], $tmpAsField[1], $tmpAsField[4], $tmpAsField[5], $expId) . "\n";
		}
		close FMXE;

	}

	#output experiment information
	print EXPERIMENT join("\t", $expId, $status, $runNum, $runId, $library, $layout, $phredScore, $readLength, 
		$spotNum, $alignPer, $novelSpliceNum,
		$totalAS{"A5SS"}, $totalAS{"A3SS"}, $totalAS{"SE"}, $totalAS{"RI"}, $totalAS{"MXE"},
		$novelAS{"A5SS"}, $novelAS{"A3SS"}, $novelAS{"SE"}, $novelAS{"RI"}, $novelAS{"MXE"},
		$JCECAS{"A5SS"}, $JCECAS{"A3SS"}, $JCECAS{"SE"}, $JCECAS{"RI"}, $JCECAS{"MXE"},
		$JCAS{"A5SS"}, $JCAS{"A3SS"}, $JCAS{"SE"}, $JCAS{"RI"}, $JCAS{"MXE"}) . "\n";
}
close EXPERIMENT;
close JCECA5SS;
close JCECA3SS;
close JCECSE;
close JCECRI;
close JCECMXE;
close JCA5SS;
close JCA3SS;
close JCSE;
close JCRI;
close JCMXE;


# output total AS
open TOTALA5SS, ">$outputTotalA5SS";
print TOTALA5SS join("\t", "ASID", "GeneID", "geneSymbol", "chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE") . "\n";

open TOTALA3SS, ">$outputTotalA3SS";
print TOTALA3SS join("\t", "ASID", "GeneID", "geneSymbol", "chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE") . "\n";

open TOTALSE, ">$outputTotalSE";
print TOTALSE join("\t", "ASID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE") . "\n";

open TOTALRI, ">$outputTotalRI";
print TOTALRI join("\t", "ASID", "GeneID", "geneSymbol", "chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES","downstreamEE") . "\n";

open TOTALMXE, ">$outputTotalMXE";
print TOTALMXE join("\t", "ASID", "GeneID", "geneSymbol", "chr", "strand", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE") . "\n";

my (@asPosition, @tttt);
@asPosition = ();
@asPosition = keys(%catalogA5ss);
@asPosition = sort@asPosition;
foreach $asPosition(@asPosition){
	@tttt = ();
	@tttt = split(/#/, $asPosition);
	print TOTALA5SS join("\t", $catalogA5ss{$asPosition}, @tttt) . "\n";
}

@asPosition = ();
@asPosition = keys(%catalogA3ss);
@asPosition = sort@asPosition;
foreach $asPosition(@asPosition){
	@tttt = ();
	@tttt = split(/#/, $asPosition);
	print TOTALA3SS join("\t", $catalogA3ss{$asPosition}, @tttt) . "\n";
}

@asPosition = ();
@asPosition = keys(%catalogSe);
@asPosition = sort@asPosition;
foreach $asPosition(@asPosition){
	@tttt = ();
	@tttt = split(/#/, $asPosition);
	print TOTALSE join("\t", $catalogSe{$asPosition}, @tttt) . "\n";
}

@asPosition = ();
@asPosition = keys(%catalogRi);
@asPosition = sort@asPosition;
foreach $asPosition(@asPosition){
	@tttt = ();
	@tttt = split(/#/, $asPosition);
	print TOTALRI join("\t", $catalogRi{$asPosition}, @tttt) . "\n";
}

@asPosition = ();
@asPosition = keys(%catalogMxe);
@asPosition = sort@asPosition;
foreach $asPosition(@asPosition){
	@tttt = ();
	@tttt = split(/#/, $asPosition);
	print TOTALMXE join("\t", $catalogMxe{$asPosition}, @tttt) . "\n";
}

close TOTALA5SS;
close TOTALA3SS;
close TOTALSE;
close TOTALRI;
close TOTALMXE;

#obtain the num of AS with at least one read support in IJC or SJC
sub getNumWithReadSupport{
	my $asFile = $_[0];
	my (@tt, $asNum);
	open ASFF, "<$asFile";	
	<ASFF>;
	while(my $line = <ASFF>){
		#ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen
		#0       4       21                      65      49
		@tt = ();
		@tt = split(/\t/, $line);
		if($tt[1]!=0 or $tt[2]!=0){
			$asNum++;
		}
	}
	close ASFF;
	return $asNum;
}

sub checkMinReadLen{
	my ($readLengthX, $filterMinReadLenX) = @_;
	my (@tt, $i);
	@tt = split(/,/, $readLengthX);
	foreach my $len(@tt){
		#print "enter judege read length\nreadLen:$len\tfilterMinReadLen:$filterMinReadLenX\n";
		return 0 if($len < $filterMinReadLenX);
	}
	#print "meet read length critera\n";
	return 1;
}

sub checkLibLayout{
	my ($layoutX, $filterLibLayoutX) = @_;
	my (@tt, $i);
	@tt = split(/,/, $layoutX);
	foreach my $tmpLayout(@tt){
		#print "enter judege layout\nlayout:$tmpLayout\tfilterLibLayout:$filterLibLayoutX\n";
		return 0 if(index($filterLibLayoutX, $tmpLayout)<0);
	}
	#print "meet layout critera\n";
	return 1;
}
sub checkLibType{
	my ($libraryX, $filterLibTypeX) = @_;
	#print "Enter check lib type\n";
	#print "LibraryType:" . $libraryX . "xxx";
	#<STDIN>;
	#print "FilterLibType:$filterLibTypeX" . "xxx\n";
	#<STDIN>;
	my (@tt, $i);
	@tt = ();
	@tt = split(/,/, $libraryX);
	foreach my $libType(@tt){
		#print "enter judge library type.\nlibraryType:$libType\tfilterLibType:$filterLibTypeX\n";
		return 0 if(index($filterLibTypeX, $libType)<0);
	}

	#print "meet lib type critera\n";
	return 1;
}
