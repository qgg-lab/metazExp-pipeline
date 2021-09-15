#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV < 0){
	print "This script is used to gather psi of AS in experiment and experiment sequencing information\n\n";
	print "\t perl $0 \\\n" . 
		"\t\t --psiOutputDir    ../008-pickup-psi-of-ASs-in-all-expers/psiOutputDir \\\n" . 
		"\t\t --outputExperimentInfo experimentInfo.tsv \n";
	exit(0);
}

my ($psiOutputDir, $outputExperimentInfo);

GetOptions(
        'psiOutputDir=s'=>\$psiOutputDir,
        'outputExperimentInfo=s'=>\$outputExperimentInfo,
);

# detect psiOutputDir
if(not -e $psiOutputDir){
	print STDERR "$psiOutputDir doesn't exist!\n";
	exit;
}

my ($expDirText, @expDir);
# process psiOutputDir to ensure this dir with "/" in the end.
$expDirText = "";
$psiOutputDir = $psiOutputDir . "/" if(substr($psiOutputDir, length($psiOutputDir)-1, 1) ne "/");

$expDirText = `find $psiOutputDir -type d -maxdepth 1`;
 @expDir = split(/\n/, $expDirText);
shift(@expDir);

# gather information from experiment dir
# 从每个 experiment中收集 RNAseq建库信息(library type, layout, read lenth等)、比对信息(比对率)、以及检测的AS数量等信息
my ($experimentDir);
my ($status, $runNum, $runId, $library, $layout, $phredScore);
my ($spotNum, $readLength, $alignPer);
my ($novelSpliceNum);
my ($annoCovGeneNum, $annoCovTranscriptNum);
my (%ASNUM, %novelASNUM, %JCECASNUMwithreads, %JCASNUMwithreads);
my (@tmpArr, @tmpAsField, $expId);

#############################
# 收集信息产生16个文件
# (1) experimentInfo : A)测序和比对信息; B)实验中识别的各类总AS数量; C)novel AS的数量; D)jcec和jc登记的各类AS的数量(与AS总数是相等的)
# (2) 5个JCEC: ASID已经重新统一编过号，编号来自于(4)
# (3) 5个JC: ASID已经重新编过号，编号来自于(4)
# (4) 5个总的AS位置信息:
# ###########################
open EXPERIMENT, ">$outputExperimentInfo";
print EXPERIMENT join("\t", "expId", "status", 
		"runNum", "runId", "library", "layout", "phredScore", "readLength", "spotNum(M)", "alignPer(%)", "novelSpliceNum",
                "totalA5SS", "totalA3SS", "totalSE", "totalRI", "totalMXE",
		"novelA5SS", "novelA3SS", "novelSE", "novelRI", "novelMXE",
                "jcecA5SS", "jcecA3SS", "jcecSE", "jcecRI", "jcecMXE",
                "jcA5SS", "jcA3SS", "jcSE", "jcRI", "jcMXE") . "\n";

# scan each experiment to gather psi and experiment information
for(my $i=0; $i<=$#expDir; $i++){
	# 初始化每个experiment的信息
	$status = "OK";
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

	# 每个experiment中登记的AS数量	
	%ASNUM = ();
	%novelASNUM = ();
	%JCECASNUMwithreads = ();
	%JCASNUMwithreads =();

	@tmpArr = ();
	@tmpArr = split(/\//, $expDir[$i]);
	$expId = $tmpArr[$#tmpArr];

	$experimentDir = $expDir[$i] . "/";

	# detect running status
	if(-e $experimentDir . "/" . "fromGTF.A5SS.txt"){
		$status = "OK";
	}else{
		$status = "ERROR";
	}

	# 获得新检测的splice的数量（在gtf注释之外的splice）
        $novelSpliceNum = "";
        if(-e $experimentDir . "/" . $expId . ".novel.splicesite.from.hisat.txt"){
                my $cmd = "wc -l " . $experimentDir . "/" . $expId . ".novel.splicesite.from.hisat.txt";
                $novelSpliceNum = `$cmd`;
                chomp($novelSpliceNum);
                if($novelSpliceNum=~/(\d+) .*/){
                        $novelSpliceNum = $1;
                }
        }

	# 获得experiment中的run ID，数量，以及整体experiment的比对信息
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

		#run数量大于1时，去掉最后的逗号,
		if($runNum > 0){
			$runId = substr($runId, 0, length($runId) -1);
			$layout =substr($layout, 0, length($layout) -1);
			$readLength = substr($readLength, 0, length($readLength) -1);
			$library = substr($library, 0, length($library) - 1);
			$phredScore = substr($phredScore, 0, length($phredScore) -1 );
		}
	}

	if(-e $experimentDir . "/" . "fromGTF.A5SS.txt"){

		############################################################################
		#									   #
		# 1. 获得该实验中报告的5种类型AS的数量，其中包括利用gtf注释获得但是RNAseq  #
		#    mapping未检测的结果    						   #
		#									   #
		############################################################################

		my $cmd = "wc -l " . $experimentDir . "/" . "fromGTF.A5SS.txt";
		$ASNUM{"A5SS"} = `$cmd`;
		chomp($ASNUM{"A5SS"});
		if($ASNUM{"A5SS"}=~/(\d+) .*/){
			$ASNUM{"A5SS"} = $1 - 1;
		}
	
		$cmd = "wc -l " . $experimentDir . "/" . "fromGTF.A3SS.txt";
		$ASNUM{"A3SS"} = `$cmd`;
		chomp($ASNUM{"A3SS"});
		if($ASNUM{"A3SS"}=~/(\d+) .*/){
			$ASNUM{"A3SS"} = $1 - 1;
		}
	
		$cmd = "wc -l " . $experimentDir . "/" . "fromGTF.SE.txt";
		$ASNUM{"SE"} = `$cmd`;
		chomp($ASNUM{"SE"});
		if($ASNUM{"SE"}=~/(\d+) .*/){
			$ASNUM{"SE"} = $1 - 1;
		}

		$cmd = "wc -l " . $experimentDir . "/" . "fromGTF.RI.txt";
		$ASNUM{"RI"} = `$cmd`;
		chomp($ASNUM{"RI"});
		if($ASNUM{"RI"}=~/(\d+) .*/){
			$ASNUM{"RI"} = $1 -1;
		}

		$cmd = "wc -l " . $experimentDir . "/" . "fromGTF.MXE.txt";
		$ASNUM{"MXE"} = `$cmd`;
		chomp($ASNUM{"MXE"});
		if($ASNUM{"MXE"}=~/(\d+) .*/){
			$ASNUM{"MXE"} = $1 - 1;
		}
	}

	# 获得novel AS的数量
	if(-e $experimentDir . "/" . "fromGTF.novelEvents.A5SS.txt"){

		$novelASNUM{"A5SS"} = `cat $experimentDir/fromGTF.novelEvents.A5SS.txt | wc -l`;
		$novelASNUM{"A5SS"} = $novelASNUM{"A5SS"} - 1;

		$novelASNUM{"A3SS"} = `cat $experimentDir/fromGTF.novelEvents.A3SS.txt | wc -l`;
		$novelASNUM{"A3SS"} = $novelASNUM{"A3SS"} - 1;

		$novelASNUM{"SE"} = `cat $experimentDir/fromGTF.novelEvents.SE.txt | wc -l`;
		$novelASNUM{"SE"} = $novelASNUM{"SE"} - 1;

		$novelASNUM{"RI"} = `cat $experimentDir/fromGTF.novelEvents.RI.txt | wc -l`;
		$novelASNUM{"RI"} = $novelASNUM{"RI"} - 1;

		$novelASNUM{"MXE"} = `cat $experimentDir/fromGTF.novelEvents.MXE.txt | wc -l`;
		$novelASNUM{"MXE"} = $novelASNUM{"MXE"} - 1;
	
	}

	############################################################################
	#									   #
	# 3. 获得各类AS的inclusion和exclusion 有mapped reads支持的AS数量           #
	#									   #
	############################################################################

        #extract num of AS with supported reads in JCEC
	if(-e $experimentDir . "/" . "JCEC.raw.input.A5SS.txt"){
		$JCECASNUMwithreads{"A5SS"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.A5SS.txt");
		$JCECASNUMwithreads{"A3SS"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.A3SS.txt");
		$JCECASNUMwithreads{"SE"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.SE.txt");
		$JCECASNUMwithreads{"RI"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.RI.txt");
		$JCECASNUMwithreads{"MXE"} = &getNumWithReadSupport($experimentDir . "/" . "JCEC.raw.input.MXE.txt");
	}
        
	#extract num of AS with supported reads in JC
	if(-e $experimentDir . "/" . "JC.raw.input.A5SS.txt"){
		$JCASNUMwithreads{"A5SS"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.A5SS.txt");
		$JCASNUMwithreads{"A3SS"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.A3SS.txt");
		$JCASNUMwithreads{"SE"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.SE.txt");
		$JCASNUMwithreads{"RI"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.RI.txt");
		$JCASNUMwithreads{"MXE"} = &getNumWithReadSupport($experimentDir . "/" . "JC.raw.input.MXE.txt");
	}

	
	#output experiment information
	print EXPERIMENT join("\t", $expId, $status, $runNum, $runId, $library, $layout, $phredScore, $readLength, 
		$spotNum, $alignPer, $novelSpliceNum,
		$ASNUM{"A5SS"}, $ASNUM{"A3SS"}, $ASNUM{"SE"}, $ASNUM{"RI"}, $ASNUM{"MXE"},
		$novelASNUM{"A5SS"}, $novelASNUM{"A3SS"}, $novelASNUM{"SE"}, $novelASNUM{"RI"}, $novelASNUM{"MXE"},
		$JCECASNUMwithreads{"A5SS"}, $JCECASNUMwithreads{"A3SS"}, $JCECASNUMwithreads{"SE"}, $JCECASNUMwithreads{"RI"}, $JCECASNUMwithreads{"MXE"},
		$JCASNUMwithreads{"A5SS"}, $JCASNUMwithreads{"A3SS"}, $JCASNUMwithreads{"SE"}, $JCASNUMwithreads{"RI"}, $JCASNUMwithreads{"MXE"}) . "\n";
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
