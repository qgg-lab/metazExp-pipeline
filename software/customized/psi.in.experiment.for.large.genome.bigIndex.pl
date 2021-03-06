#!/usr/bin/perl
use Getopt::Long;
use List::Util qw/max min/;
use strict;
if($#ARGV<0){
	print "$0 \n" . 
		" \t\t--expId ERX2073443\n" . 
		" \t\t--runIds ERR2013772,ERR2013771,ERR2013774,ERR2013773 \\\n" . 
		" \t\t--genomeDb /mnt/home/liujind1/workAS/01-cattle/database/genome \\\n" .
		" \t\t--cDNAdb /mnt/home/liujind1/workAS/01-cattle/database/cDNA \\\n" .
		" \t\t--gtfFile /mnt/home/liujind1/workAS/01-cattle/database/annotation.gtf \\\n" .
		" \t\t--rMats /mnt/home/liujind1/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS2/rmats.py \\\n" .
		" \t\t--fastqOutputDir /mnt/scratch/liujind1/fastqOutputDir \\\n" .
		" \t\t--psiOutputDir /mnt/scratch/liujind1/psiOutputDir \\\n" .
		" \t\t--logOutputDir /mnt/scratch/liujind1/logOutputDir \\\n" .
		" \t\t--manualDownloadDir /mnt/home/liujind1/workAS/01-cattle/003-reAsmbl-failed-goodExps/ \\\n" .
		" \t\t--maxIntronSize 12000 \\\n" .
		" \t\t--minIntronSize 10 \\\n" .
		" \t\t--switchPsi    on \\\n"  .
		" \t\t--switchAssemble off \\\n" . 
		" \t\t--threadNum 18 \n";
	exit(0);
}

my ($expId, $runIdList, $fastqOutputDir, $cDNAdb, $genomeDb, $psiOutputDir, $logOutputDir, $manualDownloadDir);
my ($gtfFile, $rmatsPy, $threadNum, $cmd, $switchPsi, $switchAssemble, $maxIntronSize, $minIntronSize);

GetOptions(
        'expId=s'=>\$expId,
        'runIds=s'=>\$runIdList,
	'genomeDb=s'=>\$genomeDb,
	'cDNAdb=s'=>\$cDNAdb,
	'gtfFile=s'=>\$gtfFile,
	'fastqOutputDir=s'=>\$fastqOutputDir,
	'psiOutputDir=s'=>\$psiOutputDir,
	'logOutputDir=s'=>\$logOutputDir,
	'manualDownloadDir=s'=>\$manualDownloadDir,
	'rmatsPy=s'=>\$rmatsPy,
	'maxIntronSize=s'=>\$maxIntronSize,
	'minIntronSize=s'=>\$minIntronSize,
	'switchPsi=s'=>\$switchPsi,
	'switchAssemble=s'=>\$switchAssemble,
	'threadNum=i'=>\$threadNum,
);

#prepare workpath
$fastqOutputDir = substr($fastqOutputDir, 0, length($fastqOutputDir)-1) if(substr($fastqOutputDir, length($fastqOutputDir)-1, 1) eq "/");
$psiOutputDir = substr($psiOutputDir, 0, length($psiOutputDir)-1) if(substr($psiOutputDir, length($psiOutputDir)-1, 1) eq "/");
$logOutputDir = substr($logOutputDir, 0, length($logOutputDir)-1) if(substr($logOutputDir, length($logOutputDir)-1, 1) eq "/");

$fastqOutputDir = $fastqOutputDir . "/" . $expId if(index($fastqOutputDir, $expId)<0);
$cmd = "rm -rf $fastqOutputDir";
print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
system($cmd);

$cmd = "mkdir -p $fastqOutputDir";
print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
system($cmd);

$psiOutputDir = $psiOutputDir . "/" . $expId if(index($psiOutputDir, $expId)<0);
$cmd = "rm -rf $psiOutputDir";
print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
system($cmd);

$cmd = "mkdir -p $psiOutputDir";
print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
system($cmd);

$logOutputDir = $logOutputDir . "/" . $expId if(index($logOutputDir, $expId)<0);
$cmd = "rm -rf $logOutputDir";
print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
system($cmd);

$cmd = "mkdir -p $logOutputDir";
print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
system($cmd);

my (@phredScore, @phredScoreEvidence, @layout, @libraryType, @libraryTypeEvidence, @readLen);
my ($libraryTypeEvidence, $phredScoreEvidence, $phredScoreEvidence);
my (%checkPhredScoreUniq, %checkLibraryTypeUniq, %checkLayoutUniq);
my ($dumpRunDataFailedNum, %notDumpRunId);

my $expSeqInfo = $psiOutputDir . "/" . $expId . ".SeqInfo.txt";
open EXPERIMENTINFO, ">$expSeqInfo";
print EXPERIMENTINFO join("\t", "ExpId", "RunId", "Layout", "ReadLen", "PhredScore", "PhredDetail", "Library", "LibraryDetail") . "\n";

$dumpRunDataFailedNum = 0;
my @runId = split(/,/, $runIdList);
for(my $i=0; $i<=$#runId; $i++){
	# ??????fastq dump ???????????? fastq??????
	$cmd = "fasterq-dump " .
		" " .  $runId[$i] .
		" --split-files " .
		" -O $fastqOutputDir " .
		" --threads $threadNum " .
		" -t " . $fastqOutputDir . "/tmp." . $runId[$i] .
		" 1> $logOutputDir/log.o.online-dump.$runId[$i] " .
		" 2> $logOutputDir/log.e.online-dump.$runId[$i]";
	print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";

	system($cmd);

	# ???????????????fastq-dump???????????????
	if(-e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq" or -e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq" or -e $fastqOutputDir . "/" . $runId[$i] . ".fastq"){

		print "[" . &getTime() . "]: Online dump " . $runId[$i] . " successfully.\n\n";

		if(-e $fastqOutputDir . "/" . $runId[$i] . ".fastq"){
			# single end: only generate SRR120999.fastq
			my $cmd = "mv " . $fastqOutputDir . "/" . $runId[$i] . ".fastq " . $fastqOutputDir . "/" . $runId[$i] . "_1.fastq";
			print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
			system($cmd);

		}else{
			# paired end: generate SRR12099999_1.fastq, SRR120999999_2.fastq
		}

		#remove temporary directory for dump
		my $cmd = "rm -rf " . $fastqOutputDir . "/tmp." . $runId[$i];
		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
		system($cmd);	

	}else{
	# ???????????????fastq-dump??????????????????
		
		# ????????????????????????????????????????????????????????????????????????????????????????????????sra?????????$fastqOutputDir . "/" . $runId[$i] . ".sra"
		print "[" . &getTime() . "] Can't online dump " . $runId[$i] . ". Try download sra file from sradb.\n\n";

		# ????????????????????????????????????
		if(&getFileSize($manualDownloadDir . "/" . $runId[$i] . ".sra")==1){
			my $cmd = "cp " . $manualDownloadDir . "/" . $runId[$i] . ".sra " . $fastqOutputDir . "/" . $runId[$i] . ".sra";
			print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
			system($cmd);
			goto DUMPSRA;
		}

		# ?????????????????????ftp???????????????
		my $sraAddress = "ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" . 
			substr($runId[$i], 0, 3) . "/" . 
			substr($runId[$i], 0, 6) . "/" . 
			$runId[$i] . "/" . 
			$runId[$i] . ".sra";
		my $cmd = "wget -O " .	$fastqOutputDir . "/" . $runId[$i] . ".sra" . 
			" --tries=10 " . 
			" -o $logOutputDir/log.o.wget." . $runId[$i] .
			" $sraAddress " .
			" 2> " . $logOutputDir . "/log.e.wget." . $runId[$i] . ".sra";

		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
		system($cmd);

		# ??????????????????sra????????????????????????????????????DUMPSRA??????sra????????????fastq??????
		goto DUMPSRA if(&getFileSize($fastqOutputDir . "/" . $runId[$i] . ".sra")==1);


		my ($srapubrunId, $sosSraSuffix);
		############################################
		# (1)???????????????https?????????: sra-pub-run-5,4,3,2,1??????????????????5,4,3,2,1???srr??????
		for($srapubrunId=12; $srapubrunId>=1; $srapubrunId--){
			for($sosSraSuffix=9; $sosSraSuffix>=1; $sosSraSuffix--){
				my $sraAddress = "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos/sra-pub-run-" . $srapubrunId . "/" . $runId[$i] . "/" . $runId[$i] . ".$sosSraSuffix";
				my $cmd = "wget -O " .	$fastqOutputDir . "/" . $runId[$i] . ".sra" . 
				" --tries=3 " . 
				" -o $logOutputDir/log.o.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] .
				" $sraAddress " .
				" 2> " . $logOutputDir . "/log.e.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] . ".sra";

				print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
				system($cmd);
				system("rm -rf " . "$logOutputDir/log.o.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i]);
				system("rm -rf " . $logOutputDir . "/log.e.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] . ".sra");
				# ??????????????????sra????????????????????????????????????DUMPSRA??????sra????????????fastq??????
				goto DUMPSRA if(&getFileSize($fastqOutputDir . "/" . $runId[$i] . ".sra")==1);
			}

			for($sosSraSuffix=9; $sosSraSuffix>=1; $sosSraSuffix--){
				my $sraAddress = "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-" . $srapubrunId . "/" . $runId[$i] . "/" . $runId[$i] . ".$sosSraSuffix";
				my $cmd = "wget -O " .	$fastqOutputDir . "/" . $runId[$i] . ".sra" . 
				" --tries=3 " . 
				" -o $logOutputDir/log.o.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] .
				" $sraAddress " .
				" 2> " . $logOutputDir . "/log.e.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] . ".sra";

				print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
				system($cmd);
				system("rm -rf " . "$logOutputDir/log.o.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i]);
				system("rm -rf " . $logOutputDir . "/log.e.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] . ".sra");
	
				# ??????????????????sra????????????????????????????????????DUMPSRA??????sra????????????fastq??????
				goto DUMPSRA if(&getFileSize($fastqOutputDir . "/" . $runId[$i] . ".sra")==1);
			}

			for($sosSraSuffix=9; $sosSraSuffix>=1; $sosSraSuffix--){
				my $sraAddress = "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-" . $srapubrunId . "/" . $runId[$i] . "/" . $runId[$i] . ".$sosSraSuffix";
				my $cmd = "wget -O " .	$fastqOutputDir . "/" . $runId[$i] . ".sra" . 
				" --tries=3 " . 
				" -o $logOutputDir/log.o.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] .
				" $sraAddress " .
				" 2> " . $logOutputDir . "/log.e.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] . ".sra";

				print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
				system($cmd);
				system("rm -rf " . "$logOutputDir/log.o.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i]);
				system("rm -rf " . $logOutputDir . "/log.e.wget.from.sra-pub-run-" . $srapubrunId . ".$sosSraSuffix." . $runId[$i] . ".sra");
	
				# ??????????????????sra????????????????????????????????????DUMPSRA??????sra????????????fastq??????
				goto DUMPSRA if(&getFileSize($fastqOutputDir . "/" . $runId[$i] . ".sra")==1);
			}

			# ?????????????????????1???RUN
			my $sraAddress = "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos/sra-pub-run-" . $srapubrunId . "/" . $runId[$i] . "/" . $runId[$i];
			my $cmd = "wget -O " .	$fastqOutputDir . "/" . $runId[$i] . ".sra" . 
			" --tries=10 " . 
			" -o $logOutputDir/log.o.wget.from.sra-pub-run-" . $srapubrunId . "." . $runId[$i] .
			" $sraAddress " .
			" 2> " . $logOutputDir . "/log.e.wget.from.sra-pub-run-" . $srapubrunId . "." . $runId[$i] . ".sra";

			print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
			system($cmd);
			system("rm -rf " . "$logOutputDir/log.o.wget.from.sra-pub-run-" . $srapubrunId . "." . $runId[$i]);
			system("r m-rf " . $logOutputDir . "/log.e.wget.from.sra-pub-run-" . $srapubrunId . "." . $runId[$i] . ".sra");
	
			# ??????????????????sra????????????????????????????????????DUMPSRA??????sra????????????fastq??????
			goto DUMPSRA if(&getFileSize($fastqOutputDir . "/" . $runId[$i] . ".sra")==1);	
		}


DUMPSRA:
		if(not -s $fastqOutputDir . "/" . $runId[$i] . ".sra"){
		# ??????sra????????????????????????		
			print STDERR "[" . &getTime() . "]" . " Not only can't online dump " . $runId[$i] . ", but also can't download sra file from sradb.\n\n";
			$dumpRunDataFailedNum = $dumpRunDataFailedNum + 1;
			$notDumpRunId{$runId[$i]} = 1;
			next;

		}else{
		# ??????sra??????????????????????????????????????????dump ??????fastq??????
			$cmd = "fasterq-dump " .
				" " .  $fastqOutputDir . "/" . $runId[$i] . ".sra" .
				" --split-files " .
				" -O $fastqOutputDir " .
				" --threads $threadNum " .
				" -t " . $fastqOutputDir . "/tmp." . $runId[$i] .
				" 1> $logOutputDir/log.o.online-dump.$runId[$i] " .
				" 2> $logOutputDir/log.e.online-dump.$runId[$i]";
			print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
			system($cmd);
			system("rm -rf " . "$logOutputDir/log.o.online-dump.$runId[$i]");
			system("rm -rf " . "$logOutputDir/log.e.online-dump.$runId[$i]");
		
			# Paired fastq file format: SRR1209999.sra_1.fastq and SRR1209999.sra_2.fastq
			if(-e $fastqOutputDir . "/" . $runId[$i] . ".sra_1.fastq"){
				
				my $cmd = "mv " . $fastqOutputDir . "/" . $runId[$i] . ".sra_1.fastq " . $fastqOutputDir . "/" . $runId[$i] . "_1.fastq";
				print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
				system($cmd);

			}

			if(-e $fastqOutputDir . "/" . $runId[$i] . ".sra_2.fastq"){
				my $cmd = "mv " . $fastqOutputDir . "/" . $runId[$i] . ".sra_2.fastq " . $fastqOutputDir . "/" . $runId[$i] . "_2.fastq";
				print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
				system($cmd);
			}

			#Single-End fastq file format: SRR1209999.sra.fastq
			if(-e $fastqOutputDir . "/" . $runId[$i] . ".sra.fastq"){
				my $cmd = "mv " . $fastqOutputDir . "/" . $runId[$i] . ".sra.fastq " . $fastqOutputDir . "/" . $runId[$i] . "_1.fastq";
				print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
				system($cmd);
			}

			#remove temporary directory for dump
			my $cmd = "rm -rf " . $fastqOutputDir . "/tmp." . $runId[$i];
			print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
			system($cmd);
			
			#remove sra
			my $cmd = "rm -rf " . $fastqOutputDir . "/" . $runId[$i] . ".sra" ;
			print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
			system($cmd);
		}
	}
	# above: dump fastq data over.

	##############################
	# obtain phred score         #
	##############################	
	# obtain phred score in any one run data
	if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq"){
		$phredScore[$i] = &calPhredScore($fastqOutputDir . "/" . $runId[$i] . "_1.fastq", 10000, \$phredScoreEvidence);
	}else{
	# some paired sra only dump SRR1209999_2.fastq, so check _2.fastq
		$phredScore[$i] = &calPhredScore($fastqOutputDir . "/" . $runId[$i] . "_2.fastq", 10000, \$phredScoreEvidence);
	}


	#####################
	# obtain layout     #
	#####################
	$layout[$i] = "SINGLE";
	if(-e("$fastqOutputDir/$runId[$i]" . "_1.fastq") and -e ("$fastqOutputDir/$runId[$i]" . "_2.fastq")){
		$layout[$i] = "PAIRED";
	}

	###############################################
	# begin to check library type: RF, FR, F, R   #
	###############################################
	#generate tmp fastq file
	if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq"){
		my $cmd = "head -n 4000000 $fastqOutputDir" . "/" . $runId[$i] . "_1.fastq " .
		"> $fastqOutputDir" . "/" . $runId[$i] . "_1.fq.tmp " .
		"2> " . $logOutputDir .  "/log.e.head100000Read.from.1.fq." . $runId[$i];
		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
		system($cmd);
	}

	if( -e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq"){
		my $cmd = "head -n 4000000 $fastqOutputDir" . "/" . $runId[$i] . "_2.fastq " .
		"> $fastqOutputDir" . "/" . $runId[$i] . "_2.fq.tmp " .
		"2> " . $logOutputDir .  "/log.e.head100000Read.from.2.fq." . $runId[$i];
		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
		system($cmd) if( -e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq");
	}

	# for paired layout
	if($layout[$i] eq "PAIRED"){

		my $cmd = "hisat2 -x $cDNAdb" .
			  " -1 " . $fastqOutputDir . "/" . $runId[$i] . "_1.fq.tmp" . 
			  " -2 " . $fastqOutputDir . "/" . $runId[$i] . "_2.fq.tmp" . 
			  " --phred" . $phredScore[$i] . 
			  " --no-unal" .
			  " -p " . $threadNum . 
			  " -S $fastqOutputDir" . "/" . "$runId[$i].map.cDNA.sam " .  
			  " 2> " . $logOutputDir .  "/log.e.hisat2.map.cDNA." . $runId[$i];
		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
		system($cmd);

	}else{
	# for single end layout
		
		# for regular single end fastq data
		my $cmd = "hisat2 -x $cDNAdb " . 
			  " -U " . $fastqOutputDir . "/" . $runId[$i] . "_1.fq.tmp" . 
			  " --phred" . $phredScore[$i] . 
			  " -p " . $threadNum . 
			  " --no-unal" .
			  " -S $fastqOutputDir" . "/" . "$runId[$i].map.cDNA.sam" .  
			  " 2> " . $logOutputDir .  "/log.e.hisat2.map.cDNA." . $runId[$i];
		if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq"){
			print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
			system($cmd);
		}

		# for irregular single end fastq which derive from paired sra with only _2.fastq data
		my $cmd = "hisat2 -x $cDNAdb " . 
			  " -U " . $fastqOutputDir . "/" . $runId[$i] . "_2.fq.tmp" . 
			  " --phred" . $phredScore[$i] . 
			  " -p " . $threadNum . 
			  " --no-unal" .
			  " -S $fastqOutputDir" . "/" . "$runId[$i].map.cDNA.sam" .  
			  " 2> " . $logOutputDir .  "/log.e.hisat2.map.cDNA." . $runId[$i];
		if( -e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq"){
			print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
			system($cmd);
		}
	}

	# obtain percentage of various type library
	$libraryType[$i] = &obtainLibraryType(	$fastqOutputDir. "/" . $runId[$i] . ".map.cDNA.sam", 
						1000000,
						$layout[$i], 
						95, 
						\$libraryTypeEvidence);


	#######################
	# get read length     #
	#######################
	if( -e $fastqOutputDir. "/" . $runId[$i] . "_1.fastq"){
		$readLen[$i] = &obtainReadLen(	$fastqOutputDir. "/" . $runId[$i] . "_1.fastq", 1000000);
	}else{
		$readLen[$i] = &obtainReadLen(  $fastqOutputDir. "/" . $runId[$i] . "_2.fastq", 1000000);
	}


	################################
	# output run data information  #
	################################
	
	# register run information into file
	print EXPERIMENTINFO join("\t", 
				$expId, $runId[$i], 
				$layout[$i],
				$readLen[$i], 
				$phredScore[$i], 
				$phredScoreEvidence,
				$libraryType[$i], 
				$libraryTypeEvidence) . "\n";

	# output run information into command file
	print join("\t", "ExpId", "RunId", "Layout", "ReadLen", "PhredScore", "PhredDetail", "Library", "LibraryDetail") . "\n";
	print join("\t", 
				$expId, $runId[$i], 
				$layout[$i],
				$readLen[$i], 
				$phredScore[$i], 
				$phredScoreEvidence,
				$libraryType[$i], 
				$libraryTypeEvidence) . "\n\n";


	# delete temporary fastq data file
	if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fq.tmp"){
		my $cmd = "rm -rf " . "$fastqOutputDir" . "/" . $runId[$i] . "_1.fq.tmp";
		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
		system($cmd);
	}

	if(-e "$fastqOutputDir" . "/" . $runId[$i] . "_2.fq.tmp"){
		$cmd = "rm -rf " . "$fastqOutputDir" . "/" . $runId[$i] . "_2.fq.tmp";
		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
		system($cmd);
	}

	# delete temporary alignment file against to cDNA database
	my $cmd = "rm -rf " . $fastqOutputDir. "/" . $runId[$i] . ".map.cDNA.sam";
	print "CMD[" . &getTime() . "]: " . $cmd . "\n\n"; 
	system($cmd);
	
	# register phred score, library type and layout into hash
	$checkPhredScoreUniq{$phredScore[$i]} = 1;
	$checkLibraryTypeUniq{$libraryType[$i]} = 1 if($libraryType[$i] ne "NA");
	# if library type can't be determined. its library type will be considered to be same with other runs
	$checkLayoutUniq{$layout[$i]} = 1;
}

############################################
# output the num of runs without dump data #
############################################
if($dumpRunDataFailedNum > 0){
	print STDERR "Can't dump " . $dumpRunDataFailedNum . " Run data from sradb not only online but also local.\n\n";
}

############################################
# check whether all runs have same layout, #
# library type and phred                   #
############################################
my (@keyPhredScore, @keyLibraryType, @keyLayout, @notDumpRunId);
@keyPhredScore = keys(%checkPhredScoreUniq);
@keyLibraryType = keys(%checkLibraryTypeUniq);
@keyLayout = keys(%checkLayoutUniq);
@notDumpRunId = keys(%notDumpRunId);

if($#keyPhredScore != 0 ){
	print STDERR "PhredScores of Runs in Experiment are not consistent.\n\n";
}

if($#keyLayout != 0){
	print STDERR "Layout of Runs in Experiment are not consistent.\n\n";
}

if($#keyLibraryType > 1){

	# $#keyLayout eq -1 suggesting the library type of all Runs can't be determined.
	print STDERR "Library type of Runs in Experiment are not consistent.\n\n";

}elsif($#keyLibraryType ==  1){

	print STDERR "Library type of Runs in Experiment are not consistent,\n";
	print STDERR "but library type will be all changed to U or UN.\n\n";

	if($keyLayout[0] eq "PAIRED"){

		$keyLibraryType[0] = "UN";

	}elsif($keyLayout[0] eq "SINGLE"){

		$keyLibraryType[0] = "U";

	}
	
	# ?????????????????????????????????libraryType????????????????????????????????????
	close EXPERIMENTINFO;
	my $tmpLibraryType = $keyLibraryType[0];
	system("sed -i \'s/^\\(.*\\t.*\\t.*\\t.*\\t.*\\t.*\\t\\)\\(.*\\)\\(\t.*\\)\$/\\1$tmpLibraryType\\3/\' $expSeqInfo");
	open EXPERIMENTINFO, ">>$expSeqInfo";

}

if($#keyPhredScore != 0 or $#keyLibraryType > 1 or $#keyLibraryType < 0 or $#keyLayout != 0 or $#runId == $#notDumpRunId){
	# libraryType ???2?????????????????????????????????????????????????????????????????????U??????UN???
	# libraryType ???0???????????????????????????????????????????????????????????????????????????????????????
	exit;
}



#######################################
# map fastq data of all runs to genome#
#######################################

# generate phredscore, fastqfiles and libraryType parameters
my ($fastqParameter, $fastqParameter1, $fastqParameter2, $phredScoreParameter, $hisat2LibraryTypeParameter);
# phred
$phredScoreParameter = "--phred" . $keyPhredScore[0];

# rna-strandness
$hisat2LibraryTypeParameter = "";
if($#keyLibraryType == 0){

	if($keyLibraryType[0] ne "UN" and $keyLibraryType[0] ne "U"){
		$hisat2LibraryTypeParameter = " --rna-strandness " . $keyLibraryType[0];
	}

}elsif($#keyLibraryType == -1){
	# $#keyLibraryType == -1 suggesting the library type of all Runs can't be determined.
	# So hisat2LibraryTypeParameter should be "" and unstrand
	$hisat2LibraryTypeParameter = "";

}

################################
# prepare fastqFile parameters #
################################
$fastqParameter1 = "";
$fastqParameter2 = "";
for(my $i=0; $i<=$#runId; $i++){

	# discard the runs without dumpped fastq data
	if(not exists($notDumpRunId{$runId[$i]})){
		$fastqParameter1 .= $fastqOutputDir . "/" . $runId[$i] . "_1.fastq," if(-e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq");
		$fastqParameter2 .= $fastqOutputDir . "/" . $runId[$i] . "_2.fastq," if(-e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq");
	}

}

# remove last comma
$fastqParameter1 = substr($fastqParameter1, 0, length($fastqParameter1)-1) if($fastqParameter1 ne "");
$fastqParameter2 = substr($fastqParameter2, 0, length($fastqParameter2)-1) if($fastqParameter2 ne "");

# assemble final fastqParameter with fastqParameter1 and fastqParameter2
$fastqParameter = "";
if($layout[0] eq "SINGLE"){
	
	# for irregular single end fastq files which come frm paired sra with only _2.fastq
	
	if($fastqParameter1 ne "" and $fastqParameter2 eq ""){

		$fastqParameter = $fastqParameter1;

	}elsif($fastqParameter1 ne "" and $fastqParameter2 ne ""){

		$fastqParameter = $fastqParameter1 . "," . $fastqParameter2;

	}elsif($fastqParameter1 eq "" and $fastqParameter2 ne ""){
		
		$fastqParameter = $fastqParameter2;

	}

	$fastqParameter = " -U " . $fastqParameter;
}else{

	$fastqParameter = " -1 " . $fastqParameter1 . " -2 " . $fastqParameter2;	

}

##################################################################
#                                                                #
# call hisat2 and samtools to generate sorted bam alignment file #
#                                                                #
##################################################################

my $cmd = "hisat2 -x $genomeDb " . 
	$fastqParameter . " " . 
	$phredScoreParameter . " " . 
	$hisat2LibraryTypeParameter . 
	" --max-intronlen " . $maxIntronSize .
	" --min-intronlen " . $minIntronSize .
	" -p " . $threadNum . 
	" --novel-splicesite-outfile " . $psiOutputDir  . "/" . $expId . ".novel.splicesite.from.hisat.txt " .
	" --summary-file " . $psiOutputDir  . "/" . $expId . ".alignSum.from.hisat.txt " . 
	" 2>" . $logOutputDir .  "/log.e.map." . $expId . ".to.genome" . 
	" | samtools view -bS - " . 
	" 2> $logOutputDir/" . "log.e.samView." . $expId . 
	" | samtools sort - -o " . $fastqOutputDir . "/" . $expId . ".bam" . 
	" 2>" . $logOutputDir . "/log.e.samSort." . $expId;

print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
system($cmd);

# delete fastq file
for(my $i=0; $i<=$#runId; $i++){
	if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq"){
		$cmd = "rm -rf " . $fastqOutputDir . "/" . $runId[$i] . "_1.fastq";
		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
        	system($cmd);
	}

	if(-e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq"){
		$cmd = "rm -rf " . $fastqOutputDir . "/" . $runId[$i] . "_2.fastq";
		print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
		system($cmd);
	}
}


# output alignment percentage into experiment information file
my ($totalSpotNum, $alignmentPercentage);
&obtainAlignmentPercentage($psiOutputDir  . "/" . $expId . ".alignSum.from.hisat.txt", \$totalSpotNum, \$alignmentPercentage);
print EXPERIMENTINFO "TotalSpot:$totalSpotNum\tAlignPercent:$alignmentPercentage\n";
print "TotalSpot:$totalSpotNum\tAlignPercent:$alignmentPercentage\n\n";

# obtain library type for rMATs and stringTie
my $rMATsLibraryType = "";

# 20190527: if($#keyLibraryType == 0){
if($#keyLibraryType == 0 or $#keyLibraryType == 1){
	# specify library type for rMATS
	$rMATsLibraryType = $keyLibraryType[0];
	if($rMATsLibraryType eq "RF" or $rMATsLibraryType eq "R" ){

		$rMATsLibraryType = "fr-firststrand";

	}elsif($rMATsLibraryType eq "FR" or $rMATsLibraryType eq "F"){

		$rMATsLibraryType = "fr-secondstrand";

	}elsif($rMATsLibraryType eq "UN" or $rMATsLibraryType eq "U"){

		$rMATsLibraryType = "fr-unstranded";

	}

}elsif($#keyLibraryType == -1){
	# For experiment whose library type can't be determined,
	# its library type is set unstranded.
	$rMATsLibraryType = "fr-unstranded";

}

my $readLength = max @readLen;

#######################################
# assemble transcripts with stringTie #
# on the bam file                     #
####################################### 
my $stringtieLibrary = "";
if($rMATsLibraryType eq "fr-firststrand"){

	$stringtieLibrary = "--rf";

}elsif($rMATsLibraryType eq "fr-secondstrand"){

	$stringtieLibrary = "--fr";

}else{

	$stringtieLibrary = "";
}

if(uc($switchAssemble) eq "ON"){
	my $cmd = "stringtie " . $fastqOutputDir . "/" . $expId . ".bam" .
		  " -G " . $gtfFile . 
		  " -l " . $expId . 
		  " $stringtieLibrary " . 
		  " -o " . $psiOutputDir . "/" . "transcriptomeByStringtie.gtf" . 
		  " -C " . $psiOutputDir . "/" . "annotationCovByStringtie.gtf" .
		  " -A " . $psiOutputDir . "/" . "geneAbundanceByStringtie.tab" .
		  " -p " . $threadNum . 
		  " 2> " . $logOutputDir . "/log.e.stringtie.$expId";
	print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
	system($cmd);
}

#####################################
# generate b file to save bam file  #
# to run rMATS                      #
#####################################
open RMATSB, ">" . $fastqOutputDir . "/" . $expId . ".rMats.b.txt";
print RMATSB $fastqOutputDir . "/" . $expId . ".bam";
print "CMD[" . &getTime() . "]: " . "echo \"" . $fastqOutputDir . "/" . $expId . ".bam" . "\" > " . $fastqOutputDir . "/" . $expId . ".rMats.b.txt" . "\n\n";
close RMATSB;

# library parameter for rMATS
my $rMATsLayout = lc($keyLayout[0]);
if(uc($switchPsi) eq "ON"){
	my $cmd = "python " . 
	  " $rmatsPy " .
	  " --b1 " . $fastqOutputDir . "/" . $expId . ".rMats.b.txt" . 
	  " --gtf " . $gtfFile .
	  " --od " . $psiOutputDir .
	  " -t " . $rMATsLayout . 
	  " --nthread " . $threadNum . 
	  " --tstat " . $threadNum . 
	  " --libType " . $rMATsLibraryType . 
	  " --readLength " . $readLength .
	  " --statoff " . 
	  " 1> " . $logOutputDir . "/log.o.rMats.$expId" . 
	  " 2> " . $logOutputDir . "/log.e.rMats.$expId";

	print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
	system($cmd);
}

# ??????fastqOutputDir??????
$cmd = "rm -rf " . $fastqOutputDir;
print "CMD[" . &getTime() . "]: " . $cmd . "\n\n";
# ???????????????????????????????????????call.as.from.bam.pl?????????????????????sbatch???????????????
#system($cmd);

close EXPERIMENTINFO;

#Detect PhredScore
sub calPhredScore{
        my ($fastqFile, $readNum, $evidence) = @_;
        my ($ltZeroBaseCount, $read, $score, $totalBaseCount, $percentage);
        $ltZeroBaseCount = 0;
        open FPHRED, "<$fastqFile";
        for(my $i=0; $i<$readNum*4; $i++){
		last if(eof(FPHRED));
                $read = <FPHRED>;
                next if($i-int($i/4)*4 != 3);
                for(my $j=0; $j<=length($read)-2; $j++){
                        $score = ord(substr($read, $j, 1))-64;
                        $totalBaseCount++;
                        $ltZeroBaseCount++ if($score < 0);
                }
        }
	$percentage = int($ltZeroBaseCount/$totalBaseCount * 10000)/100;
	$$evidence = "NegativeScoreBaseCount(%):" . $ltZeroBaseCount . "(" . $percentage . ")";
        if($ltZeroBaseCount > 0){
                return 33;
        }else{
                return 64;
        }

	close FPHRED;
}

# detect library type such as RF, FR, UN, R, F, U
sub obtainLibraryType{
	my ($map_sam, $num, $layout, $supportPercentage, $evidence) = @_;
	my (%sam, $line, @tmp, @binarr, $binstr, $i, $line);
	my ($total_statistic_num, $RF_num, $FR_num, $R_num, $F_num, $percentage);

	goto SINGLE if($layout eq "SINGLE");

	open FSAM, "<$map_sam";
	$FR_num = 0;
	$RF_num = 0;
	while($line=<FSAM>){
		chomp($line);

		@tmp = ();
		@tmp = split(/\t/, $line);

		#It is header
		next if($#tmp<10);
		#parse alignment information
		@binarr = (0);
		$binstr=sprintf("%b", $tmp[1]);
        	for($i=length($binstr); $i>=1; $i--){
			$binarr[$#binarr+1]=substr($binstr, $i-1, 1);
	        }

		#not proper pair end reads
		next if(not($tmp[1]==99 or $tmp[1]==147 or $tmp[1]==83 or $tmp[1]==163));

		#two reads are not in a same refseq
		next if($tmp[6] ne "="); 

		#hit multiple sites
		next if(not($tmp[$#tmp]=~/NH:i:1/));
	
		#register read1 strand and position
		if($binarr[7]==1){
			${${$sam{$tmp[0]}}{"read1"}}{"position"}=$tmp[3];
			if($binarr[5]==1){
				${${$sam{$tmp[0]}}{"read1"}}{"strand"}="-";
			}else{
				${${$sam{$tmp[0]}}{"read1"}}{"strand"}="+";
			}
		}
	
		#register read2 strand and position
		if($binarr[8] == 1){
			${${$sam{$tmp[0]}}{"read2"}}{"position"}=$tmp[3];
			if($binarr[5]==1){
				${${$sam{$tmp[0]}}{"read2"}}{"strand"}="-";
			}else{
				${${$sam{$tmp[0]}}{"read2"}}{"strand"}="+";
			}
		}
		$num--;
		last if($num<=0);
	}
	close FSAM;

	#calculate numbers of RF and FR pair-ends
	@tmp = ();
	@tmp = keys%sam;

	for($i=0; $i<=$#tmp; $i++){
		if(${${$sam{$tmp[$i]}}{"read1"}}{"strand"} eq "+" and ${${$sam{$tmp[$i]}}{"read2"}}{"strand"} eq "-" and ${${$sam{$tmp[$i]}}{"read1"}}{"position"} <= ${${$sam{$tmp[$i]}}{"read2"}}{"position"}){
		$FR_num++;
		}elsif(${${$sam{$tmp[$i]}}{"read1"}}{"strand"} eq "-" and ${${$sam{$tmp[$i]}}{"read2"}}{"strand"} eq "+" and ${${$sam{$tmp[$i]}}{"read1"}}{"position"} >= ${${$sam{$tmp[$i]}}{"read2"}}{"position"}){
		$RF_num++;
		}
	}

	$total_statistic_num = $#tmp+1;
	#???????????????????????????10000?????????????????????????????????????????????UN
	if($total_statistic_num > 0){
		$percentage = int($FR_num/$total_statistic_num * 10000)/100;
		$$evidence = "FR:$FR_num(" . $percentage . "%)";
		$percentage = int($RF_num/$total_statistic_num * 10000)/100;
		$$evidence .=",RF:$RF_num(" . $percentage . "%)";
		if($FR_num/$total_statistic_num> $supportPercentage/100){
			return "FR";
		}elsif($RF_num/$total_statistic_num>$supportPercentage/100){
			return "RF";
		}else{
			return "UN";
		}
	}else{
		return "NA";
	}

SINGLE:
	open FSAM, "<$map_sam";
	$R_num = 0;
	$F_num = 0;
	while($line=<FSAM>){
		chomp($line);
		@tmp = ();
		@tmp = split(/\t/, $line);
		#It is header
		next if($#tmp<10);
		#parse alignment information
		@binarr = (0);
		$binstr=sprintf("%b", $tmp[1]);
        	for($i=length($binstr); $i>=1; $i--){
			$binarr[$#binarr+1]=substr($binstr, $i-1, 1);
	        }

		#hit multiple sites
		next if(not($tmp[$#tmp]=~/NH:i:1/));
	
		#register alignment strand
		if($binarr[5]==1){
			$R_num++;
		}else{
			$F_num++;
		}

		$num--;
		last if($num<=0);
	}
	close FSAM;

	#calculate numbers of R and F
	#????????????????????????10000??????????????????????????????????????????????????????U
	$total_statistic_num = $R_num + $F_num;
	if($total_statistic_num > 1){
		$percentage = int($F_num/$total_statistic_num * 10000)/100;
		$$evidence = "F:$F_num(" . $percentage . "%)";
		$percentage = int($R_num/$total_statistic_num * 10000)/100;
		$$evidence .=",R:$R_num(" . $percentage . "%)";
		if($F_num/$total_statistic_num>$supportPercentage/100){
			return "F";
		}elsif($R_num/$total_statistic_num>$supportPercentage/100){
			return "R";
		}else{
			return "U";
		}
	}else{
		return "NA";
	}
}

#detect sequencing read length
sub obtainReadLen{
        my ($fastqFile, $readNum) = @_;
	my ($maxReadLen, $read);
	$maxReadLen = 0;
        open FREAD, "<$fastqFile";
        for(my $i=0; $i<$readNum*4; $i++){
                last if(eof(FREAD));
		$read = <FREAD>;
		chomp($read);
                next if($i-int($i/4)*4 != 1);
		$maxReadLen = length($read) if(length($read)>$maxReadLen);
       }
	close FREAD;
	return $maxReadLen;
}

#get alignment percentage
sub obtainAlignmentPercentage{
	my ($summaryFile, $spotNum, $percentage) = @_;;
	$$percentage = 0;
	$$spotNum = 0;
	#89.01% overall alignment rate
	open SUMMARYALIGN, "<$summaryFile";
	while(my $readLine = <SUMMARYALIGN>){
		if($readLine =~/^(.*) overall alignment rate/){
			$$percentage = $1;
		}
		if($readLine =~/^(\d+) reads; of these:/){
			$$spotNum = $1;
		}
	}
	close SUMMARYALIGN;
}

sub getTime{
	my @months = qw( 01 02 03 04 05 06 07 08 09 10 11 12 );
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	return $months[$mon] . "/" . $mday . " " . $hour . ":" . $min . ":" . $sec;              
}

# ?????????????????????????????????10000B
sub getFileSize{
        my @attr = stat($_[0]);
        if($attr[7]<50000000){
                return 0;
        }else{
                return 1;
        }
}

