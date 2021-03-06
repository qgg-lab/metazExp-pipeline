#!/usr/bin/perl
use Getopt::Long;
use List::Util qw/max min/;
use strict;
if($#ARGV<0){
	print "$0 \n" . 
		" \t\t--expId ERX2073443\n" . 
		" \t\t--runIds ERR2013772,ERR2013771,ERR2013774,ERR2013773\n" . 
		" \t\t--genomeDb /mnt/home/liujind1/workAS/01-cattle/database/genome\n" .
		" \t\t--cDNAdb /mnt/home/liujind1/workAS/01-cattle/database/cDNA\n" .
		" \t\t--gtfFile /mnt/home/liujind1/workAS/01-cattle/database/annotation.gtf\n" .
		" \t\t--rMats /mnt/home/liujind1/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS2/rmats.py\n" .
		" \t\t--fastqOutputDir /mnt/scratch/liujind1/fastqOutputDir\n" .
		" \t\t--psiOutputDir /mnt/scratch/liujind1/psiOutputDir\n" .
		" \t\t--logOutputDir /mnt/scratch/liujind1/logOutputDir\n" . 
		" \t\t--threadNum 18\n";
	exit(0);
}

my ($expId, $runIdList, $fastqOutputDir, $cDNAdb, $genomeDb, $psiOutputDir, $logOutputDir, $gtfFile, $rmatsPy, $threadNum, $cmd);

GetOptions(
        'expId=s'=>\$expId,
        'runIds=s'=>\$runIdList,
	'genomeDb=s'=>\$genomeDb,
	'cDNAdb=s'=>\$cDNAdb,
	'gtfFile=s'=>\$gtfFile,
	'fastqOutputDir=s'=>\$fastqOutputDir,
	'psiOutputDir=s'=>\$psiOutputDir,
	'logOutputDir=s'=>\$logOutputDir,
	'rmatsPy=s'=>\$rmatsPy,
	'threadNum=i'=>\$threadNum,
);

#prepare workpath
$fastqOutputDir = substr($fastqOutputDir, 0, length($fastqOutputDir)-1) if(substr($fastqOutputDir, length($fastqOutputDir)-1, 1) eq "/");
$psiOutputDir = substr($psiOutputDir, 0, length($psiOutputDir)-1) if(substr($psiOutputDir, length($psiOutputDir)-1, 1) eq "/");
$logOutputDir = substr($logOutputDir, 0, length($logOutputDir)-1) if(substr($logOutputDir, length($logOutputDir)-1, 1) eq "/");

$fastqOutputDir = $fastqOutputDir . "/" . $expId if(index($fastqOutputDir, $expId)<0);
$cmd = "rm -rf $fastqOutputDir";
print "CMD: " . $cmd . "\n\n";
system($cmd);

$cmd = "mkdir -p $fastqOutputDir";
print "CMD: " . $cmd . "\n\n";
system($cmd);

$psiOutputDir = $psiOutputDir . "/" . $expId if(index($psiOutputDir, $expId)<0);
$cmd = "rm -rf $psiOutputDir";
print "CMD: " . $cmd . "\n\n";
system($cmd);

$cmd = "mkdir -p $psiOutputDir";
print "CMD: " . $cmd . "\n\n";
system($cmd);

$logOutputDir = $logOutputDir . "/" . $expId if(index($logOutputDir, $expId)<0);
$cmd = "rm -rf $logOutputDir";
print "CMD: " . $cmd . "\n\n";
system($cmd);

$cmd = "mkdir -p $logOutputDir";
print "CMD: " . $cmd . "\n\n";
system($cmd);

my (@phredScore, @phredScoreEvidence, @layout, @libraryType, @libraryTypeEvidence, @readLen);
my ($libraryTypeEvidence, $phredScoreEvidence, $phredScoreEvidence);
my (%checkPhredScoreUniq, %checkLibraryTypeUniq, %checkLayoutUniq);


my $expSeqInfo = $psiOutputDir . "/" . $expId . ".SeqInfo.txt";
open EXPERIMENTINFO, ">$expSeqInfo";
print EXPERIMENTINFO join("\t", "ExpId", "RunId", "Layout", "ReadLen", "PhredScore", "PhredDetail", "Library", "LibraryDetail") . "\n";

my @runId = split(/,/, $runIdList);
for(my $i=0; $i<=$#runId; $i++){
	#wget sra data
	my $sraAddress = "ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" . 
			substr($runId[$i], 0, 3) . "/" . 
			substr($runId[$i], 0, 6) . "/" . 
			$runId[$i] . "/" . 
			$runId[$i] . ".sra";
	
	my $cmd = "wget -O " .	$fastqOutputDir . "/" . $runId[$i] . ".sra" . 
		" --tries=10 " . 
		" -o $logOutputDir/log.o.wget." . $runId[$i] .
		" $sraAddress " .
		" 2> " . $logOutputDir . "/log.e.wget." . $runId[$i];

	print "CMD: " . $cmd . "\n\n";
	system($cmd);

	#detect sra file exists
	my @detectFileSize = stat ($fastqOutputDir . "/" . $runId[$i] . ".sra");
	if($detectFileSize[7]==0){
		print STDERR "Can't download " . $runId[$i] . ".sra from sraDB online\n\n";
		next;
	}

	#dump data
	my $cmd = "fastq-dump " . 
		  " --split-files " . 
		  " -O $fastqOutputDir " . 
		  " " .  $fastqOutputDir . "/" . $runId[$i] . ".sra" . 
		  " 1> $logOutputDir/log.o.dump.$runId[$i] " . 
		  " 2> $logOutputDir/log.e.dump.$runId[$i]";
	print "CMD: " . $cmd . "\n\n";
	system($cmd);

	#remove sra
	my $cmd = "rm -rf " . $fastqOutputDir . "/" . $runId[$i] . ".sra" ;
	print "CMD: " . $cmd . "\n\n";
	system($cmd);

	#detect dump results
	my @detectFileSize = stat ("$fastqOutputDir/log.e.dump.$runId[$i]");
	if($detectFileSize[7]!=0){
		print STDERR "Can't dump " . $runId[$i] . "fastq data from sra file.\n\n";
		next;
	}

	#obtain phred score in any one run data
	if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq"){
		$phredScore[$i] = &calPhredScore($fastqOutputDir . "/" . $runId[$i] . "_1.fastq", 10000, \$phredScoreEvidence);
	}else{
		$phredScore[$i] = &calPhredScore($fastqOutputDir . "/" . $runId[$i] . "_2.fastq", 10000, \$phredScoreEvidence);
	}

	#obtain layout
	$layout[$i] = "SINGLE";
	if(-e("$fastqOutputDir/$runId[$i]" . "_1.fastq") and -e ("$fastqOutputDir/$runId[$i]" . "_2.fastq")){
		$layout[$i] = "PAIRED";
	}

	#generate tmp fastq file
	if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq"){
		my $cmd = "head -n 4000000 $fastqOutputDir" . "/" . $runId[$i] . "_1.fastq " .
		"> $fastqOutputDir" . "/" . $runId[$i] . "_1.fq.tmp " .
		"2> " . $logOutputDir .  "/log.e.head100000Read.from.1.fq." . $runId[$i];
		print "CMD: " . $cmd . "\n\n";
		system($cmd);
	}

	if( -e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq"){
		my $cmd = "head -n 4000000 $fastqOutputDir" . "/" . $runId[$i] . "_2.fastq " .
		"> $fastqOutputDir" . "/" . $runId[$i] . "_2.fq.tmp " .
		"2> " . $logOutputDir .  "/log.e.head100000Read.from.2.fq." . $runId[$i];
		print "CMD: " . $cmd . "\n\n";
		system($cmd) if( -e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq");
	}

	if($layout[$i] eq "PAIRED"){
		my $cmd = "hisat2 -x $cDNAdb" .
			  " -1 " . $fastqOutputDir . "/" . $runId[$i] . "_1.fq.tmp" . 
			  " -2 " . $fastqOutputDir . "/" . $runId[$i] . "_2.fq.tmp" . 
			  " --phred" . $phredScore[$i] . 
			  " --no-unal" .
			  " -p " . $threadNum . 
			  " -S $fastqOutputDir" . "/" . "$runId[$i].map.cDNA.sam " .  
			  " 2> " . $logOutputDir .  "/log.e.hisat2.map.cDNA." . $runId[$i];
		print "CMD: " . $cmd . "\n\n";
		system($cmd);

	}else{
		my $cmd = "hisat2 -x $cDNAdb " . 
			  " -U " . $fastqOutputDir . "/" . $runId[$i] . "_1.fq.tmp" . 
			  " --phred" . $phredScore[$i] . 
			  " -p " . $threadNum . 
			  " --no-unal" .
			  " -S $fastqOutputDir" . "/" . "$runId[$i].map.cDNA.sam" .  
			  " 2> " . $logOutputDir .  "/log.e.hisat2.map.cDNA." . $runId[$i];
		if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq"){
			print "CMD: " . $cmd . "\n\n";
			system($cmd);
		}

		my $cmd = "hisat2 -x $cDNAdb " . 
			  " -U " . $fastqOutputDir . "/" . $runId[$i] . "_2.fq.tmp" . 
			  " --phred" . $phredScore[$i] . 
			  " -p " . $threadNum . 
			  " --no-unal" .
			  " -S $fastqOutputDir" . "/" . "$runId[$i].map.cDNA.sam" .  
			  " 2> " . $logOutputDir .  "/log.e.hisat2.map.cDNA." . $runId[$i];
		if( -e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq"){
			print "CMD: " . $cmd . "\n\n";
			system($cmd);
		}
	}

	#obtain library type
	$libraryType[$i] = &obtainLibraryType(	$fastqOutputDir. "/" . $runId[$i] . ".map.cDNA.sam", 
						1000000,
						$layout[$i], 
						85, 
						\$libraryTypeEvidence);

	#get read length
	if( -e $fastqOutputDir. "/" . $runId[$i] . "_1.fastq"){
		$readLen[$i] = &obtainReadLen(	$fastqOutputDir. "/" . $runId[$i] . "_1.fastq", 1000000);
	}else{
		$readLen[$i] = &obtainReadLen(  $fastqOutputDir. "/" . $runId[$i] . "_2.fastq", 1000000);
	}

	#output run data information
	print EXPERIMENTINFO join("\t", 
				$expId, $runId[$i], 
				$layout[$i],
				$readLen[$i], 
				$phredScore[$i], 
				$phredScoreEvidence,
				$libraryType[$i], 
				$libraryTypeEvidence) . "\n";

	print join("\t", "ExpId", "RunId", "Layout", "ReadLen", "PhredScore", "PhredDetail", "Library", "LibraryDetail") . "\n";
	print join("\t", 
				$expId, $runId[$i], 
				$layout[$i],
				$readLen[$i], 
				$phredScore[$i], 
				$phredScoreEvidence,
				$libraryType[$i], 
				$libraryTypeEvidence) . "\n\n";


	#delete temporary fastq data file	
	if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fq.tmp"){
		my $cmd = "rm -rf " . "$fastqOutputDir" . "/" . $runId[$i] . "_1.fq.tmp";
		print "CMD: " . $cmd . "\n\n";
		system($cmd);
	}

	if(-e "$fastqOutputDir" . "/" . $runId[$i] . "_2.fq.tmp"){
		$cmd = "rm -rf " . "$fastqOutputDir" . "/" . $runId[$i] . "_2.fq.tmp";
		print "CMD: " . $cmd . "\n\n";
		system($cmd);
	}

	#delete temporary alignment file against to cDNA database
	my $cmd = "rm -rf " . $fastqOutputDir. "/" . $runId[$i] . ".map.cDNA.sam";
	print "CMD: " . $cmd . "\n\n"; 
	system($cmd);
	
	#??????phredscore???libraryType, layout?????????
	$checkPhredScoreUniq{$phredScore[$i]} = 1;
	$checkLibraryTypeUniq{$libraryType[$i]} = 1 if($libraryType[$i] ne "NA");
	$checkLayoutUniq{$layout[$i]} = 1;
}

#????????????run???phred???layout???libraryType????????????
my (@keyPhredScore, @keyLibraryType, @keyLayout);
@keyPhredScore = keys(%checkPhredScoreUniq);
@keyLibraryType = keys(%checkLibraryTypeUniq);
@keyLayout = keys(%checkLayoutUniq);
if($#keyPhredScore > 0 ){
	print STDERR "PhredScores of Runs in Experiment are not consistent.\n\n";
}
if($#keyLibraryType > 0){
	print STDERR "Library type of Runs in Experiment are not consistent.\n\n";
}
if($#keyLayout > 0){
	print STDERR "Layout of Runs in Experiment are not consistent.\n\n";
}
if($#keyPhredScore > 0 or $#keyLibraryType > 0 or $#keyLayout > 0){
	exit;
}


#generate phredscore, fastqfiles and libraryType parameters
my ($fastqParameter, $fastqParameter1, $fastqParameter2, $phredScoreParameter, $libraryTypeParameter);
#phred
$phredScoreParameter = "--phred" . $keyPhredScore[0];

#rna-strandness
$libraryTypeParameter = "";
if($#keyLibraryType==0){
	if($keyLibraryType[0] ne "UN" and $keyLibraryType[0] ne "U"){
		$libraryTypeParameter = "--rna-strandness " . $keyLibraryType[0];
	}
}elsif($#keyLibraryType==-1){
	$libraryTypeParameter = "";
}

#fastqFile
$fastqParameter1 = "";
$fastqParameter2 = "";
for(my $i=0; $i<=$#runId; $i++){
	$fastqParameter1 .= $fastqOutputDir . "/" . $runId[$i] . "_1.fastq," if(-e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq");
	$fastqParameter2 .= $fastqOutputDir . "/" . $runId[$i] . "_2.fastq," if(-e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq");
}

$fastqParameter1=substr($fastqParameter1, 0, length($fastqParameter1)-1) if( -e $fastqOutputDir . "/" . $runId[0] . "_1.fastq");
$fastqParameter2=substr($fastqParameter2, 0, length($fastqParameter2)-1) if( -e $fastqOutputDir . "/" . $runId[0] . "_2.fastq");

if($layout[0] eq "SINGLE"){
	$fastqParameter = "-U " . $fastqParameter1 if( -e $fastqOutputDir . "/" . $runId[0] . "_1.fastq");
	$fastqParameter = "-U " . $fastqParameter2 if( -e $fastqOutputDir . "/" . $runId[0] . "_2.fastq");
}else{
	$fastqParameter = "-1 " . $fastqParameter1 . " -2" . $fastqParameter2;	
}


#call hisat2 and samtools to generate sorted bam alignment file
my $cmd = "hisat2 -x $genomeDb " . 
	$fastqParameter . " " . 
	$phredScoreParameter . " " . 
	$libraryTypeParameter . 
	" -p " . $threadNum . 
	" --novel-splicesite-outfile " . $psiOutputDir  . "/" . $expId . ".novel.splicesite.from.hisat.txt " .
	" --summary-file " . $psiOutputDir  . "/" . $expId . ".alignSum.from.hisat.txt " . 
	" 2>" . $logOutputDir .  "/log.e.map." . $expId . ".to.genome" . 
	" | samtools view -bS - " . 
	" 2> $logOutputDir/" . "log.e.samView." . $expId . 
	" | samtools sort - -o " . $fastqOutputDir . "/" . $expId . ".bam" . 
	" 2>" . $logOutputDir . "/log.e.samSort." . $expId;

print "CMD: " . $cmd . "\n\n";
system($cmd);

#delete fastq file
for(my $i=0; $i<=$#runId; $i++){
	if( -e $fastqOutputDir . "/" . $runId[$i] . "_1.fastq"){
		$cmd = "rm -rf " . $fastqOutputDir . "/" . $runId[$i] . "_1.fastq";
		print "CMD: " . $cmd . "\n\n";
        	system($cmd);
	}

	if(-e $fastqOutputDir . "/" . $runId[$i] . "_2.fastq"){
		$cmd = "rm -rf " . $fastqOutputDir . "/" . $runId[$i] . "_2.fastq";
		print "CMD: " . $cmd . "\n\n";
		system($cmd);
	}
}


#output alignment percentage into experiment information file
my ($totalSpotNum, $alignmentPercentage);
&obtainAlignmentPercentage($psiOutputDir  . "/" . $expId . ".alignSum.from.hisat.txt", \$totalSpotNum, \$alignmentPercentage);
print EXPERIMENTINFO "TotalSpot:$totalSpotNum\tAlignPercent:$alignmentPercentage\n";
print "TotalSpot:$totalSpotNum\tAlignPercent:$alignmentPercentage\n\n";

#library parameter for rMATS
my $expLayout = lc($keyLayout[0]);

my $expLibraryType = "";
if($#keyLibraryType == 0){
	$expLibraryType = $keyLibraryType[0];
	if($expLibraryType eq "RF" or $expLibraryType eq "R" ){
		$expLibraryType = "fr-firststrand";
	}elsif($expLibraryType eq "FR" or $expLibraryType eq "F"){
		$expLibraryType = "fr-secondstrand";
	}elsif($expLibraryType eq "UN" or $expLibraryType eq "U"){
		$expLibraryType = "fr-unstranded";
	}
}elsif($#keyLibraryType == -1){
	$expLibraryType = "fr-unstranded";
}

my $readLength = max @readLen;


#assemble transcripts with stringTie
my $stringtieLibraryParameter = "";
if($expLibraryType eq "fr-firststrand"){
	$stringtieLibraryParameter = "--rf";
}elsif($expLibraryType eq "fr-secondstrand"){
	$stringtieLibraryParameter = "--fr";
}else{
	$stringtieLibraryParameter = "";
}

my $cmd = "stringtie " . $fastqOutputDir . "/" . $expId . ".bam" .
	  " -G " . $gtfFile . 
	  " $stringtieLibraryParameter " . 
	  " -o " . $psiOutputDir . "/" . "transcriptomeByStringtie.gtf" . 
	  " -C " . $psiOutputDir . "/" . "annotationCovByStringtie.gtf" .
	  " -A " . $psiOutputDir . "/" . "geneAbundanceByStringtie.tab" .
	  " -p " . $threadNum . 
	  " 2> " . $logOutputDir . "/log.e.stringtie.$expId";
print "CMD: " . $cmd . "\n\n";
system($cmd);

#generate b file to save bam file
open RMATSB, ">" . $fastqOutputDir . "/" . $expId . ".rMats.b.txt";
print RMATSB $fastqOutputDir . "/" . $expId . ".bam";
print "CMD: " . "echo \"" . $fastqOutputDir . "/" . $expId . ".bam" . "\" > " . $fastqOutputDir . "/" . $expId . ".rMats.b.txt" . "\n\n";
close RMATSB;

my $cmd = "python " . 
	  " $rmatsPy " .
	  " --b1 " . $fastqOutputDir . "/" . $expId . ".rMats.b.txt" . 
	  " --gtf " . $gtfFile .
	  " --od " . $psiOutputDir .
	  " -t " . $expLayout . 
	  " --nthread " . $threadNum . 
	  " --tstat " . $threadNum . 
	  " --libType " . $expLibraryType . 
	  " --readLength " . $readLength .
	  " --statoff " . 
	  " 1> " . $logOutputDir . "/log.o.rMats.$expId" . 
	  " 2> " . $logOutputDir . "/log.e.rMats.$expId";

print "CMD: " . $cmd . "\n\n";
system($cmd);

$cmd = "rm -rf " . $fastqOutputDir;
print "CMD: " . $cmd . "\n\n";
system($cmd);

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

#detect library type such as RF, FR, UN, R, F, U
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

