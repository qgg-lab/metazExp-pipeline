#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw/max min/;
use strict;
if($#ARGV < 0){
	print   "\n\nThis program is used to generate sbatch scripts.\n" . 
		"The program will only generate specified experiment\n sbatch script with experimentId parameter\n" . 
		"\n\n";
	
	print "$0 \\\n" .
		"--sampleInfoTable curated.experiment.tsv \\\n" .
		"--psiScript /mnt/home/liujind1/software/customized/psi.in.experiment.pl \\\n" .
		"--rMATSscript /mnt/home/liujind1/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS2/rmats.py \\\n" . 
		"--switchPsi  ON \\\n" .
		"--genomeDb /mnt/home/liujind1/workAS/01-cattle/database/genome \\\n" . 
		"--cDNAdb /mnt/home/liujind1/workAS/01-cattle/database/cDNA \\\n" . 
		"--gtfFile /mnt/home/liujind1/workAS/01-cattle/database/annotation.gtf \\\n" . 
		"--maxIntronSize 12000 \\\n" .
		"--minIntronSize 20 \\\n" .
		"--switchAssemble on \\\n" .
		"--fastqOutputDir /mnt/scratch/liujind1/fastqOutputDir \\\n" .
		"--logOutputDir /mnt/home/liujind1/workAS/01-cattle/001-gather-experiment-psi/logOutputDir \\\n" . 
		"--psiOutputDir /mnt/home/liujind1/workAS/01-cattle/001-gather-experiment-psi/psiOutputDir \\\n" .
		"--manualDownloadDir /mnt/home/liujind1/workAS/01-cattle/003-reAsmbl-failed-goodExps/manualDownloadDir \\\n" .
		"--taskOutputDir /mnt/home/liujind1/workAS/01-cattle/001-gather-experiment-psi/taskOutputDir \\\n" . 
		"--taskRunTimeLogDir /mnt/home/liujind1/workAS/01-cattle/001-gather-experiment-psi/taskRunTime \\\n" . 

		"--nodes 1 \\\n" . 
		"--ntasksPerNode 1 \\\n" . 
		"--cpusPerTask 18 \\\n" .
		"--memorySize 64G \\\n" . 
		"--pythonDir /mnt/home/liujind1/software/Python-2.7.13/bin \\\n" . 
		"--pythonLibDir /mnt/home/liujind1/software/Python-2.7.13/lib \\\n" .
		"--hisat2Dir /mnt/home/liujind1/software/hisat2-2.1.0 \\\n" . 
		"--stringtie1Dir /mnt/research/qgg/ding/liujind2/software/stringtie-1.3.5 \\\n" . 
		"--sratoolsDir /mnt/home/liujind1/software/sraToolkit-2.9.4-1/bin \\\n" . 
		"--samtoolsDir /mnt/home/liujind1/software/samtools-1.9/bin \\\n" . 
		"--lapackLibDir /mnt/home/liujind1/software/lapack-3.8.0 \\\n" . 
		"--blasLibDir /mnt/home/liujind1/software/BLAS-3.8.0 \\\n" .
		"--runMinutesPerMillion 3 \\\n" . 

		"--sraLocalCacheDir /mnt/home/liujind1/ncbi/public/sra \\\n" .
		
		"--outputSbatchDir /mnt/home/liujind1/workAS/01-cattle/001-gather-experiment-psi/sbatchs\n";
		exit;
}
my ($sampleInfoTable, $minSpotMillion, $maxSpotMillion);
my ($expId, $taskPrefix);
my ($psiScript, $rMATSscript, $switchPsi, $switchAssemble);
my ($genomeDb, $cDNAdb, $gtfFile);
my ($subSeqIdListDir, $subSeqIdListPrefix, $subGtfDir, $subGtfPrefix, $callAsFromBam);
my ($maxIntronSize, $minIntronSize);
my ($logOutputDir, $fastqOutputDir, $psiOutputDir, $taskOutputDir, $taskRunTimeLogDir, $manualDownloadDir);
my ($nodes, $ntasksPerNode, $cpusPerTask, $memorySize);
my ($pythonDir, $hisat2Dir, $sratoolsDir, $samtoolsDir, $stringtie1Dir);
my ($pythonLibDir, $lapackLibDir, $blasLibDir);
my $runMinutesPerMillion;
my $sraLocalCacheDir;
my $outputSbatchDir;

GetOptions(
	'sampleInfoTable=s'=>\$sampleInfoTable,
	'taskPrefix=s'=>\$taskPrefix,
	'psiScript=s'=>\$psiScript,
	'rMATSscript=s'=>\$rMATSscript,
	'callAsFromBam=s'=>\$callAsFromBam,
	'switchPsi=s'=>\$switchPsi,
	'switchAssemble=s'=>\$switchAssemble,
	'genomeDb=s'=>\$genomeDb,
	'cDNAdb=s'=>\$cDNAdb,
	'gtfFile=s'=>\$gtfFile,
	'maxIntronSize=s'=>\$maxIntronSize,
	'minIntronSize=s'=>\$minIntronSize,
	'logOutputDir=s'=>\$logOutputDir,
	'fastqOutputDir=s'=>\$fastqOutputDir,
        'psiOutputDir=s'=>\$psiOutputDir,
	'taskOutputDir=s'=>\$taskOutputDir,
	'taskRunTimeLogDir=s'=>\$taskRunTimeLogDir,
	'manualDownloadDir=s'=>\$manualDownloadDir,
	'nodes=i'=>\$nodes,
	'ntasksPerNode=i'=>\$ntasksPerNode,
	'cpusPerTask=i'=>\$cpusPerTask,
	'memorySize=s'=>\$memorySize,
	'pythonDir=s'=>\$pythonDir,
	'hisat2Dir=s'=>\$hisat2Dir,
	'stringtie1Dir=s'=>\$stringtie1Dir,
	'sratoolsDir=s'=>\$sratoolsDir,
	'samtoolsDir=s'=>\$samtoolsDir,
	'pythonLibDir=s'=>\$pythonLibDir,
	'lapackLibDir=s'=>\$lapackLibDir,
	'blasLibDir=s'=>\$blasLibDir,
	'sraLocalCacheDir=s'=>\$sraLocalCacheDir,
	'runMinutesPerMillion=i'=>\$runMinutesPerMillion,
	'outputSbatchDir=s'=>\$outputSbatchDir,
	'subSeqIdListDir=s'=>\$subSeqIdListDir, 
	'subSeqIdListPrefix=s'=>\$subSeqIdListPrefix, 
	'subGtfDir=s'=>\$subGtfDir, 
	'subGtfPrefix=s'=>\$subGtfPrefix
);

#read each experiment into hash
my (%experiment, $experimentId, $i);
my ($experimentLine, @field, $field, @fieldName, $fieldName);

open FF, "<$sampleInfoTable";
# Experiment	RunList	Study	Tissue	SubTissue	Development	Treatment	Ecotype	Cultivar	Genotype	Base	Layout	SpotsNum	ReadNum	SpotLen	ReadLen	Gather	AS	Assemble
# Ecotype Cultivar        Genotype        Tissue  SubTissue       Development     Treatment       Experiment      Study   Base    Layout  SpotsNum       ReadNum SpotLen ReadLen Gather  AS      Assemble        RunList Phenotype
$experimentLine=<FF>;
chomp($experimentLine);
@fieldName = split(/\t/, $experimentLine);
my $exptId = 0;
for($exptId=0; $exptId<=$#fieldName; $exptId++){
	last if($fieldName[$exptId] eq "Experiment");
}
while($experimentLine=<FF>){
	chomp($experimentLine);
	@field = ();
	@field = split(/\t/, $experimentLine);
	for($i=0; $i<=$#field; $i++){
		${$experiment{$field[$exptId]}}{$fieldName[$i]} = $field[$i];
	}
}
close FF;

my ($sbatchText, $minutesNum);
my @experimentId = keys(%experiment);

foreach $experimentId (@experimentId){

	$sbatchText= "#!/bin/bash\n";	
	$sbatchText.= "#SBATCH --job-name=" . $taskPrefix . "_" . "$experimentId\n";
	$sbatchText.= "#SBATCH --nodes=$nodes\n";
	$sbatchText.= "#SBATCH --ntasks-per-node=$ntasksPerNode\n";
	$sbatchText.= "#SBATCH --cpus-per-task=$cpusPerTask\n";
	$sbatchText.= "#SBATCH --mem=$memorySize\n";
	$minutesNum = ${$experiment{$experimentId}}{"SpotsNum"} * $runMinutesPerMillion + 3;
	$sbatchText.= "#SBATCH --time=" . sprintf("%02d", int($minutesNum/60)) . ":" . sprintf("%02d", $minutesNum - int($minutesNum/60) * 60)  . ":00\n\n\n";
	$sbatchText.= "export PATH=" . $hisat2Dir . ":\$PATH\n";	
	$sbatchText.= "export PATH=" . $stringtie1Dir . ":\$PATH\n";
	$sbatchText.= "export PATH=" . $sratoolsDir . ":\$PATH\n";
	$sbatchText.= "export PATH=" . $samtoolsDir . ":\$PATH\n";

	$sbatchText.= "export PATH=" . $pythonDir . ":\$PATH\n";
	$sbatchText.= "export LD_LIBRARY_PATH=" . $pythonLibDir . ":\$LD_LIBRARY_PATH\n";
	$sbatchText.= "export LD_LIBRARY_PATH=" . $lapackLibDir . ":\$LD_LIBRARY_PATH\n";
	$sbatchText.= "export LAPACK=" .  $lapackLibDir . "\n";
	$sbatchText.= "export LD_LIBRARY_PATH=" . $blasLibDir . ":\$LD_LIBRARY_PATH\n";
	$sbatchText.= "export export BLAS=" . $blasLibDir . "/libfblas.a\n";
	$sbatchText.= "export outputPath=" . $fastqOutputDir . "\n";
	$sbatchText.= "export threadNum=" . $cpusPerTask . "\n\n\n";

	$sbatchText.= "date > " . $taskRunTimeLogDir . "/runTime." . $experimentId . "\n\n\n";

	# 生成启动psi的程序的命令，放在sbatch中
	$sbatchText.= "perl " . $psiScript . 
		" --expId " . $experimentId . 
		" --runIds " . ${$experiment{$experimentId}}{"RunList"} . 
		" --genomeDb " . $genomeDb . 
		" --cDNAdb " . $cDNAdb . 
		" --gtfFile " . $gtfFile . 
		" --maxIntronSize " . $maxIntronSize . 
		" --minIntronSize " . $minIntronSize . 
		" --rMats " . $rMATSscript . 
		" --switchPsi " . $switchPsi . 
		" --switchAssemble " . $switchAssemble . 
		" --fastqOutputDir " . $fastqOutputDir . 
		" --psiOutputDir " . $psiOutputDir . 
		" --logOutputDir " . $logOutputDir . 
		" --manualDownloadDir " . $manualDownloadDir . 
		" --threadNum $cpusPerTask 1> " . $taskOutputDir . "/log.o." . $experimentId . ".cmd.sh " .
		" 2> " . $taskOutputDir . "/log.e." . $experimentId . "\n\n\n";

	

	# 获得psiOutput目录下的:DRX081307.SeqInfo.txt中提取layout, readLen和libraryType等信息
	my $alignmentFile = $psiOutputDir . "/" .  $experimentId . "/" . $experimentId . ".SeqInfo.txt";

	# 生成启动call.as.from.bam.pl的命令
	$sbatchText.= "perl " . $callAsFromBam .
		" --rmats " . $rMATSscript .
		" --subSeqIdListDir " . $subSeqIdListDir .
		" --subSeqIdListPrefix " . $subSeqIdListPrefix .
		" --subGtfDir " . $subGtfDir .
		" --subGtfPrefix " . $subGtfPrefix .
		" --specifiedSubGroupRegion 0-17 " .
		" --threadNum " . $cpusPerTask . 
		" --inputBamFile " . $fastqOutputDir . "/" . $experimentId . "/" . $experimentId . ".bam" .
		" --outputDir " . $psiOutputDir . "/" . $experimentId .
		" --outputSubPrefix sub " .
		" --seqInfoFile " . $alignmentFile .
		" 1>" . $psiOutputDir . "/" . $experimentId . "/log.o.txt " .
		" 2>" . $psiOutputDir . "/" . $experimentId . "/log.e.txt\n\n\n";

	$sbatchText.= "perl " . $callAsFromBam .
		" --rmats " . $rMATSscript .
		" --subSeqIdListDir " . $subSeqIdListDir .
		" --subSeqIdListPrefix " . $subSeqIdListPrefix .
		" --subGtfDir " . $subGtfDir .
		" --subGtfPrefix " . $subGtfPrefix .
		" --specifiedSubGroupRegion 18-34 " .
		" --threadNum " . $cpusPerTask . 
		" --inputBamFile " . $fastqOutputDir . "/" . $experimentId . "/" . $experimentId . ".bam" .
		" --outputDir " . $psiOutputDir . "/" . $experimentId .
		" --outputSubPrefix sub " .
		" --seqInfoFile " . $alignmentFile .
		" 1>" . $psiOutputDir . "/" . $experimentId . "/log.o.txt " .
		" 2>" . $psiOutputDir . "/" . $experimentId . "/log.e.txt\n\n\n";

	$sbatchText.= "perl " . $callAsFromBam .
		" --rmats " . $rMATSscript .
		" --subSeqIdListDir " . $subSeqIdListDir .
		" --subSeqIdListPrefix " . $subSeqIdListPrefix .
		" --subGtfDir " . $subGtfDir .
		" --subGtfPrefix " . $subGtfPrefix .
		" --specifiedSubGroupRegion 35-49 " .
		" --threadNum " . $cpusPerTask . 
		" --inputBamFile " . $fastqOutputDir . "/" . $experimentId . "/" . $experimentId . ".bam" .
		" --outputDir " . $psiOutputDir . "/" . $experimentId .
		" --outputSubPrefix sub " .
		" --seqInfoFile " . $alignmentFile .
		" 1>" . $psiOutputDir . "/" . $experimentId . "/log.o.txt " .
		" 2>" . $psiOutputDir . "/" . $experimentId . "/log.e.txt\n\n\n";



# 收集分散在sub中的可变剪接事件
	$sbatchText.= "perl /mnt/home/liujind1/software/customized/gather.as.from.subgroupRegion.pl " .
		" --subSeqIdListDir " . $subSeqIdListDir .
		" --subSeqIdListPrefix " . $subSeqIdListPrefix .
		" --subGtfDir " . $subGtfDir .
		" --subGtfPrefix " . $subGtfPrefix .
		" --outputDir " . $psiOutputDir . "/" . $experimentId .
		" --outputSubPrefix sub " .
		" --seqInfoFile " . $alignmentFile .
		" --subGroupRegion 0-49 " .
		" 1>" . $psiOutputDir . "/" . $experimentId . "/log.o.txt " .
		" 2>" . $psiOutputDir . "/" . $experimentId . "/log.e.txt\n\n\n";




	# 删除fastqOutpuDir
	$sbatchText.= "rm -rf " .  $fastqOutputDir . "/" . $experimentId . "\n\n\n";

	#delete run sra cache file
	my @runId = ();
	@runId = split(/,/, ${$experiment{$experimentId}}{"runIdList"});

	#If cache sra file
	for(my $j = 0; $j<=$#runId; $j++){
		$sbatchText.= "rm -rf " . $sraLocalCacheDir . "/" . $runId[$j] . ".sra\n";
	}

	$sbatchText.= "\ndate >> " . $taskRunTimeLogDir . "/runTime." . $experimentId . "\n";

	open WW, ">" . $outputSbatchDir . "/" . $experimentId;
	print WW $sbatchText;
	close WW;	
}

