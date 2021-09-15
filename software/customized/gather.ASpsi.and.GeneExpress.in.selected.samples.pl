#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV < 0){
	print "This script is used to gather psi of AS in experiment and experiment sequencing information\n\n";
	print "\t perl $0 \\\n" . 
		"\t\t --psiOutputDir    ../008-pickup-psi-of-ASs-in-all-expers/psiOutputDir \\\n" . 
		"\t\t --inputExperimentList  final.selected.sample.list.tsv \\\n" .
		"\t\t --spciesAbbr  BTAU \\\n" .

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
		"\t\t --outputJcMXE      jcMXE.tsv \n" .

		"\t\t --outputGeneExpression expression.tsv \n\n";


	exit(0);
}

my ($speciesAbbr, $psiOutputDir, $inputExperimentList);
my ($outputTotalA5SS, $outputTotalA3SS, $outputTotalSE, $outputTotalRI, $outputTotalMXE);
my ($outputNovelA5SS, $outputNovelA3SS, $outputNovelSE, $outputNovelRI, $outputNovelMXE);
my ($outputJcecA5SS, $outputJcecA3SS, $outputJcecSE, $outputJcecRI, $outputJcecMXE);
my ($outputJcA5SS, $outputJcA3SS, $outputJcSE, $outputJcRI, $outputJcMXE);
my ($outputGeneExpression);

GetOptions(
        'psiOutputDir=s'=>\$psiOutputDir,
	'inputExperimentList=s'=>\$inputExperimentList,
	'speciesAbbr=s'=>\$speciesAbbr,
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
	'outputGeneExpression=s'=>\$outputGeneExpression,
);

my ($expId, $line, @fields);
my (%totalAsA5ssHash, %totalAsA3ssHash, %totalAsSeHash, %totalAsRiHash, %totalAsMxeHash);
my (%catalogA5ss, %catalogA3ss, %catalogRi, %catalogMxe, %catalogSe);
my ($a5ssNum, $a3ssNum, $riNum, $seNum, $mxeNum);
my ($asPosition, $asId);
my (@tmpAsField);

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

open WEXPRESS, ">$outputGeneExpression";
print WEXPRESS join("\t", "geneID", "geneName", "Chromo", "Strand", "Start", "End", "Coverage", "FPKM", "TPM", "ExperimentId") . "\n";
open FF, "<$inputExperimentList";
<FF>;
while($line=<FF>){
	@fields = ();
	@fields = split(/\t/, $line);
	$expId = $fields[0];

	# 将基因表达量收集到outputGeneExpression
	open FEXPRESS, "<" . $psiOutputDir . "/" . $expId . "/geneAbundanceByStringtie.tab";
	<FEXPRESS>;
	while($line=<FEXPRESS>){
		chomp($line);
		print WEXPRESS $line . "\t" . $expId . "\n";
	}
	close FEXPRESS;
	
	#register total AS and extract total AS num
	%totalAsA5ssHash = ();
	%totalAsA3ssHash = ();
	%totalAsSeHash = ();
	%totalAsRiHash = ();
	%totalAsMxeHash = ();

	############################################################################
	# 2. register AS into %catalogA5ss, %catalogA3ss ..., %catalogSe, 	   #
	# （1）根据位置信息将每个AS登记到目录哈希中，这个目录哈希是针对所有实验设置#
	#      用于登记 位置和统一编号
	# （2）对每个实验另设哈希，用于登记每个实验内部编号 和 位置关系
	############################################################################

	# read AS and save it into %catalogA5ss and %totalAsA5ssHash
	open FA5SS, "<" . $psiOutputDir . "/" . $expId . "/fromGTF.A5SS.txt";
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


	open FA3SS, "<" . $psiOutputDir . "/" . $expId . "/fromGTF.A3SS.txt";
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

	open FSE, "<" . $psiOutputDir . "/" . $expId . "/fromGTF.SE.txt";
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

	open FRI, "<" . $psiOutputDir . "/" . $expId . "/fromGTF.RI.txt";
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

	open FMXE, "<" . $psiOutputDir . "/" . $expId . "/fromGTF.MXE.txt";
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
	
	open FA5SS, "<" . $psiOutputDir . "/" . $expId . "/JCEC.raw.input.A5SS.txt";
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

	open FA3SS, "<" . $psiOutputDir . "/" . $expId . "/JCEC.raw.input.A3SS.txt";
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

	open FSE, "<" . $psiOutputDir . "/" . $expId . "/JCEC.raw.input.SE.txt";
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

	open FRI, "<" . $psiOutputDir . "/" . $expId . "/JCEC.raw.input.RI.txt";
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

	open FMXE, "<" . $psiOutputDir . "/" . $expId . "/JCEC.raw.input.MXE.txt";
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
        
	open FA5SS, "<" . $psiOutputDir . "/" . $expId . "/JC.raw.input.A5SS.txt";
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

	open FA3SS, "<" . $psiOutputDir . "/" . $expId . "/JC.raw.input.A3SS.txt";
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

	open FSE, "<" . $psiOutputDir . "/" . $expId . "/JC.raw.input.SE.txt";
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

	open FRI, "<" . $psiOutputDir . "/" . $expId . "/JC.raw.input.RI.txt";
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

	open FMXE, "<" . $psiOutputDir . "/" . $expId . "/JC.raw.input.MXE.txt";
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
close WEXPRESS;
