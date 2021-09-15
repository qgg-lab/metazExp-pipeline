#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV <0){
	print "\n\tperl $0 \\\n" . 
	" \t\t--cmpRltListFile /mnt/home/liujind1/workAS/01-cattle/" . 
		"004-filter-trspts-with-Exp.Study/" . 
		"01.trsptMinLen200TrsptMinCov2ExonMinCov0/compared.Rlt.list.txt \\\n" .
	" \t\t--sampleInfoFile /mnt/home/liujind1/workAS/01-cattle/" . 
		"000-prepare-local-datasource/sample.information.tsv \\\n" . 
	" \t\t--trsptOriginFile /mnt/home/liujind1/workAS/01-cattle/" . 
		"004-filter-trspts-with-Exp.Study/" . 
		"01.trsptMinLen200TrsptMinCov2ExonMinCov0/originNameMappingTrsptId.tsv \\\n" .
	" \t\t--outputListOfTrsptFile /mnt/home/liujind1/workAS/01-cattle/" . 
		"004-filter-trspts-with-Exp.Study/01.trsptMinLen200TrsptMinCov2ExonMinCov0/" . 
		"combined.trspts.with.experiment.and.study.info.tsv \n";

	print "\n\t\tThis script is used to tag experiment and study num with combined trspts.\n" .
		"\t\tExperiment and study information from cmpRltListFile and sampleInfoFile.\n" . 
		"\t\ttrsptOriginFile is used to tag origin for downstream filtering operation.\n\n";
	exit(0);
}

my ($cmpRltListFile, $sampleInfoFile, $trsptOriginFile, $outputListOfTrsptFile);

GetOptions(
	'cmpRltListFile=s'=>\$cmpRltListFile, 
	'sampleInfoFile=s'=>\$sampleInfoFile,
	'trsptOriginFile=s'=>\$trsptOriginFile,
	'outputListOfTrsptFile=s'=>\$outputListOfTrsptFile,
);

#read trspt from trsptOriginFile into hash
my (%trspt, @trsptLine);
open TRSPT, "<$trsptOriginFile";
#StringTie       SRX2160812.1.1
@trsptLine=<TRSPT>;
close TRSPT;

foreach my $trsptLine(@trsptLine){
	chomp($trsptLine);
	my @tmp = ();
	@tmp = split(/\t/, $trsptLine);
	${$trspt{$tmp[1]}}{"origin"} = $tmp[0];
}

# rea dexperimeni and study mapping information from sample information file and then load into hasu
my (%experiment);
my ($runLine, @field);
open SAMPLE, "<$sampleInfoFile";
# $7|___|$10|___|$11|___|$19}
while($runLine=<SAMPLE>){
        @field = ();
        @field = split(/\|___\|/, $runLine);
        $experiment{$field[18]} = $field[57];
}
close FF;

# gather each cmpRlt file. 
my (@cmpRltList, $cmpRltFile, $experimentId);
open CMPRLTLIST, "<$cmpRltListFile";
@cmpRltList=<CMPRLTLIST>;
close CMPRLTLIST;
#/mnt/home/liujind1/workAS/01-cattle/004-filter-trspts-with-Exp.Study/01.trsptMinLen200TrsptMinCov2ExonMinCov0/assemblyDir/cmp.DRX001353.tmap
foreach $cmpRltFile(@cmpRltList){

	# Identify experiment Id from cmpRlt file name.
	my @tmp = ();
	@tmp = split(/\//, $cmpRltFile);
	if($tmp[$#tmp]=~/cmp\.(.*)\.tmap/){
		$experimentId = $1;
	}

	# open cmpRlt file to tag experiment and study information into trspt hash
	# FAM243A ENSBTAT00000018387      o       SRX1149755.1    SRX1149755.1.1  100     0.497800 ...
	open CMPRLT, "<$cmpRltFile";
	while(my $cmpRltLine=<CMPRLT>){
		my @tt = ();
		@tt = split(/\t/, $cmpRltLine);
		
		# discard any line without "=" code
		next if($tt[2] ne "=");
		
		# For this trspt, experiment num increase 1
		${$trspt{$tt[1]}}{"experimentNum"}++;
		# Add experiment Id into this trspt
		${$trspt{$tt[1]}}{"experimentList"}.=$experimentId . ",";
		
		# If the study of this experiment has been added into this trspt, nothing should to be do.
		# Only new study should be added into trspt.
		if(index(${$trspt{$tt[1]}}{"studyList"}, $experiment{$experimentId}) < 0){
			${$trspt{$tt[1]}}{"studyNum"}++;
			${$trspt{$tt[1]}}{"studyList"}.=$experiment{$experimentId} . ",";
		}		
	}
	close CMPRLT;
}

# Output trspt with origin, experiment and study num, as well as experiment and study list
open WW, ">$outputListOfTrsptFile";
my @trsptId = keys(%trspt);
foreach my $trsptId(@trsptId){
	print WW join("\t", $trsptId, ${$trspt{$trsptId}}{"origin"}, ${$trspt{$trsptId}}{"experimentNum"}, ${$trspt{$trsptId}}{"studyNum"}, ${$trspt{$trsptId}}{"experimentList"}, ${$trspt{$trsptId}}{"studyList"}) . "\n";
}
close WW;
