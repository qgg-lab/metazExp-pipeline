#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV < 0){
	print "This script is used to gather transcriptome to build triplet annotationGTF\n\n";
	print "perl $0 \\\n" . 
		"\t\t --psiOutputDir      /mnt/home/liujind1/workAS/01-cattle/001-gather-experiment-psi/psiOutputDir \\\n" . 
		"\t\t --expIdList         /mnt/home/liujind1/workAS/01-cattle/002-build-tripletAnnotation/experimentIdList.txt \\\n" . 
		"\t\t --totalSopts        40 \\\n" .
		"\t\t --alignPer          50 \\\n" .
		"\t\t --mappedSpots       40 \\\n" .
		"\t\t --layout            PAIRED,SINGLE \\\n" .
		"\t\t --minReadLength        50 \\\n" .
		"\t\t --libraryType       RF,FR,UN,R,F,U \\\n" . 
		"\t\t --stringtie         /mnt/home/liujind1/software/stringtie-1.3.5/stringtie \\\n" .
		"\t\t --basicRefAnno       /mnt/home/liujind1/workAS/01-cattle/database/ensembl.annotation.gtf \\\n" .
		"\t\t --secondRefAnno       /mnt/home/liujind1/workAS/01-cattle/database/refSeq.annotation.with_ensemblSeqId.gtf \\\n" .
		"\t\t --mergeCov          5 \\\n" .
		"\t\t --outputSelectedGtfList /mnt/home/liujind1/workAS/01-cattle/002-build-tripletAnnotation/selectedGtfList.txt \\n" .
		"\t\t --outputTripletAnno /mnt/home/liujind1/workAS/01-cattle/002-build-tripletAnnotation/triplet.annotation.gtf \n";
	exit(0);
}

my ($psiOutputDir, $expIdList, $outputExperimentInfo, $outputSelectedGtfList);
my ($basicRefAnno, $secondRefAnno, $totalSpots, $alignPer, $mappedSpots);
my ($layout, $minReadLength, $outputTripletAnno, $libraryType);
my ($mergeCov, $stringtie);
my ($expId, @expDir);
my (%layoutType, %readLength, %libraryType);
my ($totalSpotsInExp, $alignPerInExp);


GetOptions(
        'psiOutputDir=s'=>\$psiOutputDir,
        'expIdList=s'=>\$expIdList,
	'totalSopts=i'=>\$totalSpots,
	'alignPer=i'=>\$alignPer,
	'mappedSpots=i'=>\$mappedSpots,
	'layout=s'=>\$layout,
	'minReadLength=i'=>\$minReadLength,
	'libraryType=s'=>\$libraryType,
	'basicRefAnno=s'=>\$basicRefAnno,
	'secondRefAnno=s'=>\$secondRefAnno,
	'stringtie=s'=>\$stringtie,
	'mergeCov=s'=>\$mergeCov,
        'outputTripletAnno=s'=>\$outputTripletAnno,
	'outputSelectedGtfList=s'=>\$outputSelectedGtfList,
);

#检测psiOutputDir目录是否存在
if(not -e $psiOutputDir){
	print STDERR "$psiOutputDir doesn't exist!\n";
	exit;
}

#extract experiment Id into array
my $expDirText = "";
$psiOutputDir = $psiOutputDir . "/" if(substr($psiOutputDir, length($psiOutputDir)-1, 1) ne "/");
if($expIdList ne "" and -e $expIdList){
	open FF, "<$expIdList";
	while(my $expId=<FF>){
		chomp($expId);
		if(not -e $psiOutputDir . $expId . "/"){
			print STDERR $psiOutputDir . $expId . " doesn't exist!\n";
		}else{
			$expDir[$#expDir+1] = $psiOutputDir . $expId;
			print STDOUT "Gather experiment information from " . $psiOutputDir . $expId . "/\n";
		}
	}
	close FF;	
}else{
	$expDirText = `find $psiOutputDir -type d -maxdepth 1`;
	@expDir = split(/\n/, $expDirText);
	shift(@expDir);
	for(my $i=0; $i<=$#expDir; $i++){
		print STDOUT "Gather experiment information from " . $expDir[$i] . "/\n";
	}
	
}


#gather information from experiment dir
my ($experimentDir);
my (@tmpArr);
my ($mergeCmdString);
my (%libraryType, %readLength, %layout);
$mergeCmdString = $stringtie . " --merge -G " . $basicRefAnno . " -c " . $mergeCov .
			" -o " . $outputTripletAnno . 
			" -l Trip " .
			$secondRefAnno;

open SELECTEDEXPERIMENT, ">$outputSelectedGtfList";
for(my $i=0; $i<=$#expDir; $i++){

	@tmpArr = ();
	@tmpArr = split(/\//, $expDir[$i]);
	$expId = $tmpArr[$#tmpArr];

	$experimentDir = $expDir[$i] . "/";
	#detect running status
	if(-e $experimentDir . "/" . "transcriptomeByStringtie.gtf"){
		print STDOUT "Detected transcriptome assembled on experiment " . $expId . "\n";
	}else{
		print STDERR "No transcriptome assembled on experiment " . $expId . "\n";
		next;
	}
	
	#detect sequencing information
        if(-e $experimentDir . "/" . $expId . ".SeqInfo.txt"){
                open FF, "<" . $experimentDir . "/" . $expId . ".SeqInfo.txt";
                <FF>;
		%libraryType = ();
		%readLength = ();
		%layout = ();
                while(my $line = <FF>){
                        chomp($line);
                        @tmpArr = ();
                        @tmpArr = split(/\t/, $line);
                        if($#tmpArr==7){
                                $layout{$tmpArr[2]}= 1;
                                $readLength{$tmpArr[3]}= 1;
                                $libraryType{$tmpArr[6]}=1;
                        }elsif($line=~/TotalSpot:(.*)\tAlignPercent:(.*)%/){
                                $totalSpotsInExp = int(($1/1000000)*100)/100;
                                $alignPerInExp = $2;
                        }
                }
                close FF;
		my @libraryArr = keys(%libraryType);
		my @readLengthArr = keys(%readLength);
		my @layoutArr = keys(%layout);
		if(index($libraryType, $libraryArr[0]) < 0 or 
			index($layout, $layoutArr[0]) < 0 or 
			$readLengthArr[0] <$minReadLength or
			$totalSpotsInExp < $totalSpots or
			$alignPerInExp < $alignPer or
			$totalSpotsInExp * $alignPerInExp / 100 < $mappedSpots
		){
			print STDERR "Experiment " . $expId . " can't meet filter demand.\n";
			next;
		}

		$mergeCmdString .= " " . $experimentDir . "/transcriptomeByStringtie.gtf";
		print SELECTEDEXPERIMENT $experimentDir . "/transcriptomeByStringtie.gtf\n";
        }
}
print STDOUT "\n\nCMD: " . $mergeCmdString . "\n";
close SELECTEDEXPERIMENT;
system($mergeCmdString);;
