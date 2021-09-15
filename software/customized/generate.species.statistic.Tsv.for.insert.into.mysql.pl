#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--species \"Capra hircus\"\\\n" .
		"--taxon 9925 \\\n" .
		"--genomeAssemblyName genomeAssemblyName\\\n" .
		"--genomeAssemblyAddress genomeAssemblyAddress.gtf.or.gff3\\\n" .
		"--ensemblAnnotationName release-96 \\\n" .
		"--ensemblAnnotationAddress ensemblAnnotationAddress.gtf.or.gff3 \\\n" .
		"--refseqAnnotationName ARR1 \\\n" .
		"--refseqAnnotationAddress refseqAnnotationAddress.gff3 \\\n" .
		"--jgiAnnotationName jgiAnno \\\n" .
		"--jgiAnnotationAddress jgiAnnotationAddress.gff3 \\\n" .
		"--sampleListForIsoformAnno ../004-combine-assemblies-and-annos/combined.transcriptome.list \\\n". 
		"--sampleLisftForPickPsi ./samples.filtered.by.sequencing.and.mapping.tsv \\\n" .
		"--isoformTsvFile 9925.isoform.tsv \\\n" .
		"--asTsvFile 9925.as.tsv \\\n" .
		"--outputSpeciesTsvFile 9925.speciesTable.tsv\n\n\n";
	exit;
}

my ($species, $taxon, $genomeAssemblyName, $genomeAssemblyAddress, $ensemblAnnotationName, 
	$ensemblAnnotationAddress, $refseqAnnotationName, $refseqAnnotationAddress, $jgiAnnotationName, $isoformTsvFile,
	$jgiAnnotationAddress, $sampleListForIsoformAnno, $sampleLisftForPickPsi, $asTsvFile, $outputSpeciesTsvFile);

GetOptions(
	'species=s'=>\$species,
	'taxon=s'=>\$taxon,
	'genomeAssemblyName=s'=>\$genomeAssemblyName,
	'genomeAssemblyAddress=s'=>\$genomeAssemblyAddress,
	'ensemblAnnotationName=s'=>\$ensemblAnnotationName,
	'ensemblAnnotationAddress=s'=>\$ensemblAnnotationAddress,
	'refseqAnnotationName=s'=>\$refseqAnnotationName,
	'refseqAnnotationAddress=s'=>\$refseqAnnotationAddress,
	'jgiAnnotationName=s'=>\$jgiAnnotationName,
	'jgiAnnotationAddress=s'=>\$jgiAnnotationAddress,
	'sampleListForIsoformAnno=s'=>\$sampleListForIsoformAnno,
	'sampleLisftForPickPsi=s'=>\$sampleLisftForPickPsi,
	'isoformTsvFile=s'=>\$isoformTsvFile,
	'asTsvFile=s'=>\$asTsvFile,
	'outputSpeciesTsvFile=s'=>\$outputSpeciesTsvFile,
);


print &getLocalTime();
print ": Obtaion isoform numbers from ensembl, refseq, jgi and RNAseq and total gene num.\n";
my (%geneId, @geneId, $geneNum, $isoformNum);
my ($ensemblIsoformNum, $refseqIsoformNum, $jgiIsoformNum, $RNAseqIsoformNum);

($geneNum, $isoformNum, $ensemblIsoformNum, $refseqIsoformNum, $jgiIsoformNum, $RNAseqIsoformNum) = (0, 0, 0, 0, 0, 0);

my ($line, @tmp, $fieldNameString, @fieldName, $fieldName, $valueString, @value, $value, %isoform);
open FF, "<$isoformTsvFile";
while($line=<FF>){

        chomp($line);
        @tmp = ();
        @tmp = split(/_____/, $line);
        $fieldNameString = $tmp[0];
        @fieldName = ();
        @fieldName = split(/, /, $fieldNameString);

        $valueString = $tmp[1];
        @value = ();
        @value = split(/, /, $valueString);

        %isoform = ();

        for(my $i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $isoform{$fieldName[$i]}=$value[$i];
        }

	$geneId{$isoform{"geneId"}}=1;
	$isoformNum++;
	if(uc($isoform{"annoOrigin"}) eq "ENSEMBL"){
		$ensemblIsoformNum++
	}elsif(uc($isoform{"annoOrigin"}) eq "REFSEQ"){
		$refseqIsoformNum++;
	}elsif(uc($isoform{"annoOrigin"}) eq "JGI"){
		$jgiIsoformNum++;
	}elsif(uc($isoform{"annoOrigin"}) eq "RNASEQ"){
		$RNAseqIsoformNum++;
	}
}
@geneId = keys(%geneId);
$geneNum = $#geneId + 1;


print &getLocalTime();
print ": Obtain RNAseq experiment ID and number for annotation isoforms.\n";
my ($annoBasedRNAseqExpNum, $annoBasedRNAseqExpList);
$annoBasedRNAseqExpNum=0;
$annoBasedRNAseqExpList="";
open FF, "<$sampleListForIsoformAnno";
#../002-assemble-trsptome-on-goodExps/psiOutputDir/SRX1509344/transcriptomeByStringtie.gtf
while(my $line=<FF>){
	$annoBasedRNAseqExpNum++;
	if($line=~/psiOutputDir\/(.*?)\/transcriptomeByStringtie/){
		$annoBasedRNAseqExpList .= $1 . ",";
	}
}
$annoBasedRNAseqExpList = substr($annoBasedRNAseqExpList, 0, length($annoBasedRNAseqExpList) - 1);
close FF;


print &getLocalTime();
print ": Obtain total RNAseq experiment ID and number for pickup psi.\n";
my ($totalRNAseqExpNum, $totalRNAseqExpList, $totalStudyNum, $totalStudyList);
$totalRNAseqExpNum=0;
$totalRNAseqExpList="";
$totalStudyNum=0;
$totalStudyList="";
my (%studyId, @studyArr);
open FF, "<$sampleLisftForPickPsi";
# SRX1605387      OK      1       SRR3194788      UN      PAIRED  33      100     22.63   96.24%  104777  7099    9067    28091   19291   2301    0
#        0       4385    1213    419     4748    6116    18938   14582   1541    4594    6010    18582   13976   1502
<FF>;
while(my $line=<FF>){
	$totalRNAseqExpNum++;
	if($line=~/^(.*?)\tOK\t\d+\t(.*?)\t/){
		$totalRNAseqExpList .= $1 . ",";
	}
}
$totalRNAseqExpList = substr($totalRNAseqExpList, 0, length($totalRNAseqExpList) - 1);
close FF;




print &getLocalTime();
print ": Obtain num of as from total AS, ByEnsembl, ByEnsembl_Refseq, ByEnsembl_Refseq_RNAseq, ByNovel.\n";

my ($A5ssPercentage, $A3ssPercentage, $SePercentage, $RiPercentage, $MxePercentage) = (0, 0, 0, 0, 0);
my ($asNum, $A5ssNum, $A3ssNum, $SeNum, $RiNum, $MxeNum) = (0, 0, 0, 0, 0, 0);
my ($asNumByEnsembl, $asA5ssNumByEnsembl, $asA3ssNumByEnsembl, $asSeNumByEnsembl, $asRiNumByEnsembl, $asMxeNumByEnsembl) = (0, 0, 0, 0, 0, 0);
my ($asNumByEnsembl_Refseq, $asA5ssNumByEnsembl_Refseq, $asA3ssNumByEnsembl_Refseq, $asSeNumByEnsembl_Refseq, $asRiNumByEnsembl_Refseq, $asMxeNumByEnsembl_Refseq) = (0, 0, 0, 0, 0, 0);
my ($asNumByEnsembl_Refseq_RNAseq, $asA5ssNumByEnsembl_Refseq_RNAseq, $asA3ssNumByEnsembl_Refseq_RNAseq, $asSeNumByEnsembl_Refseq_RNAseq, $asRiNumByEnsembl_Refseq_RNAseq, $asMxeNumByEnsembl_Refseq_RNAseq) = (0, 0, 0, 0, 0, 0);
my ($asNumByNovel, $asA5ssNumByNovel, $asA3ssNumByNovel, $asSeNumByNovel, $asRiNumByNovel, $asMxeNumByNovel) = (0, 0, 0, 0, 0, 0);

my (%as);
open FF, "<$asTsvFile";
while($line=<FF>){
        chomp($line);
        @tmp = ();
        @tmp = split(/_____/, $line);
        $fieldNameString = $tmp[0];
        @fieldName = ();
        @fieldName = split(/, /, $fieldNameString);

        $valueString = $tmp[1];
        @value = ();
        @value = split(/, /, $valueString);

        %as = ();

        for(my $i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $as{$fieldName[$i]}=$value[$i];
        }

	if(uc($as{"asType"}) eq "A5SS"){
		$A5ssNum++;
		if(uc($as{"discoveryApproach"}) eq "ENSEMBL"){
			$asA5ssNumByEnsembl++;
		}elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ"){
			$asA5ssNumByEnsembl_Refseq++; 
		}elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ_RNASEQ"){
			$asA5ssNumByEnsembl_Refseq_RNAseq++;
		}elsif(uc($as{"discoveryApproach"}) eq "NOVEL"){
			$asA5ssNumByNovel++;
		}
	}elsif(uc($as{"asType"}) eq "A3SS"){
                $A3ssNum++;
                if(uc($as{"discoveryApproach"}) eq "ENSEMBL"){
                        $asA3ssNumByEnsembl++;
                }elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ"){
                        $asA3ssNumByEnsembl_Refseq++;
                }elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ_RNASEQ"){
                        $asA3ssNumByEnsembl_Refseq_RNAseq++;
                }elsif(uc($as{"discoveryApproach"}) eq "NOVEL"){
                        $asA3ssNumByNovel++;
                }
	}elsif(uc($as{"asType"}) eq "SE"){
		$SeNum++;
                if(uc($as{"discoveryApproach"}) eq "ENSEMBL"){
                        $asSeNumByEnsembl++;
                }elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ"){
                        $asSeNumByEnsembl_Refseq++;
                }elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ_RNASEQ"){
                        $asSeNumByEnsembl_Refseq_RNAseq++;
                }elsif(uc($as{"discoveryApproach"}) eq "NOVEL"){
                        $asSeNumByNovel++;
                }
	}elsif(uc($as{"asType"}) eq "RI"){
		$RiNum++;
                if(uc($as{"discoveryApproach"}) eq "ENSEMBL"){
                        $asRiNumByEnsembl++;
                }elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ"){
                        $asRiNumByEnsembl_Refseq++;
                }elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ_RNASEQ"){
                        $asRiNumByEnsembl_Refseq_RNAseq++;
                }elsif(uc($as{"discoveryApproach"}) eq "NOVEL"){
                        $asRiNumByNovel++;
                }
	}elsif(uc($as{"asType"}) eq "MXE"){
		$MxeNum++;
                if(uc($as{"discoveryApproach"}) eq "ENSEMBL"){
                        $asMxeNumByEnsembl++;
                }elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ"){
                        $asMxeNumByEnsembl_Refseq++;
                }elsif(uc($as{"discoveryApproach"}) eq "ENSEMBL_REFSEQ_RNASEQ"){
                        $asMxeNumByEnsembl_Refseq_RNAseq++;
                }elsif(uc($as{"discoveryApproach"}) eq "NOVEL"){
                        $asMxeNumByNovel++;
                }
	}
}
close FF;

$asNum = $A5ssNum+$A3ssNum+$SeNum+$RiNum+$MxeNum;
$A5ssPercentage = sprintf("%.2f",$A5ssNum/$asNum);
$A3ssPercentage = sprintf("%.2f", $A3ssNum/$asNum);
$SePercentage = sprintf("%.2f", $SeNum/$asNum);
$RiPercentage = sprintf("%.2f", $RiNum/$asNum);
$MxePercentage = sprintf("%.2f", $MxeNum/$asNum);

my $isoformAnnoFileName = "";
my $totalStudyNum = 0;
my $totalStudyList = "";
open WW, ">$outputSpeciesTsvFile";
print WW "species, taxon, genomeAssemblyName, genomeAssemblyAddress, geneNum, isoformNum, isoformAnnoFileName, ensemblAnnotationName, ensemblAnnotationAddress, refseqAnnotationName, refseqAnnotationAddress, jgiAnnotationName, jgiAnnotationAddress, annoBasedRNAseqExpNum, annoBasedRNAseqExpList, ensemblIsoformNum, refseqIsoformNum, jgiIsoformNum, RNAseqIsoformNum, totalRNAseqExpNum, totalRNAseqExpList, totalStudyNum, totalStudyList, asNum, A5ssPercentage, A3ssPercentage, SePercentage, RiPercentage, MxePercentage, asNumByEnsembl, asA5ssNumByEnsembl, asA3ssNumByEnsembl, asSeNumByEnsembl, asRiNumByEnsembl, asMxeNumByEnsembl, asNumByEnsembl_Refseq, asA5ssNumByEnsembl_Refseq, asA3ssNumByEnsembl_Refseq, asSeNumByEnsembl_Refseq, asRiNumByEnsembl_Refseq, asMxeNumByEnsembl_Refseq, asNumByEnsembl_Refseq_RNAseq, asA5ssNumByEnsembl_Refseq_RNAseq, asA3ssNumByEnsembl_Refseq_RNAseq, asSeNumByEnsembl_Refseq_RNAseq, asRiNumByEnsembl_Refseq_RNAseq, asMxeNumByEnsembl_Refseq_RNAseq, asNumByNovel, asA5ssNumByNovel, asA3ssNumByNovel, asSeNumByNovel, asRiNumByNovel, asMxeNumByNovel" . "_____" . "\"". $species . "\", \"" . $taxon . "\", \"" . $genomeAssemblyName . "\", \"" . $genomeAssemblyAddress . "\", " . $geneNum . ", " . $isoformNum . ", \"" . $isoformAnnoFileName . "\", \"" . $ensemblAnnotationName . "\", \"" . $ensemblAnnotationAddress . "\", \"" . $refseqAnnotationName . "\", \"" . $refseqAnnotationAddress . "\", \"" . $jgiAnnotationName . "\", \"" . $jgiAnnotationAddress . "\", " . $annoBasedRNAseqExpNum . ", \"" . $annoBasedRNAseqExpList . "\", " . $ensemblIsoformNum . ", " . $refseqIsoformNum . ", " . $jgiIsoformNum . ", " . $RNAseqIsoformNum . ", " . $totalRNAseqExpNum . ", \"" . $totalRNAseqExpList . "\", " . $totalStudyNum . ", \"" . $totalStudyList . "\", " . $asNum . ", " . $A5ssPercentage . ", " . $A3ssPercentage . ", " . $SePercentage . ", " . $RiPercentage . ", " . $MxePercentage . ", " . $asNumByEnsembl . ", " . $asA5ssNumByEnsembl . ", " . $asA3ssNumByEnsembl . ", " . $asSeNumByEnsembl . ", " . $asRiNumByEnsembl . ", " . $asMxeNumByEnsembl . ", " . $asNumByEnsembl_Refseq . ", " . $asA5ssNumByEnsembl_Refseq . ", " . $asA3ssNumByEnsembl_Refseq . ", " . $asSeNumByEnsembl_Refseq . ", " . $asRiNumByEnsembl_Refseq . ", " . $asMxeNumByEnsembl_Refseq . ", " . $asNumByEnsembl_Refseq_RNAseq . ", " . $asA5ssNumByEnsembl_Refseq_RNAseq . ", " . $asA3ssNumByEnsembl_Refseq_RNAseq . ", " . $asSeNumByEnsembl_Refseq_RNAseq . ", " . $asRiNumByEnsembl_Refseq_RNAseq . ", " . $asMxeNumByEnsembl_Refseq_RNAseq . ", " . $asNumByNovel . ", " . $asA5ssNumByNovel . ", " . $asA3ssNumByNovel . ", " . $asSeNumByNovel . ", " . $asRiNumByNovel . ", " . $asMxeNumByNovel
;
close WW;

sub getLocalTime{
        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
        return sprintf("%02d:%02d:%02d %02d-%02d-%2d", $hour, $min, $sec, $mon, $mday, $year);
}

