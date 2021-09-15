#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD njaulsas2019 \\\n" .
		"--species \"Ovis aries\"\\\n" .
		"--taxon 9940 \\\n" .
		"--asTsvFile 9940.as.tsv \\\n" .
		"--speciesTsvFile 9940.species.tsv \\\n" .
		"--psiTsvFile 9940.psi.tsv \\\n" .
		"--variantTsvFile 9940.variant.tsv \\\n" .
		"--isoformTsvFile 9940.isoform.tsv \\\n" .
		"--experimentTsvFile 9940.experiment.tsv \n\n";
	exit;
}

my ($dbName, $dbUser, $dbPWD, $species, $taxon);
my ($asTsvFile, $speciesTsvFile, $psiTsvFile, $variantTsvFile, $isoformTsvFile, $experimentTsvFile);
GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'species=s'=>\$species,
	'taxon=s'=>\$taxon,
	'asTsvFile=s'=>\$asTsvFile,
	'speciesTsvFile=s'=>\$speciesTsvFile,
	'psiTsvFile=s'=>\$psiTsvFile,
	'variantTsvFile=s'=>\$variantTsvFile,
	'isoformTsvFile=s'=>\$isoformTsvFile,
	'experimentTsvFile=s'=>\$experimentTsvFile,
);

my $dbh = DBI->connect("DBI:mysql:database=$dbName;", $dbUser, $dbPWD);
$dbh->{mysql_auto_reconnect} = 1;

my ($line, $fieldNameString, @fieldName, $fieldName, $valueString, @value, $value);
my (@tmp, $i, %record, $sql, $insert, $delete);


####################### asTable ##############################

print &getLocalTime();
print ": Begin delete alternative splicing from mysql table.\n";
$sql = "delete from asTable where species=\"" . $species . "\"";
$delete = $dbh->prepare($sql);
$delete->execute();
print &getLocalTime();
print ": Finish delete alternative splicing from mysql table.\n";

print &getLocalTime();
print ": Begin load alternative splicing into mysql table.\n";
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

	%record = ();
        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $record{$fieldName[$i]}=$value[$i];
        }

	$sql = "insert into asTable (asId, asType, species, chr, strand, start, end, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base, downstreamEE, downstreamES, exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE, shortES, upstreamEE, upstreamES, geneID, geneSymbol, firstAltExonSeries, firstIsoformIdList, firstIsoformNum, secondAltExonSeries, secondIsoformIdList, secondIsoformNum, pfam, go, variantPointNum, variantPointTypeCmb, asOrthId, conservedSpeciesNum, jcecExperimentNum, jcExperimentNum, discoveryApproach ) values (\"" . $record{"asId"} . "\", \"" . $record{"asType"} . "\", \"" . $record{"species"} . "\", \"" . $record{"chr"} . "\", \"" . $record{"strand"} . "\", " . $record{"start"} . ", " . $record{"end"} . ", " . $record{"1stExonEnd"} . ", " . $record{"1stExonStart_0base"} . ", " . $record{"2ndExonEnd"} . ", " . $record{"2ndExonStart_0base"} . ", " . $record{"downstreamEE"} . ", " . $record{"downstreamES"} . ", " . $record{"exonEnd"} . ", " . $record{"exonStart_0base"} . ", " . $record{"flankingEE"} . ", " . $record{"flankingES"} . ", " . $record{"longExonEnd"} . ", " . $record{"longExonStart_0base"} . ", " . $record{"riExonEnd"} . ", " . $record{"riExonStart_0base"} . ", " . $record{"shortEE"} . ", " . $record{"shortES"} . ", " . $record{"upstreamEE"} . ", " . $record{"upstreamES"} . ", \"" . $record{"geneID"} . "\", \"" . $record{"geneSymbol"} . "\", \"" . $record{"firstAltExonSeries"} . "\", \"" . $record{"firstIsoformIdList"} . "\", " . $record{"firstIsoformNum"} . ", \"" . $record{"secondAltExonSeries"} . "\", \"" . $record{"secondIsoformIdList"} . "\", " . $record{"secondIsoformNum"} . ", \"" . $record{"pfam"} . "\", \"" . $record{"go"} . "\", " . $record{"variantPointNum"} . ", \"" . $record{"variantPointTypeCmb"} . "\", \"" . $record{"asOrthId"} . "\", " . $record{"conservedSpeciesNum"} . ", " . $record{"jcecExperimentNum"} . ", " . $record{"jcExperimentNum"} . ", \"" . $record{"discoveryApproach"} . "\")";
	
	$insert = $dbh->prepare($sql);
	$insert->execute();

}
close FF;

print &getLocalTime();
print ": Finish load alternative splicing into mysql table.\n";


########################## isoformTable ###############################
print &getLocalTime();
print ": Begin delete isoforms from mysql table.\n";
$sql = "delete from isoformTable where species=\"" . $species . "\"";
$delete = $dbh->prepare($sql);
$delete->execute();
print &getLocalTime();
print ": Finish delete isoforms from mysql table.\n";

print &getLocalTime();
print ": Begin load isoform into mysql table.\n";
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

        %record = ();
        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $record{$fieldName[$i]}=$value[$i];
        }

	$sql = "insert into isoformTable (isoformId, species, isoformSymbol, geneId, geneSymbol, chr, strand, start, end, exonSeries, annoOrigin, cDNAseq, DNAseq, pepSeq, isoformPfam, isoformGo, proteinId) values (\"" . $record{"isoformId"} . "\", \"" . $record{"species"} . "\", \"" . $record{"isoformSymbol"} . "\", \"" . $record{"geneId"} . "\", \"" . $record{"geneSymbol"} . "\", \"" . $record{"chr"} . "\", \"" . $record{"strand"} . "\", " . $record{"start"} . ", " . $record{"end"} . ", \"" . $record{"exonSeries"} . "\", \"" . $record{"annoOrigin"} . "\", \"" . $record{"cDNAseq"} . "\", \"" . $record{"DNAseq"} . "\", \"" . $record{"pepSeq"} . "\", \"" . $record{"isoformPfam"} . "\", \"" . $record{"isoformGo"} . "\", \"" . $record{"proteinId"} . "\")";;


        $insert = $dbh->prepare($sql);
        $insert->execute();
}
close FF;

print &getLocalTime();
print ": Finish load isoform into mysql table.\n";



############################# experimentTable ################################
print &getLocalTime();
print ": Begin delete experiments from mysql table.\n";
$sql = "delete from experimentTable where species=\"" . $species . "\"";
$delete = $dbh->prepare($sql);
$delete->execute();
print &getLocalTime();
print ": Finish delete experiments from mysql table.\n";

print &getLocalTime();
print ": Begin load experiments into mysql table.\n";
open FF, "<$experimentTsvFile";
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

        %record = ();
        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $record{$fieldName[$i]}=$value[$i];
        }

	$sql = "insert into experimentTable(alignPercent, mappedSpots, experimentId, species, libraryType, libraryLayout, readLen, phredScore, totalSpots, runIdList, runNum, studyId, jcecTotalAsNum, jcecA5ssPercentage, jcecA3ssPercentage, jcecSePercentage, jcecRiPercentage, jcecMxePercentage, jcTotalAsNum, jcA5ssPercentage, jcA3ssPercentage, jcSePercentage, jcRiPercentage, jcMxePercentage) values (" . $record{"alignPercent"} . ", " . $record{"mappedSpots"} . ", \"" . $record{"experimentId"} . "\", \"" . $record{"species"} . "\", \"" . $record{"libraryType"} . "\", \"" . $record{"libraryLayout"} . "\", \"" . $record{"readLen"} . "\", \"" . $record{"phredScore"} . "\", " . $record{"totalSpots"} . ", \"" . $record{"runIdList"} . "\", " . $record{"runNum"} . ", \"" . $record{"studyId"} . "\", " . $record{"jcecTotalAsNum"} . ", " . $record{"jcecA5ssPercentage"} . ", " . $record{"jcecA3ssPercentage"} . ", " . $record{"jcecSePercentage"} . ", " . $record{"jcecRiPercentage"} . ", " . $record{"jcecMxePercentage"} . ", " . $record{"jcTotalAsNum"} . ", " . $record{"jcA5ssPercentage"} . ", " . $record{"jcA3ssPercentage"} . ", " . $record{"jcSePercentage"} . ", " . $record{"jcRiPercentage"} . ", " . $record{"jcMxePercentage"} . ")";


        $insert = $dbh->prepare($sql);
        $insert->execute();
}
close FF;

print &getLocalTime();
print ": Finish load experiments into mysql table.\n";


###########################  psiTable ##############################
print &getLocalTime();
print ": Begin delete psi from mysql table.\n";
$sql = "delete from psiTable where species=\"" . $species . "\"";
$delete = $dbh->prepare($sql);
$delete->execute();
print &getLocalTime();
print ": Finish delete psi from mysql table.\n";

print &getLocalTime();
print ": Begin load psi into mysql table.\n";
open FF, "<$psiTsvFile";
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

        %record = ();
        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $record{$fieldName[$i]}=$value[$i];
        }

	$sql = "insert into psiTable (asId, species, experiment, studyId, JCECI, JCECS, JCECIncFormLen, JCECSkipFormLen, JCECpsi, JCI, JCS, JCIncFormLen, JCSkipFormLen, JCpsi ) values(\"" . $record{"asId"} . "\", \"" . $record{"species"} . "\", \"" . $record{"experiment"} . "\", \"" . $record{"studyId"} . "\", " . $record{"JCECI"} . ", " . $record{"JCECS"} . ", " . $record{"JCECIncFormLen"} . ", " . $record{"JCECSkipFormLen"} . ", " . $record{"JCECpsi"} . ", " . $record{"JCI"} . ", " . $record{"JCS"} . ", " . $record{"JCIncFormLen"} . ", " . $record{"JCSkipFormLen"} . ", " . $record{"JCpsi"} . ")";

        $insert = $dbh->prepare($sql);
        $insert->execute();
}
close FF;

print &getLocalTime();
print ": Finish load psi into mysql table.\n";



#######################  variantTable #######################
print &getLocalTime();
print ": Begin delete variants from mysql table.\n";
$sql = "delete from variantTable where species=\"" . $species . "\"";
$delete = $dbh->prepare($sql);
$delete->execute();
print &getLocalTime();
print ": Finish delete variants from mysql table.\n";

print &getLocalTime();
print ": Begin load variant into mysql table.\n";
open FF, "<$variantTsvFile";
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

        %record = ();
        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $record{$fieldName[$i]}=$value[$i];
        }

	$sql ="insert into variantTable ( species, chr, pos, ref, alt, accId, varType ) values (\"" . $record{"species"} . "\", \"" . $record{"chr"} . "\", " . $record{"pos"} . ", \"" . $record{"ref"} . "\", \"" . $record{"alt"} . "\", \"" . $record{"accId"} . "\", \"" . $record{"varType"} . "\")";

        $insert = $dbh->prepare($sql);
        $insert->execute();
}
close FF;

print &getLocalTime();
print ": Finish load psi into mysql table.\n";



############################  speciesTable #########################
print &getLocalTime();
print ": Begin delete species from mysql table.\n";
$sql = "delete from speciesTable where species=\"" . $species . "\"";
$delete = $dbh->prepare($sql);
$delete->execute();
print &getLocalTime();
print ": Finish delete species from mysql table.\n";

print &getLocalTime();
print ": Begin load species into mysql table.\n";
open FF, "<$speciesTsvFile";
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

        %record = ();
        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $record{$fieldName[$i]}=$value[$i];
        }

	$sql = "insert into speciesTable (species, taxon, genomeAssemblyName, genomeAssemblyAddress, geneNum, isoformNum, isoformAnnoFileName, ensemblAnnotationName, ensemblAnnotationAddress, refseqAnnotationName, refseqAnnotationAddress, jgiAnnotationName, jgiAnnotationAddress, annoBasedRNAseqExpNum, annoBasedRNAseqExpList, ensemblIsoformNum, refseqIsoformNum, jgiIsoformNum, RNAseqIsoformNum, totalRNAseqExpNum, totalRNAseqExpList, totalStudyNum, totalStudyList, asNum, A5ssPercentage, A3ssPercentage, SePercentage, RiPercentage, MxePercentage, asNumByEnsembl, asA5ssNumByEnsembl, asA3ssNumByEnsembl, asSeNumByEnsembl, asRiNumByEnsembl, asMxeNumByEnsembl, asNumByEnsembl_Refseq, asA5ssNumByEnsembl_Refseq, asA3ssNumByEnsembl_Refseq, asSeNumByEnsembl_Refseq, asRiNumByEnsembl_Refseq, asMxeNumByEnsembl_Refseq, asNumByEnsembl_Refseq_RNAseq, asA5ssNumByEnsembl_Refseq_RNAseq, asA3ssNumByEnsembl_Refseq_RNAseq, asSeNumByEnsembl_Refseq_RNAseq, asRiNumByEnsembl_Refseq_RNAseq, asMxeNumByEnsembl_Refseq_RNAseq, asNumByNovel, asA5ssNumByNovel, asA3ssNumByNovel, asSeNumByNovel, asRiNumByNovel, asMxeNumByNovel) values(\"" . $record{"species"} . "\", \"" . $record{"taxon"} . "\", \"" . $record{"genomeAssemblyName"} . "\", \"" . $record{"genomeAssemblyAddress"} . "\", " . $record{"geneNum"} . ", " . $record{"isoformNum"} . ", \"" . $record{"isoformAnnoFileName"} . "\", \"" . $record{"ensemblAnnotationName"} . "\", \"" . $record{"ensemblAnnotationAddress"} . "\", \"" . $record{"refseqAnnotationName"} . "\", \"" . $record{"refseqAnnotationAddress"} . "\", \"" . $record{"jgiAnnotationName"} . "\", \"" . $record{"jgiAnnotationAddress"} . "\", " . $record{"annoBasedRNAseqExpNum"} . ", \"" . $record{"annoBasedRNAseqExpList"} . "\", " . $record{"ensemblIsoformNum"} . ", " . $record{"refseqIsoformNum"} . ", " . $record{"jgiIsoformNum"} . ", " . $record{"RNAseqIsoformNum"} . ", " . $record{"totalRNAseqExpNum"} . ", \"" . $record{"totalRNAseqExpList"} . "\", " . $record{"totalStudyNum"} . ", \"" . $record{"totalStudyList"} . "\", " . $record{"asNum"} . ", " . $record{"A5ssPercentage"} . ", " . $record{"A3ssPercentage"} . ", " . $record{"SePercentage"} . ", " . $record{"RiPercentage"} . ", " . $record{"MxePercentage"} . ", " . $record{"asNumByEnsembl"} . ", " . $record{"asA5ssNumByEnsembl"} . ", " . $record{"asA3ssNumByEnsembl"} . ", " . $record{"asSeNumByEnsembl"} . ", " . $record{"asRiNumByEnsembl"} . ", " . $record{"asMxeNumByEnsembl"} . ", " . $record{"asNumByEnsembl_Refseq"} . ", " . $record{"asA5ssNumByEnsembl_Refseq"} . ", " . $record{"asA3ssNumByEnsembl_Refseq"} . ", " . $record{"asSeNumByEnsembl_Refseq"} . ", " . $record{"asRiNumByEnsembl_Refseq"} . ", " . $record{"asMxeNumByEnsembl_Refseq"} . ", " . $record{"asNumByEnsembl_Refseq_RNAseq"} . ", " . $record{"asA5ssNumByEnsembl_Refseq_RNAseq"} . ", " . $record{"asA3ssNumByEnsembl_Refseq_RNAseq"} . ", " . $record{"asSeNumByEnsembl_Refseq_RNAseq"} . ", " . $record{"asRiNumByEnsembl_Refseq_RNAseq"} . ", " . $record{"asMxeNumByEnsembl_Refseq_RNAseq"} . ", " . $record{"asNumByNovel"} . ", " . $record{"asA5ssNumByNovel"} . ", " . $record{"asA3ssNumByNovel"} . ", " . $record{"asSeNumByNovel"} . ", " . $record{"asRiNumByNovel"} . ", " . $record{"asMxeNumByNovel"} . ")";

        $insert = $dbh->prepare($sql);
        $insert->execute();
}
close FF;

print &getLocalTime();
print ": Finish load psi into mysql table.\n";




$dbh->disconnect();
sub getLocalTime{
        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
        return sprintf("%02d:%02d:%02d %02d-%02d-%2d", $hour, $min, $sec, $mon, $mday, $year);
}

