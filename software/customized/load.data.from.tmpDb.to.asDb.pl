#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--tmpDbName tmpdb \\\n" .
                "--tmpDbUser lsas \\\n" .
                "--tmpDbPWD lsasnjau2019 \\\n" .
                "--asDbName asdb \\\n" .
                "--asDbUser lsas \\\n" .
                "--asDbPWD lsasnjau2019 \\\n" .
		"--species \"Equus caballus\" \n";
	exit;
}

my ($tmpDbName, $tmpDbUser, $tmpDbPWD, $asDbName, $asDbUser, $asDbPWD, $species);
GetOptions(
        'tmpDbName=s'=>\$tmpDbName,
        'tmpDbUser=s'=>\$tmpDbUser,
        'tmpDbPWD=s'=>\$tmpDbPWD,
        'asDbName=s'=>\$asDbName,
        'asDbUser=s'=>\$asDbUser,
        'asDbPWD=s'=>\$asDbPWD,
	'species=s'=>\$species,
);

my ($dbhTmpDb, $dbhAsDb, $sqlTmpDb, $queryTmpDb, $rowTmpDb, $sqlAsDb, $insertAsDb, $deleteAsDb, $deleteTmpDb);
$dbhTmpDb = DBI->connect("DBI:mysql:database=$tmpDbName;", $tmpDbUser, $tmpDbPWD);
$dbhTmpDb->{mysql_auto_reconnect} = 1;
$dbhAsDb = DBI->connect("DBI:mysql:database=$asDbName;", $asDbUser, $asDbPWD);
$dbhAsDb->{mysql_auto_reconnect} = 1;

$sqlTmpDb = "delete from speciesTable where species!=\"" . $species . "\"";
$deleteTmpDb = $dbhTmpDb->prepare($sqlTmpDb);
$deleteTmpDb->execute();
open WW, ">./sql.txt";
# input species into asDB
$sqlTmpDb = "select * from speciesTable where species=\"$species\"";
$queryTmpDb = $dbhTmpDb->prepare($sqlTmpDb);
$queryTmpDb->execute();
while($rowTmpDb=$queryTmpDb->fetchrow_hashref()){

	$sqlAsDb = "delete from speciesTable where species=\"" . $rowTmpDb->{"species"} . "\"";
	$deleteAsDb = $dbhAsDb->prepare($sqlAsDb);
	$deleteAsDb->execute();
	print "deleted $species from speciesTable @ asdb\n";
	<STDIN>;
	$sqlAsDb = "insert into speciesTable ( species, taxon, genomeAssemblyName, genomeAssemblyAddress, geneNum, isoformNum, isoformAnnoFileName, ensemblAnnotationName, ensemblAnnotationAddress, refseqAnnotationName, refseqAnnotationAddress, jgiAnnotationName, jgiAnnotationAddress, annoBasedRNAseqExpNum, annoBasedRNAseqExpList, ensemblIsoformNum, refseqIsoformNum, jgiIsoformNum, RNAseqIsoformNum, totalRNAseqExpNum, asNum, A5ssPercentage, A3ssPercentage, SePercentage, RiPercentage, MxePercentage, asNumByEnsembl, asA5ssNumByEnsembl, asA3ssNumByEnsembl, asSeNumByEnsembl, asRiNumByEnsembl, asMxeNumByEnsembl, asNumByEnsembl_Refseq, asA5ssNumByEnsembl_Refseq, asA3ssNumByEnsembl_Refseq, asSeNumByEnsembl_Refseq, asRiNumByEnsembl_Refseq, asMxeNumByEnsembl_Refseq, asNumByEnsembl_Refseq_RNAseq, asA5ssNumByEnsembl_Refseq_RNAseq, asA3ssNumByEnsembl_Refseq_RNAseq, asSeNumByEnsembl_Refseq_RNAseq, asRiNumByEnsembl_Refseq_RNAseq, asMxeNumByEnsembl_Refseq_RNAseq, asNumByNovel, asA5ssNumByNovel, asA3ssNumByNovel, asSeNumByNovel, asRiNumByNovel, asMxeNumByNovel) values( \"" . $rowTmpDb->{"species"} . "\", \"" . $rowTmpDb->{"taxon"} . "\", \"" . $rowTmpDb->{"genomeAssemblyName"} . "\", \"" . $rowTmpDb->{"genomeAssemblyAddress"} . "\", " . $rowTmpDb->{"geneNum"} . ", " . $rowTmpDb->{"isoformNum"} . ", \"" . $rowTmpDb->{"isoformAnnoFileName"} . "\", \"" . $rowTmpDb->{"ensemblAnnotationName"} . "\", \"" . $rowTmpDb->{"ensemblAnnotationAddress"} . "\", \"" . $rowTmpDb->{"refseqAnnotationName"} . "\", \"" . $rowTmpDb->{"refseqAnnotationAddress"} . "\", \"" . $rowTmpDb->{"jgiAnnotationName"} . "\", \"" . $rowTmpDb->{"jgiAnnotationAddress"} . "\", " . $rowTmpDb->{"annoBasedRNAseqExpNum"} . ", \"" . $rowTmpDb->{"annoBasedRNAseqExpList"} . "\", " . $rowTmpDb->{"ensemblIsoformNum"} . ", " . $rowTmpDb->{"refseqIsoformNum"} . ", " . $rowTmpDb->{"jgiIsoformNum"} . ", " . $rowTmpDb->{"RNAseqIsoformNum"} . ", " . $rowTmpDb->{"totalRNAseqExpNum"} . ", " . $rowTmpDb->{"asNum"} . ", " . $rowTmpDb->{"A5ssPercentage"} . ", " . $rowTmpDb->{"A3ssPercentage"} . ", " . $rowTmpDb->{"SePercentage"} . ", " . $rowTmpDb->{"RiPercentage"} . ", " . $rowTmpDb->{"MxePercentage"} . ", " . $rowTmpDb->{"asNumByEnsembl"} . ", " . $rowTmpDb->{"asA5ssNumByEnsembl"} . ", " . $rowTmpDb->{"asA3ssNumByEnsembl"} . ", " . $rowTmpDb->{"asSeNumByEnsembl"} . ", " . $rowTmpDb->{"asRiNumByEnsembl"} . ", " . $rowTmpDb->{"asMxeNumByEnsembl"} . ", " . $rowTmpDb->{"asNumByEnsembl_Refseq"} . ", " . $rowTmpDb->{"asA5ssNumByEnsembl_Refseq"} . ", " . $rowTmpDb->{"asA3ssNumByEnsembl_Refseq"} . ", " . $rowTmpDb->{"asSeNumByEnsembl_Refseq"} . ", " . $rowTmpDb->{"asRiNumByEnsembl_Refseq"} . ", " . $rowTmpDb->{"asMxeNumByEnsembl_Refseq"} . ", " . $rowTmpDb->{"asNumByEnsembl_Refseq_RNAseq"} . ", " . $rowTmpDb->{"asA5ssNumByEnsembl_Refseq_RNAseq"} . ", " . $rowTmpDb->{"asA3ssNumByEnsembl_Refseq_RNAseq"} . ", " . $rowTmpDb->{"asSeNumByEnsembl_Refseq_RNAseq"} . ", " . $rowTmpDb->{"asRiNumByEnsembl_Refseq_RNAseq"} . ", " . $rowTmpDb->{"asMxeNumByEnsembl_Refseq_RNAseq"} . ", " . $rowTmpDb->{"asNumByNovel"} . ", " . $rowTmpDb->{"asA5ssNumByNovel"} . ", " . $rowTmpDb->{"asA3ssNumByNovel"} . ", " . $rowTmpDb->{"asSeNumByNovel"} . ", " . $rowTmpDb->{"asRiNumByNovel"} . ", " . $rowTmpDb->{"asMxeNumByNovel"} . "\")";
	print WW $sqlAsDb . "\n";	
	$insertAsDb = $dbhAsDb->prepare($sqlAsDb);
	$insertAsDb->execute();
	print "has insert $species into speciesTable" . "\n";
	<STDIN>;
}
$queryTmpDb->finish();
close WW;
# input as into asDB
$sqlTmpDb = "select * from asTable";
$queryTmpDb = $dbhTmpDb->prepare($sqlTmpDb);
$queryTmpDb->execute();
while($rowTmpDb=$queryTmpDb->fetchrow_hashref()){

	$sqlAsDb = "insert into asTable ( asId, asType, species, chr, strand, start, end, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base, downstreamEE, downstreamES, exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE, shortES, upstreamEE, upstreamES, geneID, geneSymbol, firstAltExonSeries, firstIsoformIdList, firstIsoformNum, secondAltExonSeries, secondIsoformIdList, secondIsoformNum, pfam, go, variantPointNum, variantPointTypeCmb, asOrthId, conservedSpeciesNum, jcecExperimentNum, jcExperimentNum, discoveryApproach ) values(\"" . $rowTmpDb->{"sId"} . "\", \"" . $rowTmpDb->{"asType"} . "\", \"" . $rowTmpDb->{"species"} . "\", \"" . $rowTmpDb->{"chr"} . "\", \"" . $rowTmpDb->{"strand"} . "\", " . $rowTmpDb->{"start"} . ", " . $rowTmpDb->{"end"} . ", " . $rowTmpDb->{"1stExonEnd"} . ", " . $rowTmpDb->{"1stExonStart_0base"} . ", " . $rowTmpDb->{"2ndExonEnd"} . ", " . $rowTmpDb->{"2ndExonStart_0base"} . ", " . $rowTmpDb->{"downstreamEE"} . ", " . $rowTmpDb->{"downstreamES"} . ", " . $rowTmpDb->{"exonEnd"} . ", " . $rowTmpDb->{"exonStart_0base"} . ", " . $rowTmpDb->{"flankingEE"} . ", " . $rowTmpDb->{"flankingES"} . ", " . $rowTmpDb->{"longExonEnd"} . ", " . $rowTmpDb->{"longExonStart_0base"} . ", " . $rowTmpDb->{"riExonEnd"} . ", " . $rowTmpDb->{"riExonStart_0base"} . ", " . $rowTmpDb->{"shortEE"} . ", " . $rowTmpDb->{"shortES"} . ", " . $rowTmpDb->{"upstreamEE"} . ", " . $rowTmpDb->{"upstreamES"} . ", \"" . $rowTmpDb->{"geneID"} . "\", \"" . $rowTmpDb->{"geneSymbol"} . "\", \"" . $rowTmpDb->{"firstAltExonSeries"} . "\", \"" . $rowTmpDb->{"firstIsoformIdList"} . "\", " . $rowTmpDb->{"firstIsoformNum"} . ", \"" . $rowTmpDb->{"secondAltExonSeries"} . "\", \"" . $rowTmpDb->{"secondIsoformIdList"} . "\", " . $rowTmpDb->{"secondIsoformNum"} . ", \"" . $rowTmpDb->{"pfam"} . "\", \"" . $rowTmpDb->{"go"} . "\", " . $rowTmpDb->{"variantPointNum"} . ", \"" . $rowTmpDb->{"variantPointTypeCmb"} . "\", \"" . $rowTmpDb->{"asOrthId"} . "\", " . $rowTmpDb->{"conservedSpeciesNum"} . ", " . $rowTmpDb->{"jcecExperimentNum"} . ", " . $rowTmpDb->{"jcExperimentNum"} . ", \"" . $rowTmpDb->{"discoveryApproach"} . "\")";
	$insertAsDb = $dbhAsDb->prepare($sqlAsDb);
	$insertAsDb->execute();

	print $sqlAsDb . "\n";
	<STDIN>;
}
$queryTmpDb->finish();

# input experiment into asDB
$sqlTmpDb = "select * from experimentTable";
$queryTmpDb = $dbhTmpDb->prepare($sqlTmpDb);
$queryTmpDb->execute();
while($rowTmpDb=$queryTmpDb->fetchrow_hashref()){

	$sqlAsDb = "insert into experimentTable ( experimentId, species, tissue, libraryType, libraryLayout, readLen, phredScore, totalSpots, alignPercent, mappedSpots, runIdList, runNum, studyId, jcecTotalAsNum, jcecA5ssPercentage, jcecA3ssPercentage, jcecSePercentage, jcecRiPercentage, jcecMxePercentage, jcTotalAsNum, jcA5ssPercentage, jcA3ssPercentage, jcSePercentage, jcRiPercentage, jcMxePercentage ) values(\"" . $rowTmpDb->{"experimentId"} . "\", \"" . $rowTmpDb->{"species"} . "\", \"" . $rowTmpDb->{"tissue"} . "\", \"" . $rowTmpDb->{"libraryType"} . "\", \"" . $rowTmpDb->{"libraryLayout"} . "\", \"" . $rowTmpDb->{"readLen"} . "\", \"" . $rowTmpDb->{"phredScore"} . "\", " . $rowTmpDb->{"totalSpots"} . ", " . $rowTmpDb->{"alignPercent"} . ", " . $rowTmpDb->{"mappedSpots"} . ", \"" . $rowTmpDb->{"runIdList"} . "\", " . $rowTmpDb->{"runNum"} . ", \"" . $rowTmpDb->{"studyId"} . "\", " . $rowTmpDb->{"jcecTotalAsNum"} . ", " . $rowTmpDb->{"jcecA5ssPercentage"} . ", " . $rowTmpDb->{"jcecA3ssPercentage"} . ", " . $rowTmpDb->{"jcecSePercentage"} . ", " . $rowTmpDb->{"jcecRiPercentage"} . ", " . $rowTmpDb->{"jcecMxePercentage"} . ", " . $rowTmpDb->{"jcTotalAsNum"} . ", " . $rowTmpDb->{"jcA5ssPercentage"} . ", " . $rowTmpDb->{"jcA3ssPercentage"} . ", " . $rowTmpDb->{"jcSePercentage"} . ", " . $rowTmpDb->{"jcRiPercentage"} . ", " . $rowTmpDb->{"jcMxePercentage"} . ")";

	$insertAsDb = $dbhAsDb->prepare($sqlAsDb);
	$insertAsDb->execute();

	print $sqlAsDb . "\n";
	<STDIN>;

}
$queryTmpDb->finish();

# input psi into asDB
$sqlTmpDb = "select * from psiTable";
$queryTmpDb = $dbhTmpDb->prepare($sqlTmpDb);
$queryTmpDb->execute();
while($rowTmpDb=$queryTmpDb->fetchrow_hashref()){

	$sqlAsDb = "insert into psiTable (asId, species, tissue, studyId, experiment, JCECI, JCECS, JCECIncFormLen, JCECSkipFormLen, JCECpsi, JCI, JCS, JCIncFormLen, JCSkipFormLen, JCpsi) values(\"" . $rowTmpDb->{"asId"} . "\", \"" . $rowTmpDb->{"species"} . "\", \"" . $rowTmpDb->{"tissue"} . "\", \"" . $rowTmpDb->{"studyId"} . "\", \"" . $rowTmpDb->{"experiment"} . "\", " . $rowTmpDb->{"JCECI"} . ", " . $rowTmpDb->{"JCECS"} . ", " . $rowTmpDb->{"JCECIncFormLen"} . ", " . $rowTmpDb->{"JCECSkipFormLen"} . ", " . $rowTmpDb->{"JCECpsi"} . ", " . $rowTmpDb->{"JCI"} . ", " . $rowTmpDb->{"JCS"} . ", " . $rowTmpDb->{"JCIncFormLen"} . ", " . $rowTmpDb->{"JCSkipFormLen"} . ", " . $rowTmpDb->{"JCpsi"} . ")";
	$insertAsDb = $dbhAsDb->prepare($sqlAsDb);
	$insertAsDb->execute();

	print $sqlAsDb . "\n";
	<STDIN>;

}
$queryTmpDb->finish();

# input variant into asDB
$sqlTmpDb = "select * from variantTable";
$queryTmpDb = $dbhTmpDb->prepare($sqlTmpDb);
$queryTmpDb->execute();
while($rowTmpDb=$queryTmpDb->fetchrow_hashref()){

	$sqlAsDb = "insert into variantTable (species, chr, pos, ref, alt, accId, varType) values(\"" . $rowTmpDb->{"species"} . "\", \"" . $rowTmpDb->{"chr"} . "\", " . $rowTmpDb->{"pos"} . ", \"" . $rowTmpDb->{"ref"} . "\", \"" . $rowTmpDb->{"alt"} . "\", \"" . $rowTmpDb->{"accId"} . "\", \"" . $rowTmpDb->{"varType"} . "\")";
	$insertAsDb = $dbhAsDb->prepare($sqlAsDb);
	$insertAsDb->execute();

	print $sqlAsDb . "\n";
	<STDIN>;

}
$queryTmpDb->finish();

# input isoform into asDB
$sqlTmpDb = "select * from isoformTable";
$queryTmpDb = $dbhTmpDb->prepare($sqlTmpDb);
$queryTmpDb->execute();
while($rowTmpDb=$queryTmpDb->fetchrow_hashref()){

	$sqlAsDb = "insert into isoformTable (isoformId, isoformSymbol, proteinId, geneId, geneSymbol, species, chr, strand, start, end, exonSeries, cDNAseq, pepSeq, DNAseq, annoOrigin, isoformPfam, isoformGo) values(\"" . $rowTmpDb->{"soformId"} . "\", \"" . $rowTmpDb->{"isoformSymbol"} . "\", \"" . $rowTmpDb->{"proteinId"} . "\", \"" . $rowTmpDb->{"geneId"} . "\", \"" . $rowTmpDb->{"geneSymbol"} . "\", \"" . $rowTmpDb->{"species"} . "\", \"" . $rowTmpDb->{"chr"} . "\", \"" . $rowTmpDb->{"strand"} . "\", " . $rowTmpDb->{"start"} . ", " . $rowTmpDb->{"end"} . ", \"" . $rowTmpDb->{"exonSeries"} . "\", \"" . $rowTmpDb->{"cDNAseq"} . "\", \"" . $rowTmpDb->{"pepSeq"} . "\", \"" . $rowTmpDb->{"DNAseq"} . "\", \"" . $rowTmpDb->{"annoOrigin"} . "\", \"" . $rowTmpDb->{"isoformPfam"} . "\", \"" . $rowTmpDb->{"isoformGo"} . "\")";
	$insertAsDb = $dbhAsDb->prepare($sqlAsDb);
	$insertAsDb->execute();

	print $sqlAsDb . "\n";
	<STDIN>;

}
$queryTmpDb->finish();
