#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD njaulsas2019 \\\n" .
		"--sampleListForIsoformAnno ../004-combine-assemblies-and-annos/combined.transcriptome.list \\\n". 
		"--sampleLisftForPickPsi ./samples.filtered.by.sequencing.and.mapping.tsv \\\n" .
		"--species \"Bubalus bubalis\" \n"; 
	exit;
}

my ($species, $sampleListForIsoformAnno, $sampleLisftForPickPsi, $sql);
my ($dbName, $dbUser, $dbPWD);
GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'species=s'=>\$species,
	'sampleListForIsoformAnno=s'=>\$sampleListForIsoformAnno,
	'sampleLisftForPickPsi=s'=>\$sampleLisftForPickPsi,
);

my $dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
my ($query, $row);

# geneNum, isoformNum, isoformAnnoFileName
# ensemblIsoformNum, refseqIsoformNum, jgiIsoformNum, RNAseqIsoformNum
my (%geneId, @geneId, $geneNum, $isoformNum);
my ($ensemblIsoformNum, $refseqIsoformNum, $jgiIsoformNum, $RNAseqIsoformNum);
$sql = "select geneId, isoformId, annoOrigin from isoformTable where species=\"" . $species . "\"";
$query = $dbh->prepare($sql);
$query->execute();
$isoformNum=0;
$ensemblIsoformNum=0;
$refseqIsoformNum=0;
$jgiIsoformNum=0;
$RNAseqIsoformNum=0;

while($row=$query->fetchrow_hashref()){
	$geneId{$row->{"geneId"}}=1;
	$isoformNum++;
	if(uc($row->{"annoOrigin"}) eq "ENSEMBL"){
		$ensemblIsoformNum++
	}elsif(uc($row->{"annoOrigin"}) eq "REFSEQ"){
		$refseqIsoformNum++;
	}elsif(uc($row->{"annoOrigin"}) eq "JGI"){
		$jgiIsoformNum++;
	}elsif(uc($row->{"annoOrigin"}) eq "RNASEQ"){
		$RNAseqIsoformNum++;
	}
}
@geneId = keys(%geneId);
$geneNum = $#geneId + 1;

# annoBasedRNAseqExpNum, annoBasedRNAseqExpList
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

# totalRNAseqExpNum, totalRNAseqExpList
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

# asNum, A5ssPercentage, A3ssPercentage, SePercentage, RiPercentage, MxePercentage
my ($asNum, $A5ssPercentage, $A3ssPercentage, $SePercentage, $RiPercentage, $MxePercentage);
my ($A5ssNum, $A3ssNum, $SeNum, $RiNum, $MxeNum);
$asNum=0;
$A5ssPercentage=0;
$A3ssPercentage=0;
$SePercentage=0;
$RiPercentage=0;
$MxePercentage=0;

$query = $dbh->prepare("select count(*) as asNum, asType from asTable where species=\"" . $species . "\" group by asType");
$query->execute();
while($row=$query->fetchrow_hashref()){
	if(uc($row->{"asType"}) eq "A5SS"){
		$A5ssNum = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "A3SS"){
		$A3ssNum = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "SE"){
		$SeNum = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "RI"){
		$RiNum = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "MXE"){
		$MxeNum = $row->{"asNum"};
	}
}
$A5ssPercentage = sprintf("%.2f",$A5ssNum/($A5ssNum+$A3ssNum+$SeNum+$RiNum+$MxeNum));
$A3ssPercentage = sprintf("%.2f", $A3ssNum/($A5ssNum+$A3ssNum+$SeNum+$RiNum+$MxeNum));
$SePercentage = sprintf("%.2f", $SeNum/($A5ssNum+$A3ssNum+$SeNum+$RiNum+$MxeNum));
$RiPercentage = sprintf("%.2f", $RiNum/($A5ssNum+$A3ssNum+$SeNum+$RiNum+$MxeNum));
$MxePercentage = sprintf("%.2f", $MxeNum/($A5ssNum+$A3ssNum+$SeNum+$RiNum+$MxeNum));
$asNum = $A5ssNum+$A3ssNum+$SeNum+$RiNum+$MxeNum;

# asNumByEnsembl, asA5ssNumByEnsembl, asA3ssNumByEnsembl, asSeNumByEnsembl, asRiNumByEnsembl, asMxeNumByEnsembl
my ($asNumByEnsembl, $asA5ssNumByEnsembl, $asA3ssNumByEnsembl, $asSeNumByEnsembl, $asRiNumByEnsembl, $asMxeNumByEnsembl);
$asNumByEnsembl = 0;
$asA5ssNumByEnsembl = 0;
$asA3ssNumByEnsembl = 0;
$asSeNumByEnsembl = 0;
$asRiNumByEnsembl = 0;
$asMxeNumByEnsembl = 0;
$query = $dbh->prepare("select count(*) as asNumByEnsembl from asTable where species=\"" . $species . "\" and discoveryApproach=\"Ensembl\"");
$query->execute();
$row=$query->fetchrow_hashref();
$asNumByEnsembl = $row->{"asNumByEnsembl"};

$query = $dbh->prepare("select count(*) as asNum, asType from asTable where species=\"" . $species . "\" and discoveryApproach=\"Ensembl\" group by asType");
$query->execute();
while($row=$query->fetchrow_hashref()){
	if(uc($row->{"asType"}) eq "A5SS"){
		$asA5ssNumByEnsembl = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "A3SS"){
		$asA3ssNumByEnsembl = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "SE"){
		$asSeNumByEnsembl = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "RI"){
		$asRiNumByEnsembl = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "MXE"){
		$asMxeNumByEnsembl = $row->{"asNum"};
	}
}

# asNumByEnsembl_Refseq, asA5ssNumByEnsembl_Refseq, asA3ssNumByEnsembl_Refseq, asSeNumByEnsembl_Refseq, asRiNumByEnsembl_Refseq, asMxeNumByEnsembl_Refseq
my ($asNumByEnsembl_Refseq, $asA5ssNumByEnsembl_Refseq, $asA3ssNumByEnsembl_Refseq, $asSeNumByEnsembl_Refseq, $asRiNumByEnsembl_Refseq, $asMxeNumByEnsembl_Refseq);
$asNumByEnsembl_Refseq = 0;
$asA5ssNumByEnsembl_Refseq = 0;
$asA3ssNumByEnsembl_Refseq = 0;
$asSeNumByEnsembl_Refseq = 0;
$asRiNumByEnsembl_Refseq = 0;
$asMxeNumByEnsembl_Refseq = 0;
$query = $dbh->prepare("select count(*) as asNumByEnsembl from asTable where species=\"" . $species . "\" and discoveryApproach=\"Ensembl_Refseq\"");
$query->execute();
$row=$query->fetchrow_hashref();
$asNumByEnsembl_Refseq = $row->{"asNumByEnsembl"};

$query = $dbh->prepare("select count(*) as asNum, asType from asTable where species=\"" . $species . "\" and discoveryApproach=\"Ensembl_Refseq\" group by asType");
$query->execute();
while($row=$query->fetchrow_hashref()){
	if(uc($row->{"asType"}) eq "A5SS"){
		$asA5ssNumByEnsembl_Refseq = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "A3SS"){
		$asA3ssNumByEnsembl_Refseq = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "SE"){
		$asSeNumByEnsembl_Refseq = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "RI"){
		$asRiNumByEnsembl_Refseq = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "MXE"){
		$asMxeNumByEnsembl_Refseq = $row->{"asNum"};
	}
}

# asNumByEnsembl_Refseq_RNAseq,asA5ssNumByEnsembl_Refseq_RNAseq, asA3ssNumByEnsembl_Refseq_RNAseq, asSeNumByEnsembl_Refseq_RNAseq, asRiNumByEnsembl_Refseq_RNAseq, asMxeNumByEnsembl_Refseq_RNAseq
my ($asNumByEnsembl_Refseq_RNAseq, $asA5ssNumByEnsembl_Refseq_RNAseq, $asA3ssNumByEnsembl_Refseq_RNAseq, $asSeNumByEnsembl_Refseq_RNAseq, $asRiNumByEnsembl_Refseq_RNAseq, $asMxeNumByEnsembl_Refseq_RNAseq);

$asNumByEnsembl_Refseq_RNAseq=0;
$asA5ssNumByEnsembl_Refseq_RNAseq=0;
$asA3ssNumByEnsembl_Refseq_RNAseq=0;
$asSeNumByEnsembl_Refseq_RNAseq=0;
$asRiNumByEnsembl_Refseq_RNAseq=0;
$asMxeNumByEnsembl_Refseq_RNAseq=0;
$query = $dbh->prepare("select count(*) as asNumByEnsembl from asTable where species=\"" . $species . "\" and discoveryApproach=\"Ensembl_Refseq_RNAseq\"");
$query->execute();
$row=$query->fetchrow_hashref();
$asNumByEnsembl_Refseq_RNAseq = $row->{"asNumByEnsembl"};

$query = $dbh->prepare("select count(*) as asNum, asType from asTable where species=\"" . $species . "\" and discoveryApproach=\"Ensembl_Refseq_RNAseq\" group by asType");
$query->execute();
while($row=$query->fetchrow_hashref()){
	if(uc($row->{"asType"}) eq "A5SS"){
		$asA5ssNumByEnsembl_Refseq_RNAseq = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "A3SS"){
		$asA3ssNumByEnsembl_Refseq_RNAseq = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "SE"){
		$asSeNumByEnsembl_Refseq_RNAseq = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "RI"){
		$asRiNumByEnsembl_Refseq_RNAseq = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "MXE"){
		$asMxeNumByEnsembl_Refseq_RNAseq = $row->{"asNum"};
	}
}



# asNumByNovel,asA5ssNumByNovel, asA3ssNumByNovel, asSeNumByNovel, asRiNumByNovel, asMxeNumByNovel
my ($asNumByNovel, $asA5ssNumByNovel, $asA3ssNumByNovel, $asSeNumByNovel, $asRiNumByNovel, $asMxeNumByNovel);
$asNumByNovel=0;
$asA5ssNumByNovel=0;
$asA3ssNumByNovel=0;
$asSeNumByNovel=0;
$asRiNumByNovel=0;
$asMxeNumByNovel=0;
$query = $dbh->prepare("select count(*) as asNumByEnsembl from asTable where species=\"" . $species . "\" and discoveryApproach=\"Novel\"");
$query->execute();
$row=$query->fetchrow_hashref();
$asNumByNovel = $row->{"asNumByEnsembl"};

$query = $dbh->prepare("select count(*) as asNum, asType from asTable where species=\"" . $species . "\" and discoveryApproach=\"Novel\" group by asType");
$query->execute();
while($row=$query->fetchrow_hashref()){
	if(uc($row->{"asType"}) eq "A5SS"){
		$asA5ssNumByNovel = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "A3SS"){
		$asA3ssNumByNovel = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "SE"){
		$asSeNumByNovel = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "RI"){
		$asRiNumByNovel = $row->{"asNum"};
	}elsif(uc($row->{"asType"}) eq "MXE"){
		$asMxeNumByNovel = $row->{"asNum"};
	}
}

my $sql = "";
$sql = "update speciesTable set geneNum=$geneNum, isoformNum=$isoformNum, ensemblIsoformNum=$ensemblIsoformNum, refseqIsoformNum=$refseqIsoformNum, jgiIsoformNum=$jgiIsoformNum, RNAseqIsoformNum=$RNAseqIsoformNum, annoBasedRNAseqExpNum=$annoBasedRNAseqExpNum, annoBasedRNAseqExpList=\"$annoBasedRNAseqExpList\", totalRNAseqExpNum=$totalRNAseqExpNum, totalRNAseqExpList=\"$totalRNAseqExpList\", totalStudyNum=$totalStudyNum, totalStudyList=\"$totalStudyList\", asNum=$asNum, A5ssPercentage=$A5ssPercentage, A3ssPercentage=$A3ssPercentage, SePercentage=$SePercentage, RiPercentage=$RiPercentage, MxePercentage=$MxePercentage, asNumByEnsembl=$asNumByEnsembl, asA5ssNumByEnsembl=$asA5ssNumByEnsembl, asA3ssNumByEnsembl=$asA3ssNumByEnsembl, asSeNumByEnsembl=$asSeNumByEnsembl, asRiNumByEnsembl=$asRiNumByEnsembl, asMxeNumByEnsembl=$asMxeNumByEnsembl, asNumByEnsembl_Refseq=$asNumByEnsembl_Refseq, asA5ssNumByEnsembl_Refseq=$asA5ssNumByEnsembl_Refseq, asA3ssNumByEnsembl_Refseq=$asA3ssNumByEnsembl_Refseq, asSeNumByEnsembl_Refseq=$asSeNumByEnsembl_Refseq, asRiNumByEnsembl_Refseq=$asRiNumByEnsembl_Refseq, asMxeNumByEnsembl_Refseq=$asMxeNumByEnsembl_Refseq, asNumByEnsembl_Refseq_RNAseq=$asNumByEnsembl_Refseq_RNAseq, asA5ssNumByEnsembl_Refseq_RNAseq=$asA5ssNumByEnsembl_Refseq_RNAseq, asA3ssNumByEnsembl_Refseq_RNAseq=$asA3ssNumByEnsembl_Refseq_RNAseq, asSeNumByEnsembl_Refseq_RNAseq=$asSeNumByEnsembl_Refseq_RNAseq, asRiNumByEnsembl_Refseq_RNAseq=$asRiNumByEnsembl_Refseq_RNAseq, asMxeNumByEnsembl_Refseq_RNAseq=$asMxeNumByEnsembl_Refseq_RNAseq, asNumByNovel=$asNumByNovel, asA5ssNumByNovel=$asA5ssNumByNovel, asA3ssNumByNovel=$asA3ssNumByNovel, asSeNumByNovel=$asSeNumByNovel, asRiNumByNovel=$asRiNumByNovel, asMxeNumByNovel=$asMxeNumByNovel where species=\"" . $species . "\"";
print $sql;
<STDIN>;
$query = $dbh->prepare($sql);
$query->execute();
