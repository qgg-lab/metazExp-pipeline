#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--dbHost localhost \\\n" .
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD njaulsas2019 \\\n" .
		"--asTsvFile total.alternative.splicing.annotation.tsv \n";
	exit;
}

my ($dbHost, $dbName, $dbUser, $dbPWD, $asTsvFile);
GetOptions(
	'dbHost=s'=>\$dbHost,
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'asTsvFile=s'=>\$asTsvFile,
);

my ($dbh, $sql, $query, $row);
my ($delTable);

$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost", "$dbUser", "$dbPWD");
$dbh->{mysql_auto_reconnect} = 1;

# 删除asTable： 首先检查是否存在，如果存在那么删除
$sql = "select count(1) from information_schema.tables where table_name=\"asTable\"";
$query = $dbh->prepare($sql);
$query->execute();
while($row=$query->fetchrow_hashref()){
	$delTable = $dbh->prepare("drop table asTable");
	$delTable->execute();
}

# 创建asTable
$sql = "create table asTable(
	asId char(50) primary key,
	asType char(5) not null default \"\",
	species varchar(50) not null default \"\",
	chr varchar(50) not null default \"\",
	strand char(1) default \"+\",
	start int default 0,
	end int default 0,
	1stExonEnd int default 0,
	1stExonStart_0base int default 0,
	2ndExonEnd int default 0,
	2ndExonStart_0base int default 0,
	downstreamEE int default 0,
	downstreamES int default 0,
	exonEnd int default 0,
	exonStart_0base int default 0,
	flankingEE int default 0,
	flankingES int default 0,
	longExonEnd int default 0,
	longExonStart_0base int default 0,
	riExonEnd int default 0,
	riExonStart_0base int default 0,
	shortEE int default 0,
	shortES int default 0,
	upstreamEE int default 0,
	upstreamES int default 0,
	geneID varchar(50) default \"\",
	geneSymbol varchar(50) default \"\",
	firstAltExonSeries varchar(200) default \"\",
	firstIsoformIdList varchar(1000) default \"\",
	firstIsoformNum int default 0,
	secondAltExonSeries varchar(200) default \"\",
	secondIsoformIdList varchar(1000) default \"\",
	secondIsoformNum int default 0,
	pfam varchar(200) default \"\",
	go varchar(200) default \"\",
	variantPointNum int default 0,
	variantPointTypeCmb varchar(5) default \"\",
	asOrthId varchar(50) default \"\",
	conservedSpeciesNum int default 0,
	jcecExperimentNum int default 0,
	jcExperimentNum int default 0,
	discoveryApproach varchar(25) default \"ensembl\",
	orthAsId varchar(20) default \"\",
	orthAsTaxon varchar(100) default \"\",
	TissueExpNumList int default 0,
	TissueNum int default 0,
	splicingExpEVID char(1) default \"Y\",
)";
# 注意：虽然asOrthId和conservedSpeciesNum失效没有具体值，但是保留到asTable表格中
 $query = $dbh->prepare($sql);
 $query -> execute();



# 将TSV文件中的数据导入
my ($line);
my ($fieldNameString, $valueString, @fieldName, @value, $field, $value);
open FF, "<$asTsvFile";
# asId, asType, species, chr, strand, start, end, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base, downstreamEE, downstreamES, exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE, shortES, upstreamEE, upstreamES, geneID, geneSymbol, firstAltExonSeries, firstIsoformIdList, firstIsoformNum, secondAltExonSeries, secondIsoformIdList, secondIsoformNum, pfam, go, variantPointNum, variantPointTypeCmb, asOrthId, conservedSpeciesNum, jcecExperimentNum, jcExperimentNum, discoveryApproach, orthAsId, orthAsTaxon, TissueExpNumList, TissueNum, splicingExpEVID_____"GGALMXE0000034850", "MXE", "Gallus gallus", "20", "-", 719031, 753202, 722489, 721992, 730499, 730443, 753202, 753039, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 719178, 719030, "ENSGALG00000032417", "CHD6", "753040..753202,721993..722489,719031..719178", "ERX1751702.14362.20,SRX893872.18044.18,SRX529337.12111.16", 3, "753040..753202,730444..730499,719031..719178", "", 0, "PF00176,PF00385,PF07533,PF00271", "GO:0005524,GO:0005515,GO:0016817", 1115, "SDI", "", 0, 7, 7, "Novel", "", "", "parietal nerve(1),embryos(1),muscle(2),embryo wing bud(1),gonads(2)", 5, "Y"
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);
	$sql = "insert into asTable (" . $fieldNameString . ") values ($valueString)";
	$query = $dbh->prepare($sql);
	$query->execute();
}

$dbh->disconnect();
