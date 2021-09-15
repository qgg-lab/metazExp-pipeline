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
		"--isoformTsvFile total.isoform.tsv \n";
	exit;
}

my ($dbHost, $dbName, $dbUser, $dbPWD, $isoformTsvFile);
GetOptions(
	'dbHost=s'=>\$dbHost,
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'isoformTsvFile=s'=>\$isoformTsvFile,
);

my ($dbh, $sql, $query, $row);
my ($delTable);

$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost", "$dbUser", "$dbPWD");
$dbh->{mysql_auto_reconnect} = 1;

# 删除experimentTable： 首先检查是否存在，如果存在那么删除
$sql = "select count(1) from information_schema.tables where table_name=\"isoformTable\"";
$query = $dbh->prepare($sql);
$query->execute();
while($row=$query->fetchrow_hashref()){
	$delTable = $dbh->prepare("drop table isoformTable");
	$delTable->execute();
}

# 创建isoformTable
$sql = "create table isoformTable(
	isoformId char(50) default \"\",
	species char(50) default \"\",
	isoformSymbol char(50) default \"\",
	geneId char(50) default \"\",
	geneSymbol char(50) default \"\",
	chr varchar(50) default \"\",
	strand char(1) default \"+\",
	start int default 0,
	end int default 0,
	exonSeries varchar(2000) default \"\",
	annoOrigin varchar(50) default \"ensembl\",
	cDNAseq text,
	DNAseq text,
	pepSeq text,
	isoformPfam varchar(500) default \"\",
	isoformGo varchar(500) default \"\",
	proteinId char(50) default \"\",
	primary key(isoformId, species)
)";

 $query = $dbh->prepare($sql);
 $query -> execute();



# 将TSV文件中的数据导入
my ($line);
my ($fieldNameString, $valueString, @fieldName, @value, @field, $field, $value);
open FF, "<$isoformTsvFile";
# isoformId, species, isoformSymbol, geneId, geneSymbol, chr, strand, start, end, exonSeries, annoOrigin, cDNAseq, DNAseq, pepSeq, isoformPfam, isoformGo, proteinId_____"SRX196373.29379.1", "Gallus gallus", "NA", "ENSGALG00000002145", "NA", "8", "-", 1814854, 1825930, "1825249..1825930,1814854..1815161", "RNAseq", ...  
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);
	$sql = "insert into isoformTable (" . $fieldNameString . ") values ($valueString)";
	$query = $dbh->prepare($sql);
	$query->execute();
}

$dbh->disconnect();
