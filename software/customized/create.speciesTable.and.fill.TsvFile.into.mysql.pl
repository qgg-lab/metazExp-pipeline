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
                "--speciesTsvFile total.species.info.tsv \n";
        exit;
}

my ($dbHost, $dbName, $dbUser, $dbPWD, $speciesTsvFile);
GetOptions(
        'dbHost=s'=>\$dbHost,
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'speciesTsvFile=s'=>\$speciesTsvFile,
);

my ($dbh, $sql, $query, $row);
my ($delTable);

$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost", "$dbUser", "$dbPWD");
$dbh->{mysql_auto_reconnect} = 1;

# 删除experimentTable： 首先检查是否存在，如果存在那么删除
$sql = "select count(1) from information_schema.tables where table_name=\"speciesTable\"";
$query = $dbh->prepare($sql);
$query->execute();
while($row=$query->fetchrow_hashref()){
        $delTable = $dbh->prepare("drop table speciesTable");
        $delTable->execute();
}

$sql = "create table speciesTable(
        species varchar(50),
        taxon char(50) primary key,
        genomeAssemblyName varchar(100) default \"\",
        genomeAssemblyAddress varchar(500) default \"\",
        geneNum int default 0,
        isoformNum int default 0,
        isoformAnnoFileName varchar(100) default \"\",
        ensemblAnnotationName varchar(100) default \"\",
        ensemblAnnotationAddress varchar(500) default \"\",
        refseqAnnotationName varchar(100) default \"\",
        refseqAnnotationAddress varchar(500) default \"\",
        jgiAnnotationName varchar(100) default \"\",
        jgiAnnotationAddress varchar(500) default \"\",
        annoBasedRNAseqExpNum int default 0,
        annoBasedRNAseqExpList text,
        ensemblIsoformNum int default 0,
        refseqIsoformNum int default 0,
        jgiIsoformNum int default 0,
        RNAseqIsoformNum int default 0,
        totalRNAseqExpNum int  default 0,
        totalRNAseqExpList text,
        totalStudyNum int default 0,
        totalStudyList text,
        asNum int default 0,
        A5ssPercentage float default 0,
        A3ssPercentage float default 0,
        SePercentage float default 0,
        RiPercentage float default 0,
        MxePercentage float default 0,
        asNumByEnsembl int default 0,
        asA5ssNumByEnsembl int default 0,
        asA3ssNumByEnsembl int default 0,
        asSeNumByEnsembl int default 0,
        asRiNumByEnsembl int default 0,
        asMxeNumByEnsembl int default 0,
        asNumByEnsembl_Refseq int default 0,
        asA5ssNumByEnsembl_Refseq int default 0,
        asA3ssNumByEnsembl_Refseq int default 0,
        asSeNumByEnsembl_Refseq int default 0,
        asRiNumByEnsembl_Refseq int default 0,
        asMxeNumByEnsembl_Refseq int default 0,
        asNumByEnsembl_Refseq_RNAseq int default 0,
        asA5ssNumByEnsembl_Refseq_RNAseq int default 0,
        asA3ssNumByEnsembl_Refseq_RNAseq int default 0,
        asSeNumByEnsembl_Refseq_RNAseq int default 0,
        asRiNumByEnsembl_Refseq_RNAseq int default 0,
        asMxeNumByEnsembl_Refseq_RNAseq int default 0,
        asNumByNovel int default 0,
        asA5ssNumByNovel int default 0,
        asA3ssNumByNovel int default 0,
        asSeNumByNovel int default 0,
        asRiNumByNovel int default 0,
        asMxeNumByNovel int default 0
)";
 $query = $dbh->prepare($sql);
 $query -> execute();



# 将TSV文件中的数据导入
my ($line);
my ($fieldNameString, $valueString, @fieldName, @value, @field, $field, $value);
open FF, "<$speciesTsvFile";

while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);
        $sql = "insert into speciesTable (" . $fieldNameString . ") values ($valueString)";
#       print $sql;
#       <STDIN>;
        $query = $dbh->prepare($sql);
        $query->execute();
}

$dbh->disconnect();
