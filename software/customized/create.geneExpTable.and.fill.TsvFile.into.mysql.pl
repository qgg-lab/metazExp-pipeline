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
                "--geneExpTsvFile total.GeneExp.tsv \n";
        exit;
}

my ($dbHost, $dbName, $dbUser, $dbPWD, $geneExpTsvFile);
GetOptions(
        'dbHost=s'=>\$dbHost,
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'geneExpTsvFile=s'=>\$geneExpTsvFile,
);

my ($dbh, $sql, $query, $row);
my ($delTable);

$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost", "$dbUser", "$dbPWD");
$dbh->{mysql_auto_reconnect} = 1;

# 删除experimentTable： 首先检查是否存在，如果存在那么删除
$sql = "select * from information_schema.tables where table_name=\"geneExpTable\"";
$query = $dbh->prepare($sql);
$query->execute();
while($row=$query->fetchrow_hashref()){
        $delTable = $dbh->prepare("drop table geneExpTable");
        $delTable->execute();
}

# 创建asTable
$sql = "create table geneExpTable(
        species varchar(50) default \"\",
        taxon varchar(10) default \"\",
        geneId varchar(50) default \"\",
        geneName varchar(50) default \"\",
        chromo char(50) default \"\",
        strand char(1) default \"+\",
        start int default 0,
        end int default 0,
        coverage float default 0,
        fpkm float default 0,
        tpm float default 0,
        experimentId varchar(50) not null,
        primary key(taxon, geneId, experimentId)
)";
 $query = $dbh->prepare($sql);
 $query -> execute();



# 将TSV文件中的数据导入
my ($line);
my ($fieldNameString, $valueString, @fieldName, @value, @field, $field, $value);
open FF, "<$geneExpTsvFile";
# species, taxon, geneId, geneName, chromo, strand, start, end, coverage, fpkm, tpm, experimentId_____"Gallus gallus", "9031", "ENSGALG00000045335", "-", "W", "-", 2988277, 2989653, 1.899027, 0.142906, 0.415449, "SRX1036607"
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);

        $sql = "insert into geneExpTable (" . $fieldNameString . ") values ($valueString)";
        $query = $dbh->prepare($sql);
        $query->execute();
}

$dbh->disconnect();
