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
                "--orthExonTsvFile orthExon.tsv \n";
        exit;
}

my ($dbHost, $dbName, $dbUser, $dbPWD, $orthExonTsvFile);
GetOptions(
        'dbHost=s'=>\$dbHost,
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'orthExonTsvFile=s'=>\$orthExonTsvFile,
);

my ($dbh, $sql, $query, $row);
my ($delTable);

$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost", "$dbUser", "$dbPWD");
$dbh->{mysql_auto_reconnect} = 1;

# 删除experimentTable： 首先检查是否存在，如果存在那么删除
$sql = "select count(1) from information_schema.tables where table_name=\"orthExonTable\"";
$query = $dbh->prepare($sql);
$query->execute();
while($row=$query->fetchrow_hashref()){
        $delTable = $dbh->prepare("drop table orthExonTable");
        $delTable->execute();
}

# 创建asTable
$sql = "create table orthExonTable(
        taxonId char(10),
        orthExonId char(50),
        chr char(50),
        strand char(1) default \"+\",
        start int default 0,
        stop int default 0,
        primary key(taxonId, orthExonId)
)";

 $query = $dbh->prepare($sql);
 $query -> execute();

# 将TSV文件中的数据导入
my ($line);
my ($fieldNameString, $valueString, @fieldName, @value, @field, $field, $value);
open FF, "<$orthExonTsvFile";
# taxonId, orthExonId, chr, strand, start, stop_____"9031", "ORTH0000154419", "1", "-", 6754, 6829
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);

        $sql = "insert into orthExonTable (" . $fieldNameString . ") values ($valueString)";
        $query = $dbh->prepare($sql);
        $query->execute();
}

$dbh->disconnect();
