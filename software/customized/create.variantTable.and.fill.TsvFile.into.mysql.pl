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
                "--variantTsvFile total.variant.tsv \n";
        exit;
}

my ($dbHost, $dbName, $dbUser, $dbPWD, $variantTsvFile);
GetOptions(
        'dbHost=s'=>\$dbHost,
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'variantTsvFile=s'=>\$variantTsvFile,
);

my ($dbh, $sql, $query, $row);
my ($delTable);

$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost", "$dbUser", "$dbPWD");
$dbh->{mysql_auto_reconnect} = 1;

# 删除experimentTable： 首先检查是否存在，如果存在那么删除
$sql = "select count(1) from information_schema.tables where table_name=\"variantTable\"";
$query = $dbh->prepare($sql);
$query->execute();
while($row=$query->fetchrow_hashref()){
        $delTable = $dbh->prepare("drop table variantTable");
        $delTable->execute();
}

# 创建asTable
$sql = "create table variantTable(
        species char(50) default \"\",
        chr char(50) default \"\",
        pos int default 0,
        ref varchar(100),
        alt varchar(100),
        accId varchar(50),
        varType varchar(50),
        hitElement int default 0,
        primary key(species, chr, pos, ref, alt, accId)
)";
 $query = $dbh->prepare($sql);
 $query -> execute();



# 将TSV文件中的数据导入
my ($line);
my ($fieldNameString, $valueString, @fieldName, @value, @field, $field, $value);
open FF, "<$variantTsvFile";
# species, chr, pos, ref, alt, accId, varType, hitElement_____"Gallus gallus", "1", 194431, "T", "C", "rs1057644598", "SNV", 0
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);
        $sql = "insert into variantTable (" . $fieldNameString . ") values ($valueString)";
        $query = $dbh->prepare($sql);
        $query->execute();
}

$dbh->disconnect();

