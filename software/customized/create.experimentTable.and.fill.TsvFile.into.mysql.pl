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
		"--experimentTsvFile experiment.PSI.tsv \n";
	exit;
}

my ($dbHost, $dbName, $dbUser, $dbPWD, $experimentTsvFile);
GetOptions(
	'dbHost=s'=>\$dbHost,
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'experimentTsvFile=s'=>\$experimentTsvFile,
);

my ($dbh, $sql, $query, $row);
my ($delTable);

$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost", "$dbUser", "$dbPWD");
$dbh->{mysql_auto_reconnect} = 1;

# 删除experimentTable： 首先检查是否存在，如果存在那么删除
$sql = "select count(1) from information_schema.tables where table_name=\"experimentTable\"";
$query = $dbh->prepare($sql);
$query->execute();
while($row=$query->fetchrow_hashref()){
	$delTable = $dbh->prepare("drop table experimentTable");
	$delTable->execute();
}

# 创建asTable
$sql = "create table experimentTable(
	alignPercent float default 0,
	mappedSpots float default 0,
	experimentId char(10) primary key,
	species char(50) default \"\",
	libraryType varchar(50) default \"\",
	libraryLayout varchar(50) default \"\",
	readLen varchar(50) default \"\",
	phredScore varchar(50) default \"\",
	totalSpots int default 0,
	runIdList varchar(100) default \"\",
	runNum int default 0,
	studyId varchar(10) default \"\",
	jcecTotalAsNum int default 0,
	jcecA5ssPercentage float default 0,
	jcecA3ssPercentage float default 0,
	jcecSePercentage float default 0,
	jcecRiPercentage float default 0,
	jcecMxePercentage float default 0,
	jcTotalAsNum int default 0,
	jcA5ssPercentage float default 0,
	jcA3ssPercentage float default 0,
	jcSePercentage float default 0,
	jcRiPercentage float default 0,
	jcMxePercentage float default 0,
	tissue varchar(50) default \"\"
)";

 $query = $dbh->prepare($sql);
 $query -> execute();



# 将TSV文件中的数据导入
my ($line);
my ($fieldNameString, $valueString, @fieldName, @value, @field, $field, $value);
open FF, "<$experimentTsvFile";
# alignPercent, mappedSpots, experimentId, species, libraryType, libraryLayout, readLen, phredScore, totalSpots, runIdList, runNum, studyId, jcecTotalAsNum, jcecA5ssPercentage, jcecA3ssPercentage, jcecSePercentage, jcecRiPercentage, jcecMxePercentage, jcTotalAsNum, jcA5ssPercentage, jcA3ssPercentage, jcSePercentage, jcRiPercentage, jcMxePercentage, tissue_____94.41, 66.05, "SRX1036607", "Gallus gallus", "UN", "PAIRED", "101", "33", 69.97, "SRR2037196", 1, "SRP058621", 68536, 11.27, 13.67, 42.43, 28.47, 4.14, 68007, 11.25, 13.71, 42.47, 28.41, 4.13, "embryos"
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);

	@value = ();
	@value = split(/, /, $valueString);
	# 处理library, layout, readLen
	$value[4] = substr($value[4], 1, length($value[4])-2);
	@field = ();
	@field = split(/,/, $value[4]);
	$value[4] = "\"" . $field[0] . "\"";

	$value[5] = substr($value[5], 1, length($value[5])-2);
	@field = ();
	@field = split(/,/, $value[5]);
	$value[5] = "\"" . $field[0] . "\"";

	$value[6] = substr($value[6], 1, length($value[6])-2);
	@field = ();
	@field = split(/,/, $value[6]);
	$value[6] = "\"" . $field[0] . "\"";

        $value[7] = substr($value[7], 1, length($value[7])-2);
        @field = ();
        @field = split(/,/, $value[7]);
        $value[7] = "\"" . $field[0] . "\"";

	$valueString = join(", ", @value);
	$sql = "insert into experimentTable (" . $fieldNameString . ") values ($valueString)";
	$query = $dbh->prepare($sql);
	$query->execute();
}

$dbh->disconnect();
