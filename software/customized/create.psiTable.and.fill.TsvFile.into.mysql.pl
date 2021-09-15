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
		"--psiTsvFile total.PSI.tsv \n";
	exit;
}

my ($dbHost, $dbName, $dbUser, $dbPWD, $psiTsvFile);
GetOptions(
	'dbHost=s'=>\$dbHost,
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'psiTsvFile=s'=>\$psiTsvFile,
);

my ($dbh, $sql, $query, $row);
my ($delTable);

$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost", "$dbUser", "$dbPWD");
$dbh->{mysql_auto_reconnect} = 1;

# 删除psiTable： 首先检查是否存在，如果存在那么删除
$sql = "select count(1) from information_schema.tables where table_name=\"psiTable\"";
$query = $dbh->prepare($sql);
$query->execute();
while($row=$query->fetchrow_hashref()){
	$delTable = $dbh->prepare("drop table psiTable");
	$delTable->execute();
}

# 创建asTable
$sql = "create table psiTable(
	asId char(50) not null,
	species char(50) default \"\",
	experiment char(10) default \"\",
	studyId char(10) default \"\",
	JCECI int default 0,
	JCECS int default 0,
	JCECIncFormLen int default 0,
	JCECSkipFormLen int default 0,
	JCECpsi float default 0,
	JCI int default 0,
	JCS int default 0,
	JCIncFormLen int default 0,
	JCSkipFormLen int default 0,
	JCpsi float default 0,
	primary key(asId, experiment)
)";
# 注意：虽然asOrthId和conservedSpeciesNum失效没有具体值，但是保留到asTable表格中
 $query = $dbh->prepare($sql);
 $query -> execute();



# 将TSV文件中的数据导入
my ($line);
my ($fieldNameString, $valueString, @fieldName, @value, $field, $value);
open FF, "<$psiTsvFile";
# asId, species, experiment, studyId, JCECI, JCECS, JCECIncFormLen, JCECSkipFormLen, JCECpsi, JCI, JCS, JCIncFormLen, JCSkipFormLen, JCpsi_____"GGALMXE0000034850", "Gallus gallus", "SRX2406013", "SRP094786", 9, 627, 206, 647, 0.0431381264584954, 9, 135, 205, 300, 0.0888888888888889 
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);
        @value = ();
        @value = split(/, /, $valueString);
        for(my $i=4; $i<=$#value; $i++){
                if($value[$i] eq ""){
                        $value[$i] = 0;
                }
        }
        $valueString = join(", ", @value);

	$sql = "insert into psiTable (" . $fieldNameString . ") values ($valueString)";
	$query = $dbh->prepare($sql);
	$query->execute();
}

$dbh->disconnect();
