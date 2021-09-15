#!/usr/bin/perl
use strict;
use Getopt::Long;
 
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--host   localhost \\\n" .
		"--dbName asdb  \\\n" .
		"--dbUser lsas	\\\n" .
		"--dbPwd  njaulsas2019 \\\n\n";
	exit;
}
my ($host, $dbName, $dbUser, $dbPwd);

GetOptions(
	'host=s'=>\$host,
	'dbName=s'=>\$dbName,
	'dbUser=s'=>\$dbUser,
	'dbPwd=s'=>\$dbPwd,
);

my ($dbh, $sql, $query);
$dbh = DBI->connect("DBI:mysql:database=$dbName;host=$host", $dbUser, $dbPwd);
=cut;
############ -- 1: asTable #############
#
# 	create asTable
# 					
########################################
$sql = "create table asTable(
	asId char(15) primary key,
	1stExonEnd int default 0,
	1stExonStart_0base int default 0,
	2ndExonEnd int default 0,
	2ndExonStart_0base int default 0,
	chr int default 0,
	downstreamEE int default 0,
	downstreamES int default 0,
	exonEnd int default 0,
	exonStart_0base int default 0,
	flankingEE int default 0,
	flankingES int default 0,
	GeneID int default 0,
	geneSymbol int default 0,
	longExonEnd int default 0,
	longExonStart_0base int default 0,
	riExonEnd int default 0,
	riExonStart_0base int default 0,
	shortEE int default 0,
	shortES int default 0,
	strand char(1) default \"+\",
	upstreamEE int default 0,
	upstreamES int default 0
)";

$query = $dbh->prepare($sql);
$query->execute();

