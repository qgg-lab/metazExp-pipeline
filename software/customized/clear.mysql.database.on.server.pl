#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD njaulsas2019 \\\n" .
		"--species \"all\" \n\n\n";
	exit;
}

my ($dbName, $dbUser, $dbPWD, $species);
GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'species=s'=>\$species,
);

my ($dbh, $query, $sql);
$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;

if(uc($species) eq "ALL"){
	$sql = "delete from asTable";
}else{
	$sql = "delete from asTable where species=\"$species\"";
}
$query = $dbh->prepare($sql);
$query->execute();

if(uc($species) eq "ALL"){
	$sql = "delete from experimentTable";
}else{
	$sql = "delete from experimentTable where species=\"$species\"";
}
$query = $dbh->prepare($sql);
$query->execute();

if(uc($species) eq "ALL"){
	$sql = "delete from isoformTable";
}else{
	$sql = "delete from isoformTable where species=\"$species\"";
}
$query = $dbh->prepare($sql);
$query->execute();

if(uc($species) eq "ALL"){
	$sql = "delete from psiTable";
}else{
	$sql = "delete from psiTable where species=\"$species\"";
}
$query = $dbh->prepare($sql);
$query->execute();

if(uc($species) eq "ALL"){
	$sql = "delete from variantTable";
}else{
	$sql = "delete from variantTable where species=\"$species\"";
}
$query = $dbh->prepare($sql);
$query->execute();

$dbh->disconnect();
