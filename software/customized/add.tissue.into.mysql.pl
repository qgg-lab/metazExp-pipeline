#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD njaulsas2019 \\\n" .
		"--inputSampleFile sample.tsv \n";
	exit;
}

my ($dbName, $dbUser, $dbPWD, $inputSampleFile);
GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'inputSampleFile=s'=>\$inputSampleFile,
);

my $dbh = DBI->connect("DBI:mysql:database=$dbName", $dbUser, $dbPWD);
$dbh->{mysql_auto_reconnect} = 1;

my ($line, $experimentId, $tissue, $update);
open FF, "<$inputSampleFile";
while($line=<FF>){
	chomp($line);
	($experimentId,  $tissue) = split(/\t/, $line);
#	$update = $dbh->prepare("update psiTable set tissue=\"" . $tissue . "\" where experiment=\"" . $experimentId . "\"");
#	$update->execute();
	$update = $dbh->prepare("update experimentTable set tissue=\"" . $tissue . "\" where experimentId=\"" . $experimentId . "\"");
	$update->execute();

}
$dbh->disconnect();
