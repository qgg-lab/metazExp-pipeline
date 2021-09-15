#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD njaulsas2019 \\\n" .
		"--asFileList  new.A3SS.catalog,new.A5SS.catalog,new.MXE.catalog,new.RI.catalog,new.SE.catalog\n";
	exit;
}

my ($dbName, $dbUser, $dbPWD);
my ($asFileList, @asFile, $asFile, $line, @fields, $asId, $symbol, $sql);
GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
	'asFileList=s'=>\$asFileList,
);

my ($dbh, $update);
$dbh = DBI->connect("DBI:mysql:database=$dbName", $dbUser, $dbPWD);
$dbh->{mysql_auto_reconnect} = 1;
@asFile = split(/\t/, $asFileList);
foreach $asFile(@asFile){
	open FF, "<$asFile";
	while($line=<FF>){
		#ASID                    GeneID                  geneSymbol      chr
		#GGALA3SS0000004913      "ENSGALG00000000003"    "PANX2"         chr1
		@fields = ();
		@fields = split(/\t/, $line);

		$asId = $fields[0];

		next if($asId eq "ASID");
		$symbol = $fields[2];
		if($symbol=~/"(.*)"/){
			$symbol = $1;
		}
		
		$sql = "update asTable set geneSymbol=\"" . $symbol . "\" where asId=\"" . $asId . "\"";
		$update = $dbh->prepare($sql);
		$update->execute();
	}	
	close FF;
}


