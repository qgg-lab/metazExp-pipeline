#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--dbName asdb \\\n" . 
		"--dbUser lsas \\\n" . 
		"--dbPWD njaulsas2019 \\\n" . 
		"--species \"Bos taurus\" \\\n" .
		"--inputASfile  A5SS.catalog \\\n" . 
		"--inputAStype  A5SS \n";
	exit;
}

my ($dbName, $dbUser, $dbPWD);
my ($species);
my ($inputASfile, $inputAStype);

GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'species=s'=>\$species,
	'inputASfile=s'=>\$inputASfile, 
	'inputAStype=s'=>\$inputAStype, 
);

my ($line, @fields, $tableFields, $values, $sql);
#my $dbh = DBI->connect("DBI:mysql:database=$dbName", $dbUser, $dbPWD);

my $dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
my $query;
# ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
open FF, "<$inputASfile";
	$line = <FF>;
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	shift(@fields);
	shift(@fields);
	shift(@fields);

$tableFields = "asId, asType, species, geneId, geneSymbol, " . join(",", @fields);

while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	if($fields[2] eq "NA"){
		$values = "\"" . $fields[0] . "\", \"" . $inputAStype . "\", \"" . $species . "\", " . $fields[1] . ", \"" . $fields[2] . "\", \"" . substr($fields[3], 3) . "\", \"" . $fields[4] . "\", ";
	}else{		
		$values = "\"" . $fields[0] . "\", \"" . $inputAStype . "\", \"" . $species . "\", " . $fields[1] . ", " . $fields[2] . ", \"" . substr($fields[3], 3) . "\", \"" . $fields[4] . "\", ";
	}
	for(my $i=5; $i<=$#fields; $i++){
		$values .=$fields[$i] . ",";
	}
	$values = substr($values, 0, length($values)-1);
	$sql = "insert into asTable (" . $tableFields . ") values(" . $values . ")";
	$query = $dbh->prepare($sql);
	$query->execute();
}

close FF;
