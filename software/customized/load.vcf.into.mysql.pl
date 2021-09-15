#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--dbName asdb \\\n" .
                "--dbUser lsas \\\n" .
                "--dbPWD njaulsas2019 \\\n" .
		"--species \"Equus caballus\" \\\n" .
		"--genomeFile  genome.fa \\\n" . 
		"--inputVcfFile genome.vcf \n";
	exit;
}

my ($genomeFile, $inputVcfFile, $outputVcfFile, $species);
my ($dbName, $dbUser, $dbPWD);
GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'genomeFile=s'=>\$genomeFile,
	'species=s'=>\$species,
        'inputVcfFile=s'=>\$inputVcfFile,
);

# read genome sequence into hash
my (%genomeSeq, $line, @tt, $id);
open FF, "<$genomeFile";
while($line = <FF>){
        chomp($line);
        if($line=~/>/){
                @tt = ();
                @tt = split(/ /, $line);
                if($tt[0]=~/>ref\|(.*?)\|/ or $tt[0]=~/>gi.*\|ref\|(.*?)\|/ or  $tt[0]=~/>(.*)/){
                        $id = $1;
                }
        }else{
                $genomeSeq{$id} .=uc($line);
        }
}
close FF;

print "Finish load genome sequence into hash.\n";

my $dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
my $query;
my $sql = "select chr, start, end from asTable where species=\"$species\"";
$query = $dbh->prepare($sql);
$query->execute();
my $sblen;
while(my $row=$query->fetchrow_hashref()){
	$sblen = $row->{"end"} - $row->{"start"} + 1;
	substr($genomeSeq{$row->{"chr"}}, $row->{"start"}-1, $sblen) = "Z"x$sblen;
}
$query->finish();
$dbh->disconnect();

print "Finish mask as region.\n";


$dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
# filter vcf position
my ($chr, $pos, $ref, $alt, $accId, $varType);
open FF, "<" . $inputVcfFile;
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#1       207381  rs719849645     C       T       .       .       dbSNP_150;TSA=SNV
while($line=<FF>){
	next if($line=~/#/);
	@tt = ();
	@tt = split(/\t/, $line);
	($chr, $pos, $accId, $ref, $alt, $varType) = ($tt[0], $tt[1], $tt[2], $tt[3], $tt[4], $tt[7]);
	if($varType=~/=SNV/){
		$varType="SNV";
	}elsif($varType=~/=deletion/){
		$varType="deletion";
	}elsif($varType=~/=insertion/){
		$varType="insertion";
	}else{
		next;	
	}

	if(substr($genomeSeq{$tt[0]}, $tt[1]-1, 1) eq "Z"){
		$sql = "insert into variantTable (species, chr, pos, ref, alt, accId, varType) values (\"" . $species . "\", \"" . $chr . "\", $pos, \"$ref\", \"$alt\", \"$accId\", \"$varType\")";

		$query = $dbh->prepare($sql);
		$query->execute();
	}
}
close FF;
$dbh->disconnect();
