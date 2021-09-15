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
		"--totalSampleInfoFile sample.information.tsv \\\n" .
		"--filterSampleInfoFile samples.filtered.by.sequencing.and.mapping.tsv \n";
	exit;
}

my ($dbName, $dbUser, $dbPWD);
my ($species);
my ($totalSampleInfoFile, $filterSampleInfoFile);
my (%expToStudy, @fields);

GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'species=s'=>\$species,
	'totalSampleInfoFile=s'=>\$totalSampleInfoFile,
	'filterSampleInfoFile=s'=>\$filterSampleInfoFile
);

my $dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
my $query;

# load experiment => studyId into hash
my ($line);
open FF, "<$totalSampleInfoFile";
while($line=<FF>){
	@fields = ();
	@fields = split(/\|___\|/, $line);
	$expToStudy{$fields[18]}=$fields[57];
}
close FF;

my ($sql);
my ($expId, $status, $runNum, $runId, $library, $layout, $phredScore, $readLength, $spotNum, $alignPer, $novelSpliceNum, $totalA5SS, $totalA3SS, $totalSE, $totalRI, $totalMXE,$novelA5SS, $novelA3SS, $novelSE, $novelRI, $novelMXE, $jcecA5SS, $jcecA3SS, $jcecSE, $jcecRI, $jcecMXE, $jcA5SS, $jcA3SS, $jcSE, $jcRI, $jcMXE);
my ($jcecTotalAsNum, $jcecA5ssPer, $jcecA3ssPer, $jcecSePer, $jcecRiPer, $jcecMxePer);
my ($jcTotalAsNum, $jcA5ssPer, $jcA3ssPer, $jcSePer, $jcRiPer, $jcMxePer);
my ($mappedSpots);

open FF, "<$filterSampleInfoFile";
#expId   status  runNum  runId   library layout  phredScore      readLength      spotNum(M)      alignPer(%)     novelSpliceNum  totalA5SS       totalA3SS       totalSE totalRI totalMXE        novelA5SS       novelA3SS       novelSE novelRI novelMXE        jcecA5SS        jcecA3SS        jcecSE  jcecRI  jcecMXE jcA5SS  jcA3SS  jcSE    jcRI    jcMXE
<FF>;
while($line=<FF>){
	chomp($line);

	($expId, $status, $runNum, $runId, $library, $layout, $phredScore, $readLength, $spotNum, $alignPer, $novelSpliceNum, $totalA5SS, $totalA3SS, $totalSE, $totalRI, $totalMXE,$novelA5SS, $novelA3SS, $novelSE, $novelRI, $novelMXE, $jcecA5SS, $jcecA3SS, $jcecSE, $jcecRI, $jcecMXE, $jcA5SS, $jcA3SS, $jcSE, $jcRI, $jcMXE) = split(/\t/, $line);
	
	$alignPer = substr($alignPer, 0, length($alignPer) -1);
	$mappedSpots = int(($spotNum * $alignPer/100)*100)/100;

#	print join("\t", $alignPer, $mappedSpots, $expId, $status, $runNum, $runId, $library, $layout, $phredScore, $readLength, $spotNum, $alignPer, $novelSpliceNum, $totalA5SS, $totalA3SS, $totalSE, $totalRI, $totalMXE,$novelA5SS, $novelA3SS, $novelSE, $novelRI, $novelMXE, $jcecA5SS, $jcecA3SS, $jcecSE, $jcecRI, $jcecMXE, $jcA5SS, $jcA3SS, $jcSE, $jcRI, $jcMXE) . "\n";
#	<STDIN>;
	
	$jcecTotalAsNum = $jcecA5SS + $jcecA3SS + $jcecSE + $jcecRI + $jcecMXE;
	$jcecA5ssPer = int(($jcecA5SS/$jcecTotalAsNum)*10000)/100;
	$jcecA3ssPer = int(($jcecA3SS/$jcecTotalAsNum)*10000)/100;
	$jcecSePer = int(($jcecSE/$jcecTotalAsNum)*10000)/100;
	$jcecRiPer = int(($jcecRI/$jcecTotalAsNum)*10000)/100;
	$jcecMxePer = int(($jcecMXE/$jcecTotalAsNum)*10000)/100;

	$jcTotalAsNum = $jcA5SS + $jcA3SS + $jcSE + $jcRI + $jcMXE;
	$jcA5ssPer = int(($jcA5SS/$jcTotalAsNum)*10000)/100;
	$jcA3ssPer = int(($jcA3SS/$jcTotalAsNum)*10000)/100;
	$jcSePer = int(($jcSE/$jcTotalAsNum)*10000)/100;
	$jcRiPer = int(($jcRI/$jcTotalAsNum)*10000)/100;
	$jcMxePer = int(($jcMXE/$jcTotalAsNum)*10000)/100;

	$sql = "insert into experimentTable (alignPercent, mappedSpots, experimentId, species, libraryType, libraryLayout, readLen, phredScore, totalSpots, runIdList, runNum, studyId, jcecTotalAsNum, jcecA5ssPercentage, jcecA3ssPercentage, jcecSePercentage, jcecRiPercentage, jcecMxePercentage, jcTotalAsNum, jcA5ssPercentage, jcA3ssPercentage, jcSePercentage, jcRiPercentage, jcMxePercentage) values (" . $alignPer . ", " . $mappedSpots  . ", \"" . $expId . "\", \"" . $species . "\", \"" . $library . "\", \"" . $layout . "\", \"" . $readLength . "\", \"" . $phredScore . "\", " . $spotNum . ", \"" . $runId . "\", " . $runNum . ", \"" . $expToStudy{$expId} . "\", " . $jcecTotalAsNum . ", " . $jcecA5ssPer . ", " . $jcecA3ssPer . ", " . $jcecSePer  . ", " . $jcecRiPer . ", " . $jcecMxePer . "," . $jcTotalAsNum . ", " . $jcA5ssPer . ", " . $jcA3ssPer . ", " . $jcSePer  . ", " . $jcRiPer . ", " . $jcMxePer . ")";
#	print $sql;
#	<STDIN>;
	$query = $dbh->prepare($sql);
	$query->execute();
}
close FF;

