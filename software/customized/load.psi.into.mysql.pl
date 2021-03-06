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
		"--sampleInfoFile sample.information.tsv \\\n" .
		"--inputJCECfileList  jcec.A3SS.tsv,jcec.A5SS.tsv,jcec.MXE.tsv,jcec.RI.tsv,jcec.SE.tsv \\\n" . 
		"--inputJCfileList jc.A3SS.tsv,jc.A5SS.tsv,jc.MXE.tsv,jc.RI.tsv,jc.SE.tsv \n";
	exit;
}

my ($dbName, $dbUser, $dbPWD);
my ($species);
my ($inputJCECfileList, $inputJCfileList, $sampleInfoFile);

GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'species=s'=>\$species,
	'sampleInfoFile=s'=>\$sampleInfoFile,
	'inputJCECfileList=s'=>\$inputJCECfileList, 
	'inputJCfileList=s'=>\$inputJCfileList, 
);

my ($line, @fields, %as, @asId, $asId, @experimentId, $experimentId, $sql, %experiment, $jcecPsi, $jcPsi);
#my $dbh = DBI->connect("DBI:mysql:database=$dbName", $dbUser, $dbPWD);


# load experiment => studyId into hash
open FF, "<$sampleInfoFile";
while($line=<FF>){
	@fields = ();
	@fields = split(/\|___\|/, $line);
	$experiment{$fields[18]}=$fields[57];
}
close FF;


my $dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
my $query;

# as->exp1 .. exp3 
#     expi -> ijcec, sjcec, incformlenjcec, skipformlenjcec, ijc, sjc, incformlenjc, skipformlenjc

# jcec
my @fileName = split(/,/, $inputJCECfileList);
foreach my $fileName(@fileName){
	open FF, "<$fileName";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@fields = split(/\t/, $line);
		#      asId         experiment
		${$as{$fields[0]}}{$fields[5]} = join("\t", $fields[1], $fields[2], $fields[3], $fields[4]);
	}
	close FF;
}

# jc
my @fileName = split(/,/, $inputJCfileList);
foreach my $fileName(@fileName){
	open FF, "<$fileName";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@fields = split(/\t/, $line);
		if(exists(${$as{$fields[0]}}{$fields[5]})){
			${$as{$fields[0]}}{$fields[5]} = ${$as{$fields[0]}}{$fields[5]} . "\t" . join("\t", $fields[1], $fields[2], $fields[3], $fields[4]);
		}else{
			${$as{$fields[0]}}{$fields[5]} = join("\t", "0", "0", "0", "0") . "\t" . join("\t", $fields[1], $fields[2], $fields[3], $fields[4]);
		}
	}
	close FF;
}

@asId = keys(%as);
foreach $asId(@asId){
	@experimentId = ();
	@experimentId = keys(%{$as{$asId}});
	foreach $experimentId(@experimentId){

		@fields = ();
		@fields = split(/\t/, ${$as{$asId}}{$experimentId});

		if($fields[0] == 0 and $fields[1] == 0){
			$jcecPsi = 0;
		}else{
			$jcecPsi = ($fields[0]/$fields[2])/($fields[0]/$fields[2] + $fields[1]/$fields[3]);
		}

		if($fields[4] == 0 and $fields[5] == 0){
			$jcPsi = 0;
		}else{
			$jcPsi = ($fields[4]/$fields[6])/($fields[4]/$fields[6] + $fields[5]/$fields[7]);
		}

		$sql = "insert into psiTable (asId, species, experiment, studyId, JCECI, JCECS, JCECIncFormLen, JCECSkipFormLen, JCECpsi, JCI, JCS, JCIncFormLen, JCSkipFormLen, JCpsi) values(\"" . $asId . "\", \"" . $species . "\", \"" . $experimentId . "\", \"" . $experiment{$experimentId} . "\", " . $fields[0] . ", " . $fields[1] . ", " . $fields[2] . ", " . $fields[3] . ", " . $jcecPsi . ", " . $fields[4] . ", " . $fields[5] . ", " . $fields[6] . ", " . $fields[7] . ", " . $jcPsi . ")";

		$query = $dbh->prepare($sql);
		$query->execute();
	}
}
