#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--species \"Bos taurus\" \\\n" .
		"--sampleInfoFile sample.information.tsv \\\n" .
		"--inputJCECfileList  jcec.A3SS.tsv,jcec.A5SS.tsv,jcec.MXE.tsv,jcec.RI.tsv,jcec.SE.tsv \\\n" . 
		"--inputJCfileList jc.A3SS.tsv,jc.A5SS.tsv,jc.MXE.tsv,jc.RI.tsv,jc.SE.tsv \\\n" .
		"--psiTsvFile 89462.psi.tsv \n\n";
	exit;
}

my ($species);
my ($inputJCECfileList, $inputJCfileList, $sampleInfoFile, $psiTsvFile);

GetOptions(
        'species=s'=>\$species,
	'sampleInfoFile=s'=>\$sampleInfoFile,
	'inputJCECfileList=s'=>\$inputJCECfileList, 
	'inputJCfileList=s'=>\$inputJCfileList,
	'psiTsvFile=s'=>\$psiTsvFile,
);

my ($line, @fields, %as, @asId, $asId, @experimentId, $experimentId, $sql, %experiment, $jcecPsi, $jcPsi);
# load experiment => studyId into hash
open FF, "<$sampleInfoFile";
while($line=<FF>){
	@fields = ();
	@fields = split(/\|___\|/, $line);
	$experiment{$fields[18]}=$fields[57];
}
close FF;

# jcec
my @fileName = split(/,/, $inputJCECfileList);
foreach my $fileName(@fileName){
	open FF, "<$fileName";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@fields = split(/\t/, $line);
		#      asId         experiment
		if($fields[1] > 0 or $fields[2] > 0){
			${$as{$fields[0]}}{$fields[5]} = join("\t", $fields[1], $fields[2], $fields[3], $fields[4]);
		}
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
		if($fields[1] > 0 or $fields[2] > 0){
			if(exists(${$as{$fields[0]}}{$fields[5]})){
				${$as{$fields[0]}}{$fields[5]} = ${$as{$fields[0]}}{$fields[5]} . "\t" . join("\t", $fields[1], $fields[2], $fields[3], $fields[4]);
			}else{
				${$as{$fields[0]}}{$fields[5]} = join("\t", "0", "0", "0", "0") . "\t" . join("\t", $fields[1], $fields[2], $fields[3], $fields[4]);
			}
		}
	}
	close FF;
}

open WW, ">$psiTsvFile";
@asId = keys(%as);
foreach $asId(@asId){
	@experimentId = ();
	@experimentId = keys(%{$as{$asId}});
	foreach $experimentId(@experimentId){

		@fields = ();
		@fields = split(/\t/, ${$as{$asId}}{$experimentId});

		if($fields[0] == 0 and $fields[1] != 0){
			$jcecPsi = 0;
		}elsif($fields[0] == 0 and $fields[1] == 0){
			$jcecPsi = -1;
		}elsif($fields[0] != 0 and $fields[1] == 0){
			$jcecPsi = 1;
		}else{
			$jcecPsi = ($fields[0]/$fields[2])/($fields[0]/$fields[2] + $fields[1]/$fields[3]);
		}

		if($fields[4] == 0 and $fields[5] != 0){
			$jcPsi = 0;
		}elsif($fields[4] == 0 and $fields[5] == 0){
			$jcPsi = -1;
		}elsif($fields[4] != 0 and $fields[5] == 0){
			$jcPsi = 1;
		}else{
			$jcPsi = ($fields[4]/$fields[6])/($fields[4]/$fields[6] + $fields[5]/$fields[7]);
		}

		print WW "asId, species, experiment, studyId, JCECI, JCECS, JCECIncFormLen, JCECSkipFormLen, JCECpsi, JCI, JCS, JCIncFormLen, JCSkipFormLen, JCpsi" . "_____"  . "\"" . $asId . "\", \"" . $species . "\", \"" . $experimentId . "\", \"" . $experiment{$experimentId} . "\", " . $fields[0] . ", " . $fields[1] . ", " . $fields[2] . ", " . $fields[3] . ", " . $jcecPsi . ", " . $fields[4] . ", " . $fields[5] . ", " . $fields[6] . ", " . $fields[7] . ", " . $jcPsi . "\n";
	}
}
close WW;
