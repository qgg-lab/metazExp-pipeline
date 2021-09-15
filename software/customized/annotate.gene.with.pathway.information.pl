#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--koAnnoTsv \\\n" .
                "--keggDicTsv \\\n" .
		"--outputCompleteKeggAnnoMysqlTsv \\\n" .
		"--outputCompleteKeggAnnoTsv \n";
	exit;
}

my ($koAnnoTsv, $keggDicTsv, $outputCompleteKeggAnnoTsv, $outputCompleteKeggAnnoMysqlTsv, $line);

GetOptions(
        'koAnnoTsv=s'=>\$koAnnoTsv,
        'keggDicTsv=s'=>\$keggDicTsv,
	'outputKeggAnnoMysqlTsv=s'=>\$outputCompleteKeggAnnoMysqlTsv,
        'outputKeggAnnoTsv=s'=>\$outputCompleteKeggAnnoTsv,
);

my (%orthHash, $orthHashHref);
my ($Class1Id, $Class1Name, $Class2Id, $Class2Name, $pathwayId, $pathwayName, $orthology, $enzymeAbbr, $enzymeName, $enzymeId);
$orthHashHref=\%orthHash;
open FF, "<$keggDicTsv";
<FF>;
while($line=<FF>){
	chomp($line);
	($Class1Id, $Class1Name, $Class2Id, $Class2Name, $pathwayId, $pathwayName, $orthology, $enzymeAbbr, $enzymeName, $enzymeId) = split(/\t/, $line);
	$orthHashHref->{$orthology}->{"Class1Id"} = $Class1Id;
	$orthHashHref->{$orthology}->{"Class1Name"} = $Class1Name;
	$orthHashHref->{$orthology}->{"Class2Id"} = $Class2Id;
	$orthHashHref->{$orthology}->{"Class2Name"} = $Class2Name;
	$orthHashHref->{$orthology}->{"pathwayId"} = $pathwayId;
	$orthHashHref->{$orthology}->{"pathwayName"} = $pathwayName;
	$orthHashHref->{$orthology}->{"enzymeAbbr"} = $enzymeAbbr;
	$orthHashHref->{$orthology}->{"enzymeName"} = $enzymeName;
	$orthHashHref->{$orthology}->{"enzymeId"} = $enzymeId;
}
close FF;

my ($geneId);
open FF, "<$koAnnoTsv";
open MYSQL, ">$outputCompleteKeggAnnoMysqlTsv";
open WW, ">$outputCompleteKeggAnnoTsv";
print WW join("\t", "geneId", "orthology", "1stClassId", "Class1Name", "Class2Id", "Class2Name", "pathwayId", "pathwayName", "enzymeAbbr", "enzymeName", "enzymeId") . "\n";
while($line=<FF>){
	chomp($line);
	next if(not($line=~/\t/));
	($geneId, $orthology) = split(/\t/, $line);
	print WW join("\t", $geneId, $orthology, $orthHashHref->{$orthology}->{"Class1Id"}, $orthHashHref->{$orthology}->{"Class1Name"}, $orthHashHref->{$orthology}->{"Class2Id"}, $orthHashHref->{$orthology}->{"Class2Name"}, $orthHashHref->{$orthology}->{"pathwayId"}, $orthHashHref->{$orthology}->{"pathwayName"}, $orthHashHref->{$orthology}->{"enzymeAbbr"}, $orthHashHref->{$orthology}->{"enzymeName"}, $orthHashHref->{$orthology}->{"enzymeId"}) . "\n";
	print MYSQL join("\t", "geneId", "orthology", "1stClassId", "Class1Name", "Class2Id", "Class2Name", "pathwayId", "pathwayName", "enzymeAbbr", "enzymeName", "enzymeId") . "___" . join("\t", $geneId, $orthology, $orthHashHref->{$orthology}->{"Class1Id"}, $orthHashHref->{$orthology}->{"Class1Name"}, $orthHashHref->{$orthology}->{"Class2Id"}, $orthHashHref->{$orthology}->{"Class2Name"}, $orthHashHref->{$orthology}->{"pathwayId"}, $orthHashHref->{$orthology}->{"pathwayName"}, $orthHashHref->{$orthology}->{"enzymeAbbr"}, $orthHashHref->{$orthology}->{"enzymeName"}, $orthHashHref->{$orthology}->{"enzymeId"}) . "\n";
}
close FF;
close WW;
close MYSQL;
