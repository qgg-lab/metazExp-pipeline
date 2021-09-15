#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--keggJson \\\n" .
                "--keggMysqlTsv \n";
	exit;
}

my ($keggJson, $keggMysqlTsv);

GetOptions(
        'keggJson=s'=>\$keggJson,
        'keggMysqlTsv=s'=>\$keggMysqlTsv,
);

my ($line, $firstLevelName, $firstLevelId, $secondLevelName, $secondLevelId, $thirdLevelName, $thirdLevelId, $thirdLevelType, $thirdPathwayId);
my ($orthId, $enzymeAbbr, $enzymeName, $enzymeId);
open WW, ">$keggMysqlTsv";
open FF, "<$keggJson";
while($line=<FF>){
	if($line=~/^\t\t"name":"(\d+) (.*)",\n/){
		$firstLevelId = $1;
		$firstLevelName = $2;
		#print join("\t", $firstLevelId, $firstLevelName) . "\n";
		#<STDIN>;
	}elsif($line=~/^\t\t\t"name":"(\d+) (.*)",\n/){
		$secondLevelId = $1;
		$secondLevelName = $2;
		#print join("\t", $secondLevelId, $secondLevelName) . "\n";
		#<STDIN>;
	}elsif($line=~/^\t\t\t\t"name":"(\d+) (.*) \[(.*):(ko\d+)\]",/){
		$thirdLevelId = $1;
		$thirdLevelName = $2;
		$thirdPathwayId = $4;
		$thirdLevelType = $3;
		#print join("\t", $thirdLevelId, $thirdLevelName);
		#<STDIN>;
	}elsif($line=~/^\t\t\t\t"name":"(\d+) (.*)",/){
		$thirdLevelId = $1;
		$thirdLevelName = $2;
		$thirdPathwayId = "none";
		$thirdLevelType = "none";
	}elsif($line=~/^\t\t\t\t\t"name":"(K\d+)  .*"\n/){
		if($line=~/^\t\t\t\t\t"name":"(K\d+)  (.*); (.*) \[(.*)\]"\n/){
			$orthId = $1;
			$enzymeAbbr = $2;
			$enzymeName = $3;
			$enzymeId = $4;
			print WW join("\t", "1stClassId", "1stClassName", "2ndClassId", "2ndClassName", "3rdClassId", "3rdClassName", "3rdClassType", "pathwayId", "orthology", "enzymeAbbr", "enzymeName", "enzymeId") . "___" . join("\t", $firstLevelId, $firstLevelName, $secondLevelId, $secondLevelName, $thirdLevelId, $thirdLevelName, $thirdLevelType, $thirdPathwayId, $orthId, $enzymeAbbr, $enzymeName, $enzymeId) . "\n";
		}elsif($line=~/^\t\t\t\t\t"name":"(K\d+)  (.*); (.*)"\n/){
			$orthId = $1;
			$enzymeAbbr = $2;
			$enzymeName = $3;
			print WW join("\t", "1stClassId", "1stClassName", "2ndClassId", "2ndClassName", "3rdClassId", "3rdClassName", "3rdClassType", "pathwayId", "orthology", "enzymeAbbr", "enzymeName", "enzymeId") . "___" . join("\t", $firstLevelId, $firstLevelName, $secondLevelId, $secondLevelName, $thirdLevelId, $thirdLevelName, $thirdLevelType, $thirdPathwayId, $orthId, $enzymeAbbr, $enzymeName, "-") . "\n";
		}
	}
}
close FF;
close WW;
