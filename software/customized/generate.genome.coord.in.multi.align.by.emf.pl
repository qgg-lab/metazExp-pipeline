#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "perl $0 \\\n" .
	"--tmpDir tmpDir \\\n" .
	"--emfLinkFile xaa.link.txt \\\n" .
	"--species bos_taurus \\\n" .
	"--outputFile xaa.tsv \n\n";
	exit;
}

my ($emfLinkFile, $outputFile, $species, $tmpDir);
GetOptions(
        'emfLinkFile=s'=>\$emfLinkFile,
	'tmpDir=s'=>\$tmpDir,
	'species=s'=>\$species,
        'outputFile=s'=>\$outputFile,
);

system("mkdir -p " . $tmpDir);

my ($cmd, $link, $line, @fields, $blockId, $speciesId, $blockPos, %genomeCoord);

open WW, ">" . $outputFile;

open FF, "<$emfLinkFile";
while($link=<FF>){
	
	system("rm -rf " . $tmpDir . "/emf.gz");
	system("rm -rf " . $tmpDir . "/emf");

	chomp($link);
	$cmd = "wget " . $link . " -O " . $tmpDir . "/emf.gz " .
		"-o " . $tmpDir . "/log.o.wget.emf.gz " . 
		"2> " . $tmpDir . "/log.e.wget.emf.gz ";

	system($cmd);
	system("gunzip " . $tmpDir . "/emf.gz");

	open EMF, "<" .  $tmpDir . "/emf";
	while($line=<EMF>){

		next if($line=~/#/ or $line=~/\/\// or $line=~/^SCORE/);

		chomp($line);
		@fields = ();
		@fields = split(/ /, $line);

		# begin of a new block
		if($fields[0] eq ""){

			$speciesId = 0;
			$genomeCoord{"speciesId"} = -1;
			$genomeCoord{"emerge"} = 0;
			$genomeCoord{"chain"} = "";
			$genomeCoord{"start"} = 0;
			$genomeCoord{"stop"} = 0;			
			$genomeCoord{"chr"} = "";
			

		}elsif($fields[0] eq "SEQ"){			

			if($fields[1] eq $species){
				$genomeCoord{"speciesId"} = $speciesId;
				$genomeCoord{"chr"} = $fields[2];
				$genomeCoord{"chain"} = "+";
				$genomeCoord{"chain"} = "-" if($fields[5] eq "-1");
				$genomeCoord{"start"} = $fields[3];
				$genomeCoord{"stop"} = $fields[4];
				$genomeCoord{"emerge"} = 1;
				print WW join("\t", $species, $fields[2], $genomeCoord{"chain"}, $fields[3], $fields[4]);
			}
			$speciesId++;

		}elsif($fields[0] eq "ID"){

			$blockId = $fields[1];
			print WW "\t" . $blockId . "\n" if($genomeCoord{"emerge"} == 1);

		}elsif($fields[0] eq "DATA"){

			$blockPos = 0;

		}else{

			$blockPos++;
			if($genomeCoord{"emerge"} == 1 and substr($line, $genomeCoord{"speciesId"}, 1) ne "-"){
				print WW $blockPos . "\n";
			}
		}

	}
	close EMF;
	system("rm -rf " . $tmpDir . "/emf.gz");
	system("rm -rf " . $tmpDir . "/emf");

}
close FF;
close WW;
