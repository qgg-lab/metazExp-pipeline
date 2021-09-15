#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "perl $0 \\\n" .
	"--tmpDir tmpDir \\\n" .
	"--emfLinkFile first_54_amniotes.pecan.file.link.txt \\\n" .
	"--speciesList gallus_gallus,bos_taurus,equus_caballus,sus_scrofa,ovis_aries \\\n" .
	"--outputOrthRegFile first_54_amniotes.pecan.region.txt \\\n" .
	"--outputOrthSpeDistFile first_54_amniotes.pecan.dis.txt \n\n";
	exit;
}

my ($emfLinkFile, $outputOrthRegFile, $outputOrthSpeDistFile, $speciesList, $tmpDir);
GetOptions(
        'emfLinkFile=s'=>\$emfLinkFile,
	'tmpDir=s'=>\$tmpDir,
	'speciesList=s'=>\$speciesList,
        'outputOrthRegFile=s'=>\$outputOrthRegFile,
        'outputOrthSpeDistFile=s'=>\$outputOrthSpeDistFile,
);

system("mkdir -p " . $tmpDir);
my (%species, @species, $species);
@species = split(/,/, $speciesList);
foreach $species(@species){
	${$species{$species}}{"flag"}=1;
	${$species{$species}}{"start"}=1;
	${$species{$species}}{"stop"}=1;
	${$species{$species}}{"emerge"} = 0;
}


my ($cmd, $link, $line, @fields, $contactSpecies, $orthId);

open WREGION, ">" . $outputOrthRegFile;
print WREGION join("\t", "orthId", "species", "chromosome", "chain", "start", "stop") . "\n";

open WDIST, ">" . $outputOrthSpeDistFile;
foreach $species(@species){
	$contactSpecies .= $species . "\t";
}
$contactSpecies = "orthId\t" . substr($contactSpecies, 0, length($contactSpecies)-1);
print WDIST $contactSpecies . "\n";

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

		chomp($line);
		@fields = ();
		@fields = split(/ /, $line);

		if($fields[0] eq "SEQ"){

			if(${$species{$fields[1]}}{"flag"} == 1){
				${$species{$fields[1]}}{"chr"} = $fields[2];
				${$species{$fields[1]}}{"chain"} = "+";
				${$species{$fields[1]}}{"chain"} = "-" if($fields[5] eq "-1");
				${$species{$fields[1]}}{"start"} = $fields[3];
				${$species{$fields[1]}}{"stop"} = $fields[4];
				${$species{$fields[1]}}{"emerge"} = 1;
			}else{
				#nothing to do
			}

		}elsif($fields[0] eq "ID"){

			# output region of all species
			$orthId = $fields[1];
			foreach $species(@species){
				print WREGION join("\t", $orthId, $species, ${$species{$species}}{"chr"}, ${$species{$species}}{"chain"}, ${$species{$species}}{"start"}, ${$species{$species}}{"stop"}) . "\n" if(${$species{$species}}{"start"} - ${$species{$species}}{"stop"} !=0 and ${$species{$species}}{"emerge"} == 1);
			}

			# output species dist
			print WDIST $orthId;
			foreach $species(@species){
				if(${$species{$species}}{"start"} - ${$species{$species}}{"stop"} !=0 and ${$species{$species}}{"emerge"} == 1){
					print WDIST "\t1";
				}else{
					print WDIST "\t0";
				}
			}
			print WDIST "\n";
			
			# clear species status
			foreach $species(@species){
				${$species{$species}}{"emerge"} = 0;
			}
		}
	}
	close EMF;
}
close FF;


