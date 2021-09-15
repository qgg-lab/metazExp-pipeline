#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
                "--asFileList 9031.A3ss,9796.A3ss,9823.A3ss,9913.A3ss,9940.A3ss\\\n" .
                "--taxonIdList  9031,9796,9823,9913,9940\\\n" .
		"--asType A3SS \\\n" .
		"--orthAsTag ORTH \\\n" .
                "--outputOrthAsFile  orth.A3ss.tsv \n";
        exit;
}

my ($asFileList, $taxonIdList, $asType, $orthAsTag, $outputOrthAsFile);

GetOptions(
        'asFileList=s'=>\$asFileList,
        'taxonIdList=s'=>\$taxonIdList,
	'asType=s'=>\$asType,
	'orthAsTag=s'=>\$orthAsTag,
        'outputOrthAsFile=s'=>\$outputOrthAsFile,
);

my (@asFile, $asFile, @taxonId, $taxonId, @fields, $line, $orthFlag, $orthAsId);
my (%orthAs, @orthAsId, $orthAsId, $orthAsNum, $standOrthAsId);

@asFile = split(/,/, $asFileList);
@taxonId = split(/,/, $taxonIdList);
# read all as into hash to build whole as catalogo
# 将所有物种的AS读入hash，记录全部位置都能找到orth位置的AS，用其位置组成hash的关键字
# 将所有物种对应的AS编号放置到hash中 
for(my $i=0; $i<=$#asFile; $i++){
	$asFile = $asFile[$i];
	$taxonId = $taxonId[$i];

	open FF, "<$asFile";
	while($line=<FF>){
		chomp($line);
		@fields = ();
		@fields = split(/\t/, $line);
		next if($fields[0] eq "ASID");
		$orthFlag = 1;
		$orthAsId = "";

		if(uc($asType) eq "A5SS"){
			if($fields[4] eq "+"){
				if($fields[8] eq "-" or $fields[6] eq "-" or $fields[9] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[8] . "#" . $fields[6] . "#" . $fields[9];
				}
			}elsif($fields[4] eq "-"){
				if($fields[7] eq "-" or $fields[5] eq "-" or $fields[10] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[7] . "#" . $fields[5] . "#" . $fields[10];
				}
			}
		}elsif(uc($asType) eq "A3SS"){
			if($fields[4] eq "+"){
				if($fields[10] eq "-" or $fields[5] eq "-" or $fields[7] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[10] . "#" . $fields[5] . "#" . $fields[7];
				}
			}elsif($fields[4] eq "-"){
				if($fields[9] eq "-" or $fields[6] eq "-" or $fields[8] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[9] . "#" . $fields[6] . "#" . $fields[8];
				}
			}
		}elsif(uc($asType) eq "SE"){
			if($fields[4] eq "+"){
				if($fields[8] eq "-" or $fields[5] eq "-" or $fields[6] eq "-" or $fields[9] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[8] . "#" . $fields[5] . "#" . $fields[6] . "#" . $fields[9];
				}
			}elsif($fields[4] eq "-"){
				if($fields[9] eq "-" or $fields[6] eq "-" or $fields[5] eq "-" or $fields[8] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[9] . "#" . $fields[6] . "#" . $fields[5] . "#" . $fields[8];
				}
			}
		}elsif(uc($asType) eq "RI"){
			if($fields[4] eq "+"){
				if($fields[8] eq "-" or $fields[9] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[8] . "#" . $fields[9];
				}
			}elsif($fields[4] eq "-"){
				if($fields[9] eq "-" or $fields[8] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[9] . "#" . $fields[8];
				}		
			}
		}elsif(uc($asType) eq "MXE"){
			if($fields[4] eq "+"){
				if($fields[10] eq "-" or $fields[5] eq "-" or $fields[6] eq "-" or $fields[7] eq "-" or $fields[8] eq "-" or $fields[11] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[10] . "#" . $fields[5] . "#" . $fields[6] . "#" . $fields[7] . "#" . $fields[8] . "#" . $fields[11];
				}
			}elsif($fields[4] eq "-"){
				if($fields[11] eq "-" or $fields[8] eq "-" or $fields[7] eq "-" or $fields[6] eq "-" or $fields[5] eq "-" or $fields[10] eq "-"){
					$orthFlag = 0;
				}else{
					$orthAsId = $fields[11] . "#" . $fields[8] . "#" . $fields[7] . "#" . $fields[6] . "#" . $fields[5] . "#" . $fields[10];
				}
			}
		}
		if($orthFlag == 1){
			${$orthAs{$orthAsId}}{$taxonId} = $fields[0];
		}
	}
	close FF;		
}

# 输出orthAS矩阵
$orthAsNum = 0;
@orthAsId = keys(%orthAs);
open WW, ">$outputOrthAsFile";
print WW "orthASId";
foreach $taxonId(@taxonId){
        print WW "\t" . $taxonId;
}
print WW  "\n";

foreach $orthAsId(@orthAsId){
	$orthAsNum++;
	$standOrthAsId = $orthAsTag . $asType . sprintf("%010d", $orthAsNum);
	print WW $standOrthAsId;
	foreach $taxonId(@taxonId){
		if(exists(${$orthAs{$orthAsId}}{$taxonId})){
			print WW "\t" . ${$orthAs{$orthAsId}}{$taxonId};
		}else{
			print WW "\t-";
		}
	}
	print WW "\n";
}
close WW;
