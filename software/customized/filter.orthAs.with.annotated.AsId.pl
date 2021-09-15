#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--conservMatrix \\\n" .
                "--asIdListFileList \\\n" .
                "--outputFilteredConservMatrix \n";
	exit;
}

my ($conservMatrix, $asIdListFileList, $outputFilteredConservMatrix);

GetOptions(
        'conservMatrix=s'=>\$conservMatrix,
        'asIdListFileList=s'=>\$asIdListFileList,
        'outputFilteredConservMatrix=s'=>\$outputFilteredConservMatrix,
);

# 将所有AS读入hash
my (%asId);
my ( @asIdListFile, $asIdListFile, $asId);
@asIdListFile = split(/,/, $asIdListFileList);
foreach $asIdListFile(@asIdListFile){
	open FF, "<$asIdListFile";
	while($asId=<FF>){
		chomp($asId);
#		print $asId;
#		<STDIN>;
		$asId{$asId} = 1;
	}
	close FF;
}

my (@asId, $checkRlt);
# 将orth矩阵读入，如果任何一个物种中的as不存在于%as，那么放弃该orth
open WW, ">$outputFilteredConservMatrix";
open FF, "<$conservMatrix";
# orthAsId        9031    9796    9823    9913    		9940
# BOVSE00055548   -       -       -       BTAUSE0000189862      OARISE0000046685
# BOVSE00039639   -       -       -       BTAUSE0000089240      OARISE0000064993
<FF>;
while(my $line=<FF>){
	chomp($line);
	@asId = split(/\t/, $line);
	shift(@asId);
	$checkRlt = 1;
	foreach $asId(@asId){
		if($asId ne"-" and not exists($asId{$asId})){
			$checkRlt = 0;
			last;
		}
	}
	if($checkRlt == 1){
		print WW $line . "\n";
	}
}
close FF;
close WW;
