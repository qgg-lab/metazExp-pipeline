#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
                "--exonFileList 9031.exon.align.coord.tsv,9796.exon.align.coord.tsv,9823.exon.align.coord.tsv,9913.exon.align.coord.tsv,9940.exon.align.coord.tsv\\\n" .
                "--taxonIdList  9031,9796,9823,9913,9940\\\n" .
		"--orthExonTag ORTH \\\n" .
                "--outputOrthExonFile  orth.exon.tsv \n";
        exit;
}

my ($exonFileList, $taxonIdList, $orthExonTag, $outputOrthExonFile);

GetOptions(
        'exonFileList=s'=>\$exonFileList,
        'taxonIdList=s'=>\$taxonIdList,
	'orthExonTag=s'=>\$orthExonTag,
        'outputOrthExonFile=s'=>\$outputOrthExonFile,
);

my (@exonFile, $exonFile, @taxonId, $taxonId, @fields, $line, $orthFlag, $orthExonId);
my (%orthExon, @orthExonId, $orthExonId, $orthExonNum, $standOrthExonId);

@exonFile = split(/,/, $exonFileList);
@taxonId = split(/,/, $taxonIdList);
# read all exon into hash to build whole as catalogo
# 将所有物种的Exon读入hash，记录全部位置都能找到orth位置的Exon，用其位置组成hash的关键字
# 将所有物种对应的AS编号放置到hash中 
for(my $i=0; $i<=$#exonFile; $i++){
	$exonFile = $exonFile[$i];
	$taxonId = $taxonId[$i];
	open FF, "<$exonFile";
	while($line=<FF>){
		chomp($line);
		@fields = ();
		@fields = split(/\t/, $line);

		$orthFlag = 1;
		$orthExonId = "";
		for(my $j=1; $j<=$#fields; $j++){
			$orthFlag = 0 if($fields[$j] eq "-");
			$orthExonId .= $fields[$j] . "#";
		}

		if($orthFlag == 1){
			${$orthExon{$orthExonId}}{$taxonId} = $fields[0];
		}
	}
	close FF;		
}

# 输出orthExon矩阵
$orthExonNum = 0;
@orthExonId = keys(%orthExon);
open WW, ">$outputOrthExonFile";
print WW "orthExonId";
foreach $taxonId(@taxonId){
	print WW "\t" . $taxonId;
}
print WW  "\n";
foreach $orthExonId(@orthExonId){
	$orthExonNum++;
	$standOrthExonId = $orthExonTag . sprintf("%010d", $orthExonNum);
	print WW $standOrthExonId;
	foreach $taxonId(@taxonId){
		if(exists(${$orthExon{$orthExonId}}{$taxonId})){
			print WW "\t" . ${$orthExon{$orthExonId}}{$taxonId};
		}else{
			print WW "\t-";
		}
	}
	print WW "\n";
}
close WW;
