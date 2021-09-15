#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--origCdnaFasta \\\n" .
                "--origOrfInCdna \\\n" .
		"--origPfamNameFile \\\n" .
                "--origPfamInCdna \\\n" .
                "--origGoUniqListInCdna \\\n" .
                "--asAlteredCdnaFasta \\\n" .
                "--asAlteredOrfInCdna \\\n" .
		"--asAlteredPfamNameFile \\\n" .
                "--asAlteredPfamInCdna \\\n" .
                "--asAlteredGoUniqListInCdna \\\n" .
                "--outputFunctionAffectionTsv \n";
	exit;
}

my (
	$origCdnaFasta, $origOrfInCdna, $origPfamInCdna, $origGoUniqListInCdna, $origPfamNameFile,
	$asAlteredCdnaFasta, $asAlteredOrfInCdna, $asAlteredPfamInCdna, $asAlteredGoUniqListInCdna, $asAlteredPfamNameFile,
	$outputFunctionAffectionTsv
);

GetOptions(
        'origCdnaFasta=s'=>\$origCdnaFasta,
        'origOrfInCdna=s'=>\$origOrfInCdna,
        'origPfamInCdna=s'=>\$origPfamInCdna,
	'origPfamNameFile=s'=>\$origPfamNameFile,
        'origGoUniqListInCdna=s'=>\$origGoUniqListInCdna,
	'asAlteredCdnaFasta=s'=>\$asAlteredCdnaFasta,
	'asAlteredOrfInCdna=s'=>\$asAlteredOrfInCdna,
	'asAlteredPfamNameFile=s'=>\$asAlteredPfamNameFile,
	'asAlteredPfamInCdna=s'=>\$asAlteredPfamInCdna,
	'asAlteredGoUniqListInCdna=s'=>\$asAlteredGoUniqListInCdna,
	'outputFunctionAffectionTsv=s'=>\$outputFunctionAffectionTsv,
);


my (%asOnTrspt, $asOnTrsptHref);
my ($line, @trsptId, $trsptId, $asId, @field);
# 将origCdnaFasta读入hash
$asOnTrsptHref=\%asOnTrspt;
$asId = "NA";
open FF, "<$origCdnaFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@trsptId = split(/ /, $1);
		$trsptId = $trsptId[0];
	}else{
		$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"seq"} .= $line;
	}
}
close FF;

# 将orf位置读入
my ($orfStart, $orfStop, $startCodonStart, $startCodonStop, $stopCodonStart, $stopCodonStop);
open FF, "<$origOrfInCdna";
#AT3G53400.1     420     1817    420     422     1818    1820
#AT5G45110.2     515     2062    515     517     2063    2065
while($line=<FF>){
	chomp($line);
	($trsptId, $orfStart, $orfStop, $startCodonStart, $startCodonStop, $stopCodonStart, $stopCodonStop) = split(/\t/, $line);
	$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"orfStart"} = $orfStart;
	$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"orfStop"} = $orfStop;
	$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"startCodonStart"} = $startCodonStart;
	$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"startCodonStop"} = $startCodonStop;
	$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"stopCodonStart"} = $stopCodonStart;
	$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"stopCodonStop"} = $stopCodonStop;
}
close FF;

# 将pfam的名称读入
my (%pfamIdToName, $pfamId, $pfamName, $beginInCdna, $endInCdna);
open FF, "<$origPfamNameFile";
# trsptId pfamId  pfamName        beginInCdna     endInCdna
# AT5G45110.2     PF00651 BTB/POZ domain  680     1054
<FF>; 
while($line=<FF>){
	chomp($line);
	($trsptId, $pfamId, $pfamName, $beginInCdna, $endInCdna) = split(/\t/, $line);
	$pfamIdToName{$pfamId} = $pfamName;
}
close FF;

# pfam位置读入
my ($pfamIdPosList, @pfamIdPos, $pfamIdPos, $pos, $pfamIdNamePosList, $pfamIdNamePos);
open FF, "<$origPfamInCdna";
# AT3G11450.1     PF00226[547,786]|PF00249[2020,2160]
while($line=<FF>){
	chomp($line);
	($trsptId, $pfamIdPosList) = split(/\t/, $line);
	@pfamIdPos = ();
	@pfamIdPos = split(/\|/, $pfamIdPosList);
	$pfamIdNamePosList = "";
	foreach $pfamIdPos(@pfamIdPos){
		if($pfamIdPos=~/(PF\d+)\[(.*)\]/){
			$pfamId = $1;
			$pos = $2;
			$pfamName = $pfamIdToName{$pfamId};
			$pfamIdNamePos = $pfamId . "[" . $pfamName . "]" . "[" . $pos . "]";
			if($pfamIdNamePosList eq ""){
				$pfamIdNamePosList = $pfamIdNamePos;
			}else{
				$pfamIdNamePosList .= "|" . $pfamIdNamePos;
			}
		}
	}
	$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"pfamIdNamePosList"} = $pfamIdNamePosList;
}
close FF;

# go读入
my ($goList);
open FF, "<$origGoUniqListInCdna";
# AT3G08840.13    GO:0005524|GO:0008716|GO:0046872 
while($line=<FF>){
	chomp($line);
	($trsptId, $goList) = split(/\t/, $line);
	$asOnTrsptHref->{$trsptId . "XXX" . "NA"}->{"goList"} = $goList;
}
close FF;








###################
## 将after as CdnaFasta读入hash
open FF, "<$asAlteredCdnaFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@trsptId = split(/ /, $1);
		$trsptId = $trsptId[0];
	}else{
		$asOnTrsptHref->{$trsptId}->{"seq"} .= $line;
	}
}
close FF;



# 将orf位置读入
open FF, "<$asAlteredOrfInCdna";
#AT3G53400.1     420     1817    420     422     1818    1820
#AT5G45110.2     515     2062    515     517     2063    2065
while($line=<FF>){
	chomp($line);
	($trsptId, $orfStart, $orfStop, $startCodonStart, $startCodonStop, $stopCodonStart, $stopCodonStop) = split(/\t/, $line);
	$asOnTrsptHref->{$trsptId}->{"orfStart"} = $orfStart;
	$asOnTrsptHref->{$trsptId}->{"orfStop"} = $orfStop;
	$asOnTrsptHref->{$trsptId}->{"startCodonStart"} = $startCodonStart;
	$asOnTrsptHref->{$trsptId}->{"startCodonStop"} = $startCodonStop;
	$asOnTrsptHref->{$trsptId}->{"stopCodonStart"} = $stopCodonStart;
	$asOnTrsptHref->{$trsptId}->{"stopCodonStop"} = $stopCodonStop;
}
close FF;

# 将pfam的名称读入
open FF, "<$asAlteredPfamNameFile";
# trsptId pfamId  pfamName        beginInCdna     endInCdna
# AT5G45110.2XXXATHASE00001     PF00651 BTB/POZ domain  680     1054
<FF>; 
while($line=<FF>){
	chomp($line);
	($trsptId, $pfamId, $pfamName, $beginInCdna, $endInCdna) = split(/\t/, $line);
	$pfamIdToName{$pfamId} = $pfamName;
}
close FF;


# pfam位置读入
open FF, "<$asAlteredPfamInCdna";
# AT3G11450.1XXXATHASE00001     PF00226[547,786]|PF00249[2020,2160]
while($line=<FF>){
	chomp($line);
	($trsptId, $pfamIdPosList) = split(/\t/, $line);
	@pfamIdPos = ();
	@pfamIdPos = split(/\|/, $pfamIdPosList);
	$pfamIdNamePosList = "";
	foreach $pfamIdPos(@pfamIdPos){
		if($pfamIdPos=~/(PF\d+)\[(.*)\]/){
			$pfamId = $1;
			$pos = $2;
			$pfamName = $pfamIdToName{$pfamId};
			$pfamIdNamePos = $pfamId . "[" . $pfamName . "]" . "[" . $pos . "]";
			if($pfamIdNamePosList eq ""){
				$pfamIdNamePosList = $pfamIdNamePos;
			}else{
				$pfamIdNamePosList .= "|" . $pfamIdNamePos;
			}
		}
	}
	$asOnTrsptHref->{$trsptId}->{"pfamIdNamePosList"} = $pfamIdNamePosList;
}
close FF;

# go读入
open FF, "<$asAlteredGoUniqListInCdna";
# AT3G08840.13    GO:0005524|GO:0008716|GO:0046872 
while($line=<FF>){
	chomp($line);
	($trsptId, $goList) = split(/\t/, $line);
	$asOnTrsptHref->{$trsptId}->{"goList"} = $goList;
}
close FF;


# 生成符合mysql插入的TSV文件
my @trsptIdAndAsId = keys(%asOnTrspt);
my ($nameFieldString, $valueFieldString, $trsptIdAndAsId, $asId);
open WW, ">$outputFunctionAffectionTsv";
foreach $trsptIdAndAsId(@trsptIdAndAsId){
	($trsptId, $asId) = split(/XXX/, $trsptIdAndAsId);
	$nameFieldString = join(", ", "trsptId", "asId", "orfStart", "orfStop", "startCodonStart", "startCodonStop", "stopCodonStart", "stopCodonStop", "pfamIdNamePosList", "goList", "sequence");

	if(not(exists($asOnTrsptHref->{$trsptIdAndAsId}->{"pfamIdNamePosList"}))){
		$asOnTrsptHref->{$trsptIdAndAsId}->{"pfamIdNamePosList"} = "NA";
	}
	$asOnTrsptHref->{$trsptIdAndAsId}->{"pfamIdNamePosList"}=~s/ /_/g;
	if(not(exists($asOnTrsptHref->{$trsptIdAndAsId}->{"goList"}))){
		$asOnTrsptHref->{$trsptIdAndAsId}->{"goList"} = "NA";
	}
	$valueFieldString = join(", ", $trsptId, $asId, $asOnTrsptHref->{$trsptIdAndAsId}->{"orfStart"}, $asOnTrsptHref->{$trsptIdAndAsId}->{"orfStop"}, $asOnTrsptHref->{$trsptIdAndAsId}->{"startCodonStart"}, $asOnTrsptHref->{$trsptIdAndAsId}->{"startCodonStop"}, $asOnTrsptHref->{$trsptIdAndAsId}->{"stopCodonStart"}, $asOnTrsptHref->{$trsptIdAndAsId}->{"stopCodonStop"}, $asOnTrsptHref->{$trsptIdAndAsId}->{"pfamIdNamePosList"}, $asOnTrsptHref->{$trsptIdAndAsId}->{"goList"}, $asOnTrsptHref->{$trsptIdAndAsId}->{"seq"});
	print WW $nameFieldString . "___" . $valueFieldString . "\n";
}
close WW;
