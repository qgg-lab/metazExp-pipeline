#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--orfInCdnaTsv  \\\n" .
                "--pfamInPepTsv \\\n" .
                "--pfamInCdnaTsv \n";
	exit;
}

my ($orfInCdnaTsv, $pfamInPepTsv, $pfamInCdnaTsv);

GetOptions(
        'orfInCdnaTsv=s'=>\$orfInCdnaTsv,
        'pfamInPepTsv=s'=>\$pfamInPepTsv,
        'pfamInCdnaTsv=s'=>\$pfamInCdnaTsv,
);
my ($line, @field, %trsptPfamPosition, %trsptOrfPosition, $trsptPfamPositionHref, $trsptOrfPositionHref, $pepLen);
my ($trsptId, $pfamId, $pfamName, $pfamStartInPep, $pfamStopInPep, $pfamPositionId, @pfamPositionId);
$trsptPfamPositionHref = \%trsptPfamPosition;
$trsptOrfPositionHref = \%trsptOrfPosition;
# 读取pfam在pep中的位置 
open FF, "<$pfamInPepTsv";
# Bra015256.1     70950ca75d35e99a1eb56855a5e628e8        175     Pfam    PF04535 Domain of unknown function (DUF588)     7       154     4.1E-32 T       14-12-2019      IPR006702       Casparian strip membrane protein domain
while($line=<FF>){
#	print $line;
#	<STDIN>;
	chomp($line);
	@field = split(/\t/, $line);
	($trsptId, $pepLen, $pfamId, $pfamName, $pfamStartInPep, $pfamStopInPep) = ($field[0], $field[2], $field[4], $field[5], $field[6], $field[7]);
	# 每个pfamId用的是该pfam在pep中的起始和终止位置
	$pfamPositionId = $pfamStartInPep . "_" . $pfamStopInPep;
	$trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pepLen"} = $pepLen;
	$trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pfamId"} = $pfamId;
	$trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pfamName"} = $pfamName;
	$trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pfamStartInPep"} = $pfamStartInPep;
	$trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pfamStopInPep"} = $pfamStopInPep;
}
close FF;


# 读取orf在trspt中的位置
# 根据起始密码子在trspt中的位置以及pfam的位置，推算出pfam在trspt中的位置
my ($pfamStartInCdna, $pfamStopInCdna, $orfBegin);
open FF, "<$orfInCdnaTsv";
<FF>;
# trsptId		orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonEnd
# SRX2909953.16532.2    1378    	3195    1378    	1380    	3196    	3198
open WW, ">$pfamInCdnaTsv";
print WW join("\t", "trsptId", "pfamId", "pfamName", "beginInCdna", "endInCdna") . "\n";
# Zm00001d013026_T001	  PF09336 Vps4 C terminal oligomerisation domain 690 1200
while($line=<FF>){
#	print $line;
#	<STDIN>;
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$trsptId = $field[0];
	$orfBegin = $field[1];

	# 如果该trspt没有找到任何pfam，那么直接跳过
	next if(not exists($trsptPfamPositionHref->{$trsptId}));
#	print "trsptId:$trsptId\n";
	# 提取该trspt内所有的pfamPositonId
	@pfamPositionId = ();
	@pfamPositionId = keys(%{$trsptPfamPositionHref->{$trsptId}});
	foreach $pfamPositionId(@pfamPositionId){
		# (..- 1)为pfam之前的氨基酸个数，*3 为pfam之前核苷酸个数，再+orfBegin-1为pfam之前碱基的位置，再加1为pfam开始位置
		$pfamStartInCdna = $orfBegin + ($trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pfamStartInPep"} - 1)*3 - 1 + 1;
		$pfamStopInCdna = $orfBegin + $trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pfamStopInPep"} *3 - 1;
		print WW join("\t", $trsptId, $trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pfamId"}, $trsptPfamPositionHref->{$trsptId}->{$pfamPositionId}->{"pfamName"}, $pfamStartInCdna, $pfamStopInCdna) . "\n";
	}
}
close FF;
close WW;

