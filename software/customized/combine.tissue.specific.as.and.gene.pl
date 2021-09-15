#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--asTsvList \\\n" .
                "--specificAsRltTsv \\\n" .
                "--specificGeneRltTsv \\\n" .
		"--outputSpecificAsAndGeneTsv \n";
	exit;
}

my ($asTsvList, $specificAsRltTsv, $specificGeneRltTsv, $outputSpecificAsAndGeneTsv);

GetOptions(
	'asTsvList=s'=>\$asTsvList,
        'specificAsRltTsv=s'=>\$specificAsRltTsv,
        'specificGeneRltTsv=s'=>\$specificGeneRltTsv,
        'outputSpecificAsAndGeneTsv=s'=>\$outputSpecificAsAndGeneTsv
);


# 将as和gene之间映射关系读入hash
my (@asTsv, $asTsv, %asToGene, $asId, $geneId, @tmp, $line);
@asTsv=split(/,/, $asTsvList);
foreach $asTsv(@asTsv){
	open ASTSV, "<$asTsv";
	<ASTSV>;
	while($line=<ASTSV>){
		chomp($line);
		($asId, $geneId) = split(/\t/, $line);
		@tmp = split(/"/, $geneId);
		$geneId = $tmp[1];
		$asToGene{$asId}=$geneId;
	}
	close ASTSV;
}


# 将as特异性结果读入hash，只保留pvalue0.05为Y的asId
my (%as, $asHref, %gene, $geneHref, $tissue, $specific, $asId, $tissue, $psi, $sampleNum, $lowPvalue, $topPvalue, @tt, $minDltPsi, $dltPsi);
$asHref=\%as;
$geneHref = \%gene;
# pvalue0.01  pvalue0.03  pvalue0.05  ASID    		Regulate        Tissue=PSI|exptNum|minDltPsi      1stTissue|avgPSI|exptNum|pvalue
# N           N           N           SLYCRI0000000486  High    	fruit=0.522|4|0.123   		  seed|0.131|4|0.004      	  leaf|0.504|25|0.370
open FF, "<$specificAsRltTsv";
<FF>;
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);

	# 至少为0.05
	next if($tmp[2] eq "N");

	$asId = $tmp[3];
	$specific = $tmp[4];
	$tmp[5]=~/(.*)=(.*)\|(\d+)\|(.*)/;
	$tissue = $1;
	$psi = $2;
	$sampleNum = $3;
	$minDltPsi = $4;
	
	# 找出最高pvalue和最低pvalue
	$lowPvalue = 2;
	$topPvalue = -2;
	for(my $i=6; $i<=$#tmp; $i++){
		@tt = split(/\|/, $tmp[$i]);
		$lowPvalue = $tt[$#tt] if($lowPvalue > $tt[$#tt]);
		$topPvalue = $tt[$#tt] if($topPvalue < $tt[$#tt]);
	}

	# 将该as的特意信息登记下来
	$asHref->{$asId}->{"specific"} = $specific;
	$asHref->{$asId}->{"tissue"} = $tissue;
	$asHref->{$asId}->{"psi"} = $psi;
	$asHref->{$asId}->{"sampleNum"} = $sampleNum;
	$asHref->{$asId}->{"lowPvalue"} = $lowPvalue;
	$asHref->{$asId}->{"topPvalue"} = $topPvalue;
	$asHref->{$asId}->{"minDltPsi"} = $minDltPsi;
}
close FF;

# 将gene特异性结果读入，只保留pvalue0.05为Y的gene
# N       N       N       Solyc02g092260.3        High    root    15.313  30      leaf|13.489|26|0.131    seed|10.820|12|0.032
my ($tpm);
open FF, "<$specificGeneRltTsv";
<FF>;
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);

	# 至少为0.05
	next if($tmp[2] eq "N");

	$geneId = $tmp[3];
	$specific = $tmp[4];
	$tissue = $tmp[5];
	$tpm = $tmp[6];
	$sampleNum = $tmp[7];
	
	# 找出最高pvalue和最低pvalue
	$lowPvalue = 2;
	$topPvalue = -2;
	for(my $i=8; $i<=$#tmp; $i++){
		@tt = split(/\|/, $tmp[$i]);
		$lowPvalue = $tt[$#tt] if($lowPvalue > $tt[$#tt]);
		$topPvalue = $tt[$#tt] if($topPvalue < $tt[$#tt]);
	}

	# 将该as的特意信息登记下来
	$geneHref->{$geneId}->{"specific"} = $specific;
	$geneHref->{$geneId}->{"tissue"} = $tissue;
	$geneHref->{$geneId}->{"tpm"} = $tpm;
	$geneHref->{$geneId}->{"sampleNum"} = $sampleNum;
	$geneHref->{$geneId}->{"lowPvalue"} = $lowPvalue;
	$geneHref->{$geneId}->{"topPvalue"} = $topPvalue;
}
close FF;


my ($asSpecificPart, $geneSpecificPart);

# 将as的特异性输出
open OUTPUT, ">$outputSpecificAsAndGeneTsv";
# asId          tissue  specif  PSI     sampleNum minDltPsi low-pvalue  top pvalue  geneId  tissue  specific  TPM   low pvalue      top pvalue
# SLA3SS00001   fruit   high    0.12    0.001       0.003       sl00132 -       -         -     -               -
# -             -       -       -       -               -       sl00135 fruit   low       0.123 0.04            0.06
my @asId = keys(%as);
print OUTPUT join("\t", "asId", "tissue", "specific", "PSI" , "experiment", "minDltPsi", "low-Pvalue", "top-Pvalue", "geneId", "tissue", "specific", "TPM", "experiment", "low-Pvalue", "top-Pvalue") . "\n";
foreach $asId(@asId){

	# 如果最小平均PSI小于0.1那么就放弃
	next if($asHref->{$asId}->{"minDltPsi"}<0.1);

	$asSpecificPart = join("\t", $asId, $asHref->{$asId}->{"tissue"}, $asHref->{$asId}->{"specific"}, $asHref->{$asId}->{"psi"}, $asHref->{$asId}->{"sampleNum"}, $asHref->{$asId}->{"minDltPsi"}, $asHref->{$asId}->{"lowPvalue"}, $asHref->{$asId}->{"topPvalue"});

	$geneId = $asToGene{$asId};
	if(exists($gene{$geneId})){
		$geneSpecificPart = join("\t", $geneId, $geneHref->{$geneId}->{"tissue"}, $geneHref->{$geneId}->{"specific"}, $geneHref->{$geneId}->{"tpm"}, $geneHref->{$geneId}->{"sampleNum"}, $geneHref->{$geneId}->{"lowPvalue"}, $geneHref->{$geneId}->{"topPvalue"});
		# 将当前的gene删除掉
		delete($geneHref->{$geneId});
	}else{
		$geneSpecificPart = join("\t", "-", "-", "-", "-", "-", "-", "-");
	}
	
	# 输出
	print OUTPUT join("\t", $asSpecificPart, $geneSpecificPart) . "\n";
}


# 将gene的特异性输出
my @geneId = keys(%gene);
foreach $geneId(@geneId){
	$asSpecificPart = join("\t", "-", "-", "-", "-", "-", "-", "-", "-");
	$geneSpecificPart = join("\t", $geneId, $geneHref->{$geneId}->{"tissue"}, $geneHref->{$geneId}->{"specific"}, $geneHref->{$geneId}->{"tpm"}, $geneHref->{$geneId}->{"sampleNum"}, $geneHref->{$geneId}->{"lowPvalue"}, $geneHref->{$geneId}->{"topPvalue"});
	print OUTPUT join("\t", $asSpecificPart, $geneSpecificPart) . "\n";
}

close OUTPUT;
