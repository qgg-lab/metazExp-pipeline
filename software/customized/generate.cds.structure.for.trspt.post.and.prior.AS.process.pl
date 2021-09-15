#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtf  \\\n" .
                "--beforeOrfInCdna \\\n" .
                "--afterOrfInCdna \\\n" .
		"--relationOfAsAndTrspt \\\n" .
		"--trsptCdsStrucGtf \n";
	exit;
}

my ($gtf, $beforeOrfInCdna, $afterOrfInCdna, $relationOfAsAndTrspt, $trsptCdsStrucGtf);

GetOptions(
        'gtf=s'=>\$gtf,
        'beforeOrfInCdna=s'=>\$beforeOrfInCdna,
        'afterOrfInCdna=s'=>\$afterOrfInCdna,
        'relationOfAsAndTrspt=s'=>\$relationOfAsAndTrspt,
	'trsptCdsStrucGtf=s'=>\$trsptCdsStrucGtf,
);

# 将gtf读入到hash中
my (%trsptExonStrct, $trsptExonStrctHref, $cmd, $exonLineText, @exonLine, $exonLine, @field);
my ($geneId, $trsptId, $i);
$trsptExonStrctHref = \%trsptExonStrct;
$cmd = "grep -P \"\\texon\\t\" " . $gtf;
$exonLineText = `$cmd`;
@exonLine = split(/\n/, $exonLineText);
foreach $exonLine(@exonLine){
	@field = ();
	@field = split(/\t/, $exonLine);
	if($field[8]=~/gene_id "(.*?)"; transcript_id "(.*?)";/){
		$geneId = $1;
		$trsptId = $2;
	}

	$trsptExonStrctHref->{$trsptId}->{"strand"} = $field[6];
	$trsptExonStrctHref->{$trsptId}->{"geneId"} = $geneId;
	$trsptExonStrctHref->{$trsptId}->{"chr"} = $field[0];
	if(not exists($trsptExonStrctHref->{$trsptId}->{"exon"})){
		$trsptExonStrctHref->{$trsptId}->{"exon"} = $field[3] . ".." . $field[4];
	}else{
		$trsptExonStrctHref->{$trsptId}->{"exon"} .= "," . $field[3] . ".." . $field[4];
	}
}


# 将trspt和AS之间的映射关系读入到hash中
my (%asAltedTrspt, $asAltedTrsptHref, @titleField, @valueField, $relationLine, %tmpHash, $tmpHashHref, @longAltTrsptId, @shortAltTrsptId, $longAltTrsptIdList, $shortAltTrsptIdList, $longAltTrsptId, $shortAltTrsptId, $trsptIdAndAsId, $origExonSeries, $newExonSeries, $altBeginInExonSeries);
$tmpHashHref=\%tmpHash;
$asAltedTrsptHref=\%asAltedTrspt;
#print $relationOfAsAndTrspt;
#<STDIN>;
open FF, "<$relationOfAsAndTrspt";
# ASID    chr     strand  longAltExonSeries       longAltEnsemblTrsptIdList       longAltImprovedTrsptIdList      shortAltExonSeries      shortAltEnsemblTrsptIdList      shortAltImprovedTrsptIdList
# ATHAA3SS0000002126      1       +       317294..317439,317553..317642   AT1G01920.2,AT1G01920.3,AT1G01920.5,AT1G01920.4 SRX1182480.125.5        317294..317439,317561..317642   -       -
$relationLine = <FF>;
#print "Title:\n";
#print $relationLine;
#<STDIN>;
chomp($relationLine);
@titleField = split(/\t/, $relationLine);
while($relationLine=<FF>){
#	print $relationLine;
#	<STDIN>;
	chomp($relationLine);
	@valueField = split(/\t/, $relationLine);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHashHref->{$titleField[$i]} = $valueField[$i];
	}

	# 处理longAltTrsptIdList
	$longAltTrsptIdList = "";
	if($tmpHashHref->{"longAltEnsemblTrsptIdList"} ne "-" and $tmpHashHref->{"longAltImprovedTrsptIdList"} ne "-"){
		$longAltTrsptIdList = $tmpHashHref->{"longAltEnsemblTrsptIdList"} . "," . $tmpHashHref->{"longAltImprovedTrsptIdList"}
	}elsif($tmpHashHref->{"longAltEnsemblTrsptIdList"} ne "-" and $tmpHashHref->{"longAltImprovedTrsptIdList"} eq "-"){
		$longAltTrsptIdList = $tmpHashHref->{"longAltEnsemblTrsptIdList"};
	}elsif($tmpHashHref->{"longAltEnsemblTrsptIdList"} eq "-" and $tmpHashHref->{"longAltImprovedTrsptIdList"} ne "-"){
		$longAltTrsptIdList=$tmpHashHref->{"longAltImprovedTrsptIdList"};
	}
#	print "longAltTrsptIdList:" . $longAltTrsptIdList;
#	<STDIN>;
	if($longAltTrsptIdList ne ""){
		@longAltTrsptId = ();
		@longAltTrsptId = split(/,/, $longAltTrsptIdList);
		foreach $longAltTrsptId(@longAltTrsptId){
			$trsptIdAndAsId = $longAltTrsptId . "XXX" . $tmpHashHref->{"ASID"};
			$origExonSeries = $trsptExonStrctHref->{$longAltTrsptId}->{"exon"};
			$altBeginInExonSeries = index($origExonSeries, $tmpHashHref->{"longAltExonSeries"});
			$newExonSeries = substr($origExonSeries, 0, $altBeginInExonSeries) . $tmpHashHref->{"shortAltExonSeries"} . substr($origExonSeries, $altBeginInExonSeries + length($tmpHashHref->{"longAltExonSeries"}));
			$asAltedTrsptHref->{$trsptIdAndAsId}->{"exon"} = $newExonSeries;
			$asAltedTrsptHref->{$trsptIdAndAsId}->{"chr"} = $tmpHashHref->{"chr"};
			$asAltedTrsptHref->{$trsptIdAndAsId}->{"strand"} = $tmpHashHref->{"strand"};
			$asAltedTrsptHref->{$trsptIdAndAsId}->{"geneId"} = $trsptExonStrctHref->{$longAltTrsptId}->{"geneId"};
#			if($longAltTrsptId eq "Zm00001d006457_T003" and $tmpHashHref->{"ASID"} eq "ZMAYMXE0000000583"){
#				print "contain longAlt:" . $tmpHashHref->{"longAltExonSeries"} . "\n";
#				print $longAltTrsptId . " origExonSeries:" . $origExonSeries . "\n";;
#				print $longAltTrsptId . " newExonSeries :" . $newExonSeries;
#				<STDIN>;
#			}
		}
	}

	# 处理shortAltTrsptIdList
	$shortAltTrsptIdList = "";
	if($tmpHashHref->{"shortAltEnsemblTrsptIdList"} ne "-" and $tmpHashHref->{"shortAltImprovedTrsptIdList"} ne "-"){
		$shortAltTrsptIdList = $tmpHashHref->{"shortAltEnsemblTrsptIdList"} . "," . $tmpHashHref->{"shortAltImprovedTrsptIdList"}
	}elsif($tmpHashHref->{"shortAltEnsemblTrsptIdList"} ne "-" and $tmpHashHref->{"shortAltImprovedTrsptIdList"} eq "-"){
		$shortAltTrsptIdList = $tmpHashHref->{"shortAltEnsemblTrsptIdList"};
	}elsif($tmpHashHref->{"shortAltEnsemblTrsptIdList"} eq "-" and $tmpHashHref->{"shortAltImprovedTrsptIdList"} ne "-"){
		$shortAltTrsptIdList = $tmpHashHref->{"shortAltImprovedTrsptIdList"};
	}
#	print "shortAltTrsptIdList:" . $shortAltTrsptIdList;
#	<STDIN>;
	if($shortAltTrsptIdList ne ""){
		@shortAltTrsptId = ();
		@shortAltTrsptId = split(/,/, $shortAltTrsptIdList);
		foreach $shortAltTrsptId(@shortAltTrsptId){
			$trsptIdAndAsId = $shortAltTrsptId . "XXX" . $tmpHashHref->{"ASID"};
			$origExonSeries = $trsptExonStrctHref->{$shortAltTrsptId}->{"exon"};
			$altBeginInExonSeries = index($origExonSeries, $tmpHashHref->{"shortAltExonSeries"});
			$newExonSeries = substr($origExonSeries, 0, $altBeginInExonSeries) . $tmpHashHref->{"longAltExonSeries"} . substr($origExonSeries, $altBeginInExonSeries + length($tmpHashHref->{"shortAltExonSeries"}));
			$asAltedTrsptHref->{$trsptIdAndAsId}->{"exon"} = $newExonSeries;
			$asAltedTrsptHref->{$trsptIdAndAsId}->{"chr"} = $tmpHashHref->{"chr"};
			$asAltedTrsptHref->{$trsptIdAndAsId}->{"strand"} = $tmpHashHref->{"strand"};
			$asAltedTrsptHref->{$trsptIdAndAsId}->{"geneId"} = $trsptExonStrctHref->{$shortAltTrsptId}->{"geneId"};
#			if($shortAltTrsptId eq "Zm00001d006457_T003" and $tmpHashHref->{"ASID"} eq "ZMAYMXE0000000583"){
#				print "contain shortAlt:" . $tmpHashHref->{"shortAltExonSeries"} . "\n";
#				print $shortAltTrsptId . " origExonSeries:" . $origExonSeries . "\n";;
#				print $shortAltTrsptId . " newExonSeries :" . $newExonSeries;
#				<STDIN>;
#			}
		}
	}


}
close FF;

# 输出每个trspt对应的gtf
open WW, ">$trsptCdsStrucGtf";


################################################
# 将读取AS改变前的trspt的orf
#
my ($orfLine, $orfBegin, $orfEnd, $cdsExonSeries, $trsptCdsBegin, $trsptCdsEnd);
open FF, "<$beforeOrfInCdna";
# trsptId orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonEnd
# AT3G53400.1     420     1817    420     422     1818    1820
#<FF>;
while($orfLine=<FF>){
	@field = split(/\t/, $orfLine);
	$trsptId = $field[0];
	$orfBegin = $field[1];
	$orfEnd = $field[2];
	#print "trsptId:" . $trsptId . ", strand:" . $trsptExonStrctHref->{$trsptId}->{"strand"} . ", orfB:[" . $orfBegin . "-" . $orfEnd . "]\n";
	# 获得转录本的CDS exon 串结构
	#print "trspt exon: " . $trsptExonStrctHref->{$trsptId}->{"exon"} . "\n";
	$cdsExonSeries = &getCdsExonSeries($trsptExonStrctHref->{$trsptId}->{"exon"}, $trsptExonStrctHref->{$trsptId}->{"strand"}, $orfBegin, $orfEnd);
	#print "trspt cds exon：" . $cdsExonSeries;
	#<STDIN>;
	# 获得转录本的起始位置和结束位置
	&getMaxAndMinCoor($cdsExonSeries, \$trsptCdsBegin, \$trsptCdsEnd);
	print WW join("\t", $trsptExonStrctHref->{$trsptId}->{"chr"}, "plantAS", "transcript", $trsptCdsBegin, $trsptCdsEnd, ".", $trsptExonStrctHref->{$trsptId}->{"strand"}, ".", "gene_id \"" . $trsptExonStrctHref->{$trsptId}->{"geneId"} . "\"; transcript_id \"" . $trsptId . "\";") . "\n";
	# 生成gtf格式的文本输出
	print WW &generateGtfExonText($trsptId, $trsptExonStrctHref->{$trsptId}->{"geneId"}, $cdsExonSeries, $trsptExonStrctHref->{$trsptId}->{"chr"}, $trsptExonStrctHref->{$trsptId}->{"strand"}) . "\n";
	
}
close FF;

# 读取AS改变后的trspt的orf
open FF, "<$afterOrfInCdna";
<FF>;
# trsptId orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonEnd
# AT5G51710.2XXXATHASE0000033759  213     1907    213     215     1908    1910
# AT2G06005.1XXXATHASE0000043648  290     1354    290     292     1355    1357
while($orfLine=<FF>){
	@field = split(/\t/, $orfLine);
	$trsptId = $field[0];
	$orfBegin = $field[1];
	$orfEnd = $field[2];
#	if($trsptId eq "AT3G15390.1XXXATHASE0000035094"){
#		print "orfLine:" . $orfLine . "\n";
#		print "trsptId:" .$trsptId . "\n";
#		print "orfBegin:" . $orfBegin . "\n";
#		print "orfEnd:" . $orfEnd . "\n";
	# 获得转录本的CDS exon 串结构
#		print "exonSeries:" . $asAltedTrsptHref->{$trsptId}->{"exon"} . "\n";
#	}
	$cdsExonSeries = &getCdsExonSeries($asAltedTrsptHref->{$trsptId}->{"exon"}, $asAltedTrsptHref->{$trsptId}->{"strand"}, $orfBegin, $orfEnd);
#	if($trsptId eq "AT3G15390.1XXXATHASE0000035094"){
#		print "cdsExonSeries: " . $cdsExonSeries . "\n";
#		<STDIN>;
#	}
	# 获得转录本的起始位置和结束位置
	&getMaxAndMinCoor($cdsExonSeries, \$trsptCdsBegin, \$trsptCdsEnd);
	print WW join("\t", $asAltedTrsptHref->{$trsptId}->{"chr"}, "plantAS", "transcript", $trsptCdsBegin, $trsptCdsEnd, ".", $asAltedTrsptHref->{$trsptId}->{"strand"}, ".", "gene_id \"" . $asAltedTrsptHref->{$trsptId}->{"geneId"} . "\"; transcript_id \"" . $trsptId . "\";") . "\n";	
	# 生成gtf格式的文本输出
	my $originTrsptId = "";
	if($trsptId=~/(.*)XXX(.*)/){
		$originTrsptId = $1;
	}
	print WW &generateGtfExonText($trsptId, $asAltedTrsptHref->{$trsptId}->{"geneId"}, $cdsExonSeries, $asAltedTrsptHref->{$trsptId}->{"chr"}, $trsptExonStrctHref->{$originTrsptId}->{"strand"}) . "\n";

}
close FF;

# 获得trspt的gtf格式的外显子
sub generateGtfExonText{
	my ($trsptId, $geneId, $exonSeries, $chr, $strand) = @_;
	my (@exon, $exon, @coord, $coord, $i, $j, $gtfExonText, $exonNumber);
	@exon = split(/,/, $exonSeries);
	$exonNumber = 0;
	$gtfExonText = "";
	foreach $exon(@exon){
		@coord = split(/\.\./, $exon);
		if($gtfExonText eq ""){
			$gtfExonText = join("\t", $chr, "plantAS", "exon", $coord[0], $coord[1], ".", $strand, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\"; exon_number \"$exonNumber\";");
		}else{
			$gtfExonText .= "\n" . join("\t", $chr, "plantAS", "exon", $coord[0], $coord[1], ".", $strand, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\"; exon_number \"$exonNumber\";");
		}
		$exonNumber++;
	}
	return $gtfExonText;
}

# 根据 orf在exon串中的位置，给出CDS exon 串位置信息
sub getCdsExonSeries{
	my ($exonSeries, $strand, $orfBeginInCdna, $orfEndInCdna) = @_;
	my ($cdsExonSeries);
	my (@exon, $exon, @coor, $coor, @exonCoor, $i, $j, $exonNum);
	my ($orfBeginExonId, $orfEndExonId,  $orfBeginInExon, $orfEndInExon);
	@exon = split(/,/, $exonSeries);
	# 计算外显子长度放在exon[2]
	$exonNum = 0;
	foreach $exon(@exon){
		@coor = split(/\.\./, $exon);
		$exonCoor[$exonNum][0] = $coor[0];
		$exonCoor[$exonNum][1] = $coor[1];
		# 当前外显子的长度
		$exonCoor[$exonNum][2] = $coor[1] - $coor[0] + 1;
		#if($orfBeginInCdna == 298 and $orfEndInCdna == 1767 and $strand eq "-"){
			#print join("\t", "exon_" . $exonNum, $exonCoor[$exonNum][0], $exonCoor[$exonNum][1], $exonCoor[$exonNum][2]) . "\n";
		#}
		$exonNum++;
	}

	# 计算外显子积累长度，放在exon[3]
	$exonNum = 0;
	for($i=0; $i<=$#exonCoor; $i++){
		if($i==0){
			$exonCoor[$i][3] = $exonCoor[$i][2];
		}else{
			$exonCoor[$i][3] = $exonCoor[$i][2] + $exonCoor[$i-1][3];
		}
		if($orfBeginInCdna == 298 and $orfEndInCdna == 1767 and $strand eq "-"){
			#print "exon_$i:" . $exonCoor[$i][3] . "\n";
		}
	}

	# 获得orfBegin和orfEnd分别在外显子串中的exon编号
#	if($orfBeginInCdna == 298 and $orfEndInCdna == 1767 and $strand eq "-"){
#		print "orfBeginInCdna:" . $orfBeginInCdna . "\n";
#		print "orfEndInCdna:" . $orfEndInCdna . "\n";
#	}
	for($orfBeginExonId=0; $orfBeginExonId<=$#exonCoor; $orfBeginExonId++){
		if($orfBeginInCdna<=$exonCoor[$orfBeginExonId][3]){
			last;
		}
	}

	for($orfEndExonId=0; $orfEndExonId<=$#exonCoor; $orfEndExonId++){
		if($orfEndInCdna <= $exonCoor[$orfEndExonId][3]){
			last;
		}
	}
#	if($orfBeginInCdna == 298 and $orfEndInCdna == 1767 and $strand eq "-"){
#		print "orfBeginExonId:" . $orfBeginExonId . "\n";
#		print "orfEndExonId:" . $orfEndExonId . "\n";
#	}
	# 分正负链计算orf的起始和终止位置
	if($strand eq "+"){
		$orfBeginInExon = $exonCoor[$orfBeginExonId][1] - ($exonCoor[$orfBeginExonId][3] - $orfBeginInCdna);
		$orfEndInExon = $exonCoor[$orfEndExonId][1] - ($exonCoor[$orfEndExonId][3] - $orfEndInCdna);
		if($orfBeginExonId != $orfEndExonId){
			$cdsExonSeries = $orfBeginInExon . ".." . $exonCoor[$orfBeginExonId][1];
			for($i=$orfBeginExonId+1; $i<=$orfEndExonId-1; $i++){
				$cdsExonSeries .= "," . $exonCoor[$i][0] . ".." . $exonCoor[$i][1];
			}
			$cdsExonSeries .= "," . $exonCoor[$orfEndExonId][0] . ".." . $orfEndInExon;
		}else{
			$cdsExonSeries = $orfBeginInExon . ".." . $orfEndInExon;
		}
	}else{
		$orfBeginInExon = $exonCoor[$orfBeginExonId][0] + ($exonCoor[$orfBeginExonId][3] - $orfBeginInCdna);
		$orfEndInExon = $exonCoor[$orfEndExonId][0] + ($exonCoor[$orfEndExonId][3] - $orfEndInCdna);
		if($orfBeginExonId != $orfEndExonId){
			$cdsExonSeries = $exonCoor[$orfBeginExonId][0] . ".." . $orfBeginInExon;
			for($i=$orfBeginExonId+1; $i<=$orfEndExonId-1; $i++){
				$cdsExonSeries .= "," . $exonCoor[$i][0] . ".." . $exonCoor[$i][1];
			}
			$cdsExonSeries .= "," . $orfEndInExon . ".." . $exonCoor[$orfEndExonId][1];
		}else{
			$cdsExonSeries = $orfEndInExon . ".." . $orfBeginInExon;
		}
	}
#	if($orfBeginInCdna == 298 and $orfEndInCdna == 1767 and $strand eq "-"){
#		print "orfBeginInExon:" . $orfBeginInExon . "\n";
#		print "orfEndInExon:" . $orfEndInExon . "\n";
#		<STDIN>;
#	}
	# 生成CDS exon串. 无论正负链orBeginExonId始终小于orfEndExonId
	# 因为：正链外显子串从小到大排列，负链外显子串从大到小排列
	return $cdsExonSeries;
}

# 获得exon串中最大值和最小值
sub getMaxAndMinCoor{
	my ($exonSeries, $minCoor, $maxCoor) = @_;
	my (@exon, @coor, $exon, $coor);
	@exon = split(/,/, $exonSeries);
	$$maxCoor = 1;
	$$minCoor = 10000000000;
	foreach $exon(@exon){
		@coor = split(/\.\./, $exon);
		foreach $coor(@coor){
			if($$maxCoor < $coor){
				$$maxCoor = $coor;
			}
			if($$minCoor > $coor){
				$$minCoor = $coor;
			}
		}
	}
}
