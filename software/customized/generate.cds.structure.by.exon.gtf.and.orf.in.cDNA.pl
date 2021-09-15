#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtf  \\\n" .
                "--orfInCdna \\\n" .
		"--trsptCdsStrucGtf \n";
	exit;
}

my ($gtf, $orfInCdna, $trsptCdsStrucGtf);

GetOptions(
        'gtf=s'=>\$gtf,
        'orfInCdna=s'=>\$orfInCdna,
	'trsptCdsStrucGtf=s'=>\$trsptCdsStrucGtf,
);

# 将gtf读入到hash中
my (%trsptExonStrct, $trsptExonStrctHref, $cmd, $lineText, @line, $line, @field);
my ($geneId, $trsptId, $i);
$trsptExonStrctHref = \%trsptExonStrct;
#$cmd = "grep -P \"\\texon\\t\" " . $gtf;
$cmd = "cat " . $gtf;
$lineText = `$cmd`;
@line = split(/\n/, $lineText);
foreach $line(@line){
	@field = ();
	@field = split(/\t/, $line);
	if($field[8]=~/gene_id "(.*?)"; transcript_id "(.*?)";/){
		$geneId = $1;
		$trsptId = $2;
	}

	$trsptExonStrctHref->{$trsptId}->{"strand"} = $field[6];
	$trsptExonStrctHref->{$trsptId}->{"geneId"} = $geneId;
	$trsptExonStrctHref->{$trsptId}->{"chr"} = $field[0];

	# 将转录本及其exon行都登记下来
	$line =~s/\tStringTie\t/\tplantAS\t/g;
	$trsptExonStrctHref->{$trsptId}->{"origTrsptAndExonText"} .= $line . "\n";
	# 提取exon的坐标
	if($field[2] eq "exon"){
		if(not exists($trsptExonStrctHref->{$trsptId}->{"exon"})){
			$trsptExonStrctHref->{$trsptId}->{"exon"} = $field[3] . ".." . $field[4];
		}else{
			$trsptExonStrctHref->{$trsptId}->{"exon"} .= "," . $field[3] . ".." . $field[4];
		}
	}
}



# 输出每个trspt对应的gtf
open WW, ">$trsptCdsStrucGtf";


################################################
# 将读取trspt的orf
#
my ($orfLine, $orfBegin, $orfEnd, $cdsExonSeries, $fiveUtrExonSeries, $threeUtrExonSeries, $trsptCdsBegin, $trsptCdsEnd);
open FF, "<$orfInCdna";
# trsptId orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonEnd
# AT3G53400.1     420     1817    420     422     1818    1820
# trsptId orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonEnd    cDNAseqLen
# SRX399565.4024.3        1011    1595    1011    1013    1596    1598    1733
# SRX1660584.17147.2      80      877     80      82      878     880     1744
#<FF>;
<FF>;
while($orfLine=<FF>){
	@field = split(/\t/, $orfLine);
	$trsptId = $field[0];
	$orfBegin = $field[1];
	$orfEnd = $field[2];
	#print "trsptId:" . $trsptId . ", strand:" . $trsptExonStrctHref->{$trsptId}->{"strand"} . ", orfB:[" . $orfBegin . "-" . $orfEnd . "]\n";
	# 获得转录本的CDS exon 串结构
	#print "trspt exon: " . $trsptExonStrctHref->{$trsptId}->{"exon"} . "\n";
	($cdsExonSeries, $fiveUtrExonSeries, $threeUtrExonSeries) = ("", "", "");
	&get5utrCds3utrExonSeries($trsptExonStrctHref->{$trsptId}->{"exon"}, $trsptExonStrctHref->{$trsptId}->{"strand"}, $orfBegin, $orfEnd, \$cdsExonSeries, \$fiveUtrExonSeries, \$threeUtrExonSeries);
	#<STDIN>;
	# 生成gtf格式的文本输出
	print WW $trsptExonStrctHref->{$trsptId}->{"origTrsptAndExonText"};
	if($fiveUtrExonSeries ne ""){
		print WW &generateGtf5utrExonText($trsptId, $trsptExonStrctHref->{$trsptId}->{"geneId"}, $fiveUtrExonSeries, $trsptExonStrctHref->{$trsptId}->{"chr"}, $trsptExonStrctHref->{$trsptId}->{"strand"}) . "\n";
	}
	if($cdsExonSeries ne ""){
		print WW &generateGtfCdsExonText($trsptId, $trsptExonStrctHref->{$trsptId}->{"geneId"}, $cdsExonSeries, $trsptExonStrctHref->{$trsptId}->{"chr"}, $trsptExonStrctHref->{$trsptId}->{"strand"}) . "\n";
	}
	if($threeUtrExonSeries ne ""){
		print WW &generateGtf3utrExonText($trsptId, $trsptExonStrctHref->{$trsptId}->{"geneId"}, $threeUtrExonSeries, $trsptExonStrctHref->{$trsptId}->{"chr"}, $trsptExonStrctHref->{$trsptId}->{"strand"}) . "\n";
	}
	
	
}
close FF;

# 获得trspt的CDS格式的外显子
sub generateGtfCdsExonText{
	my ($trsptId, $geneId, $exonSeries, $chr, $strand) = @_;
	my (@exon, $exon, @coord, $coord, $i, $j, $gtfCdsExonText, $exonNumber);
	@exon = split(/,/, $exonSeries);
	$exonNumber = 1;
	$gtfCdsExonText = "";
	foreach $exon(@exon){
		@coord = split(/\.\./, $exon);
		if($gtfCdsExonText eq ""){
			$gtfCdsExonText = join("\t", $chr, "plantAS", "CDS", $coord[0], $coord[1], ".", $strand, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\";");
		}else{
			$gtfCdsExonText .= "\n" . join("\t", $chr, "plantAS", "CDS", $coord[0], $coord[1], ".", $strand, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\";");
		}
		$exonNumber++;
	}
	return $gtfCdsExonText;
}

# 获得trspt的5utr格式的外显子
sub generateGtf5utrExonText{
	my ($trsptId, $geneId, $exonSeries, $chr, $strand) = @_;
	my (@exon, $exon, @coord, $coord, $i, $j, $gtf5utrExonText, $exonNumber);
	@exon = split(/,/, $exonSeries);
	$exonNumber = 1;
	$gtf5utrExonText = "";
	foreach $exon(@exon){
		@coord = split(/\.\./, $exon);
		if($gtf5utrExonText eq ""){
			$gtf5utrExonText = join("\t", $chr, "plantAS", "5UTR", $coord[0], $coord[1], ".", $strand, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\";");
		}else{
			$gtf5utrExonText .= "\n" . join("\t", $chr, "plantAS", "5UTR", $coord[0], $coord[1], ".", $strand, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\";");
		}
		$exonNumber++;
	}
	return $gtf5utrExonText;
}

# 获得trspt的3utr格式的外显子
sub generateGtf3utrExonText{
	my ($trsptId, $geneId, $exonSeries, $chr, $strand) = @_;
	my (@exon, $exon, @coord, $coord, $i, $j, $gtf3utrExonText, $exonNumber);
	@exon = split(/,/, $exonSeries);
	$exonNumber = 1;
	$gtf3utrExonText = "";
	foreach $exon(@exon){
		@coord = split(/\.\./, $exon);
		if($gtf3utrExonText eq ""){
			$gtf3utrExonText = join("\t", $chr, "plantAS", "3UTR", $coord[0], $coord[1], ".", $strand, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\";");
		}else{
			$gtf3utrExonText .= "\n" . join("\t", $chr, "plantAS", "3UTR", $coord[0], $coord[1], ".", $strand, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\";");
		}
		$exonNumber++;
	}
	return $gtf3utrExonText;
}




# 根据 orf在exon串中的位置，给出5UTR, 3UTR 和CDS exon 串位置信息
sub get5utrCds3utrExonSeries{
	my ($exonSeries, $strand, $orfBeginInCdna, $orfEndInCdna, $rltCdsExonSeries, $rltFiveUtrExonSeries, $rltThreeUtrExonSeries) = @_;
	my ($cdsExonSeries, $fiveUtrExonSeries, $threeUtrExonSeries);
	my (@exon, $exon, @coor, $coor, @exonCoor, $i, $j, $exonNum);
	my ($orfBeginExonId, $orfEndExonId,  $orfBeginInExon, $orfEndInExon);

	#print "original exonSeries: " . $exonSeries . "\n";
	#print "strand: " . $strand . "\n";
	#print "orfBeginInCdna: " . $orfBeginInCdna . "\n";
	#print "orfEndInCdna: " . $orfEndInCdna . "\n";;

	($$rltCdsExonSeries, $$rltFiveUtrExonSeries, $$rltThreeUtrExonSeries) = ("", "", "");

	@exon = split(/,/, $exonSeries);
	# 计算外显子长度放在exonCoor[2]
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

	# 计算外显子积累长度，放在exonCoor[3]
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
	
	# exonCoor: [0] exonStart [1] exonStop [2] exonLen [3] exonAccumulatedLen


	# 获得orfBegin和orfEnd分别在外显子串中的exon编号
	# $orfBeginExonId
	for($orfBeginExonId=0; $orfBeginExonId<=$#exonCoor; $orfBeginExonId++){
		if($orfBeginInCdna<=$exonCoor[$orfBeginExonId][3]){
			last;
		}
	}
	# orfEndExonId
	for($orfEndExonId=0; $orfEndExonId<=$#exonCoor; $orfEndExonId++){
		if($orfEndInCdna <= $exonCoor[$orfEndExonId][3]){
			last;
		}
	}


	# 生成cds, 5utr, 3utr的外显子串
	if($strand eq "+"){
		
		# 计算orfBegin在orfBeginExonId外显子中的坐标
		$orfBeginInExon = $exonCoor[$orfBeginExonId][1] - ($exonCoor[$orfBeginExonId][3] - $orfBeginInCdna);
		# 计算orfEndInExon在orfEndExonId外显子中的坐标
		$orfEndInExon = $exonCoor[$orfEndExonId][1] - ($exonCoor[$orfEndExonId][3] - $orfEndInCdna);

		# 收集5'utr exon series
		$fiveUtrExonSeries = "";
		if($orfBeginExonId > 0){
			for($i=0; $i<=$orfBeginExonId-1; $i++){
				$fiveUtrExonSeries .= $exonCoor[$i][0] . ".." . $exonCoor[$i][1] . ",";
			}
		}
		if($orfBeginInExon > $exonCoor[$orfBeginExonId][0]){
		# 5UTR结束于当前orfBeginExonId外显子内部，因此需要补上当前orfBeginExonId上的剩余的左侧部分
			my $last5UtrExonEnd = $orfBeginInExon - 1;
			$fiveUtrExonSeries .= $exonCoor[$orfBeginExonId][0] . ".." . $last5UtrExonEnd;
		}elsif($orfBeginInExon == $exonCoor[$orfBeginExonId][0]){
		# 5UTR结束于orfBeginExonId之前的外显子，因此不需要补充片段外显子
			if($fiveUtrExonSeries ne ""){
				$fiveUtrExonSeries = substr($fiveUtrExonSeries, 0 , length($fiveUtrExonSeries) - 1);
			}
		}

		# 收集CDS的exon series
		if($orfBeginExonId != $orfEndExonId){
		# CDS只有1个外显子
			$cdsExonSeries = $orfBeginInExon . ".." . $exonCoor[$orfBeginExonId][1];
			for($i=$orfBeginExonId+1; $i<=$orfEndExonId-1; $i++){
				$cdsExonSeries .= "," . $exonCoor[$i][0] . ".." . $exonCoor[$i][1];
			}
			$cdsExonSeries .= "," . $exonCoor[$orfEndExonId][0] . ".." . $orfEndInExon;
		}else{
		# CDS只有1个外显子
				$cdsExonSeries = $orfBeginInExon . ".." . $orfEndInExon;
		}

		# 收集3'utr exon series
		$threeUtrExonSeries = "";
		if($orfEndExonId < $#exonCoor){
			for($i=$orfEndExonId+1; $i<=$#exonCoor; $i++){
				$threeUtrExonSeries .= "," . $exonCoor[$i][0] . ".." . $exonCoor[$i][1];
			}
		}
		if($orfEndInExon < $exonCoor[$orfEndExonId][1]){
		# 3UTR结束于当前orfEndExonId外显子内部，因此需要补上当前orfEndExonId外显子剩余的右侧部分
			my $first3UtrExonEnd = $orfEndInExon + 1;
			$threeUtrExonSeries = $first3UtrExonEnd . ".." . $exonCoor[$orfEndExonId][1] . $threeUtrExonSeries;
		}elsif($orfEndInExon == $exonCoor[$orfEndExonId][1]){
		# 第1个3UTR外显子从$orfEndExonId之后的外显子开始，因此不需要补充片段外显子
			if($threeUtrExonSeries ne ""){
			# 如果存在多个完整的3UTR外显子，那么需要掉前导","
				$threeUtrExonSeries = substr($threeUtrExonSeries, 1);
			}
		}

	}else{
		$orfBeginInExon = $exonCoor[$orfBeginExonId][0] + ($exonCoor[$orfBeginExonId][3] - $orfBeginInCdna);
		$orfEndInExon = $exonCoor[$orfEndExonId][0] + ($exonCoor[$orfEndExonId][3] - $orfEndInCdna);

		# 收集5UTR外显子串
		$fiveUtrExonSeries = "";
		if($orfBeginExonId > 0){
			for($i=0; $i<=$orfBeginExonId-1; $i++){
				$fiveUtrExonSeries .= $exonCoor[$i][0] . ".." . $exonCoor[$i][1] . ",";
			}
		}
		if($orfBeginInExon < $exonCoor[$orfBeginExonId][1]){
		# 5UTR结束于当前orfEndExonId外显子内部，因此需要补上当前外显子剩余的右侧部分
			my $last5UtrExonBegin = $orfBeginInExon + 1;
			$fiveUtrExonSeries .= $last5UtrExonBegin . ".." . $exonCoor[$orfBeginExonId][1];
		}elsif($orfBeginInExon == $exonCoor[$orfBeginExonId][1]){
		# 5UTR结束于当前orfEndExonId外显子的前面，因此不需要补充片段外显子。但是需要处理掉后缀","
			$fiveUtrExonSeries = substr($fiveUtrExonSeries, 0, length($fiveUtrExonSeries) - 1);
		}
		
		# 收集CDS外显子
		if($orfBeginExonId != $orfEndExonId){
			# CDS有多个外显子
			$cdsExonSeries = $exonCoor[$orfBeginExonId][0] . ".." . $orfBeginInExon;
			for($i=$orfBeginExonId+1; $i<=$orfEndExonId-1; $i++){
				$cdsExonSeries .= "," . $exonCoor[$i][0] . ".." . $exonCoor[$i][1];
			}
			$cdsExonSeries .= "," . $orfEndInExon . ".." . $exonCoor[$orfEndExonId][1];
		}else{
			# CDS只有1个外显子
			$cdsExonSeries = $orfEndInExon . ".." . $orfBeginInExon;
		}

		# 收集3UTR外显子串
		$threeUtrExonSeries = "";
		if($orfEndExonId < $#exonCoor){
			for($i=$orfEndExonId+1; $i<=$#exonCoor; $i++){
				$threeUtrExonSeries .= "," . $exonCoor[$i][0] . ".." . $exonCoor[$i][1];
			}
		}
		if($orfEndInExon > $exonCoor[$orfEndExonId][0]){
		# 3UTR结束于当前orfEndExonId外显子内部，则需要补上当前orfEndInExonId外显子剩余的左侧部分
			my $first3UtrExonEnd = $orfEndInExon - 1;
			$threeUtrExonSeries = $exonCoor[$orfEndExonId][0] . ".." . $first3UtrExonEnd . $threeUtrExonSeries;
		}elsif($orfEndInExon == $exonCoor[$orfEndExonId][0]){
		# 3UTR结束于当前orfEndExonId外显子之前，那么不需要补充片段外显子，但是需要去掉前导","
			$threeUtrExonSeries = substr($threeUtrExonSeries, 1);
		}
		
	}
#	print "rltFiveUtrExonSeries: " . $$rltFiveUtrExonSeries . "\n";
#	print "rltCdsExonSeries: " . $$rltCdsExonSeries . "\n";
#	print "rltThreeUtrExonSeries: " . $$rltThreeUtrExonSeries . "\n";
#	<STDIN>;
	($$rltFiveUtrExonSeries, $$rltCdsExonSeries, $$rltThreeUtrExonSeries) = ($fiveUtrExonSeries, $cdsExonSeries, $threeUtrExonSeries);
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
