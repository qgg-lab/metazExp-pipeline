#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputExonTsv \\\n" .
                "--inputCdsExonGtf \\\n" .
                "--outputMergedExonNotInUniqPep \\\n" .
                "--outputExonPosInTrspt \n";
	exit;
}

my ($inputExonTsv, $inputCdsExonGtf, $outputExonPosInCdsSeq, $outputMergedExonNotInUniqPep);

GetOptions(
        'inputExonTsv=s'=>\$inputExonTsv,
        'inputCdsExonGtf=s'=>\$inputCdsExonGtf,
        'outputExonPosInCdsSeq=s'=>\$outputExonPosInCdsSeq,
	'outputMergedExonNotInUniqPep=s'=>\$outputMergedExonNotInUniqPep
);

# 将gtf中所有的cds exon读入，构建一个以trsptId为关键字，所有exon串为值的hash
# 同时建立从exon到trsptId映射关系的hash
my (%trsptIdToExonSeries, @field, $trsptId, @attr, $attr, $cmd, @exonLine, $exonLine, $exonLineText, @tt, $exonString, @exon, $exon);
my (%exonInTrspt, $exonInTrsptHref);
$exonInTrsptHref=\%exonInTrspt;

$cmd = "grep -P \"\\texon\\t\" " . $inputCdsExonGtf;
$exonLineText = `$cmd`;
@exonLine = split(/\n/, $exonLineText);
foreach $exonLine(@exonLine){
	@field = ();
	@field = split(/\t/, $exonLine);
	@attr = ();
	@attr = split(/; /, $field[8]);
	@tt = ();
	@tt = split(/"/, $attr[1]);
	$trsptId = $tt[1];
	# 将exon追加到trsptId指向的值中
	if(not exists($trsptIdToExonSeries{$trsptId})){
		$trsptIdToExonSeries{$trsptId} = $field[3] . ".." . $field[4];
	}else{
		$trsptIdToExonSeries{$trsptId} .= "," . $field[3] . ".." . $field[4];
	}
	# 建立从exon到trsptId的映射关系
	$exon = $field[0] . "#" . $field[3] . "#" . $field[4] . "#" . $field[6];
	$exonInTrsptHref->{$exon}->{$trsptId}->{"posStart"} = 0;
	$exonInTrsptHref->{$exon}->{$trsptId}->{"posStop"} = 0;
}

# 重新扫描%exonToTrspt中的每个exon，然后找到它对应trsptId，再用trsptId找到对应的exon串，最后计算当前exon在exon串中的开始和结束位置
my ($chr, $exonStart, $exonEnd, $strand, $exonStartInCDS, $exonStopInCDS, @trsptId, $exonInCDSString, $exonSeries);
@exon = ();
@exon = keys(%exonInTrspt);
foreach $exon(@exon){
	($chr, $exonStart, $exonEnd, $strand) = split(/#/, $exon);
	# 获得该exon所在的所有trsptId
	@trsptId = ();
	@trsptId = keys(%{$exonInTrsptHref->{$exon}});
	# 获得每个trspt的exonSeries，计算exon在exonSeries中的位置，并登记到hash中
	foreach $trsptId(@trsptId){
		$exonSeries = $trsptIdToExonSeries{$trsptId};
		# 计算exon在当前trspt的CDS中的位置
		&getExonPosInExonSeries($exonSeries, $exonStart, $exonEnd, \$exonStartInCDS, \$exonStopInCDS);
		# 登记exon在cds中的位置
		$exonInTrsptHref->{$exon}->{$trsptId}->{"posStart"} = $exonStartInCDS;
		$exonInTrsptHref->{$exon}->{$trsptId}->{"posStop"} = $exonStopInCDS;
	}
}

# 读取合并后的exon，将其在所在trspt中的位置输出
open EXON, ">$outputMergedExonNotInUniqPep";
open WW, ">$outputExonPosInCdsSeq";
$cmd = "cat " . $inputExonTsv;
$exonLineText = `$cmd`;
@exonLine = ();
@exonLine = split(/\n/, $exonLineText);
my ($nondetecedMergedExonNum);
foreach $exonLine(@exonLine){
	($chr, $exonStart, $exonEnd, $strand) = split(/\t/, $exonLine);
	$exon = join("#", $chr, $exonStart, $exonEnd, $strand);
	# 检测当前合并外显子是否存在于之前的uniqPep对应的cdsExon中
	if(not exists($exonInTrsptHref->{$exon})){
		print EXON $exonLine . "\n";
		$nondetecedMergedExonNum++;
		next;
	}
	@trsptId = ();
	@trsptId = keys(%{$exonInTrsptHref->{$exon}});
	$exonInCDSString = "";
	foreach $trsptId(@trsptId){
		if($exonInCDSString eq ""){
			$exonInCDSString = $trsptId . ":" . $exonInTrsptHref->{$exon}->{$trsptId}->{"posStart"} . "-" . $exonInTrsptHref->{$exon}->{$trsptId}->{"posStop"};
		}else{
			$exonInCDSString .= "," . $trsptId . ":" . $exonInTrsptHref->{$exon}->{$trsptId}->{"posStart"} . "-" . $exonInTrsptHref->{$exon}->{$trsptId}->{"posStop"};
		}
	}
	print WW join("\t", $exonLine, $exonInCDSString) . "\n";
}
close WW;

print EXON "There are $nondetecedMergedExonNum merged exons not in uniqPep transcripts\n";
close EXON;

sub getExonPosInExonSeries{
	my ($exonSeries, $exonStart, $exonStop, $exonStartInCDS, $exonStopInCDS) = @_;
	my (@exon, $exonId, $i, $j, $exonBegin, $exonEnd);
	@exon = split(/,/, $exonSeries);
	# 检查exon为exonSeries中第几个exon
	for($exonId=0; $exonId<=$#exon; $exonId++){
		if($exonStart . ".." . $exonStop eq $exon[$exonId]){
			last;
		}
	}

	# 开始位置
	$$exonStartInCDS = 0;
	for($i=0; $i<=$exonId-1; $i++){
		($exonBegin, $exonEnd) = split(/\.\./, $exon[$i]);
		$$exonStartInCDS += $exonEnd - $exonBegin + 1;
	}
	$$exonStartInCDS = $$exonStartInCDS + 1;

	# 结束位置
	$$exonStopInCDS = 0;
	for($i=0; $i<=$exonId; $i++){
		($exonBegin, $exonEnd) = split(/\.\./, $exon[$i]);
		$$exonStopInCDS += $exonEnd - $exonBegin + 1;
	}
}
