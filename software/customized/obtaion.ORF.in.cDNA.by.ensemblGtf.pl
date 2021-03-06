#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--ensemblGtf ensembl.gtf \\\n" .
		"--cDNAseqFile cDNAseqFile \\\n" .
                "--orfInCdnaTsv orf.in.ensembl.trspt.cDNA.tsv \n";
	exit;
}

my ($ensemblGtf, $orfInCdnaTsv, $cDNAseqFile);

GetOptions(
        'ensemblGtf=s'=>\$ensemblGtf,
	'cDNAseqFile=s'=>\$cDNAseqFile,
        'orfInCdnaTsv=s'=>\$orfInCdnaTsv,
);

my (%cDNAseq);
my ($line, %trsptCoordinate, $trsptId, $trsptCoordinateHref, @tmp, @trsptId, @exon, @exonString, $exonString, $exonNum);
my ($startCodonBegin, $startCodonEnd, $stopCodonBegin, $stopCodonEnd);
$trsptCoordinateHref = \%trsptCoordinate;

# 读取cDNAseq文件，获得cDNAseq的长度
open FF, "<$cDNAseqFile";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$trsptId = $1;
		@trsptId = split(/ /, $trsptId);
		$trsptId = $trsptId[0];
	}else{
		$cDNAseq{$trsptId}.=$line;
	}
}
close FF;


# 读取gtf中转录本的外显子串和起始终止密码子位置
open FF, "<$ensemblGtf";
open WW, ">$orfInCdnaTsv";
print WW join("\t", "trsptId", "orfBegin", "orfEnd", "startCodonBegin", "startCodonEnd", "stopCodonBegin", "stopCodonEnd", "cDNAseqLen") . "\n";
while($line=<FF>){
	chomp($line);
	@tmp = ();
	@tmp = split(/\t/, $line);

	$trsptId = &getTrsptIdInAttr($tmp[8]);
	$trsptCoordinateHref->{$trsptId}->{"strand"} = $tmp[6];

	if($tmp[2] eq "CDS"){
		$trsptCoordinateHref->{$trsptId}->{"cdsExon"} .= $tmp[3] . "-" . $tmp[4] . ",";
		$trsptCoordinateHref->{$trsptId}->{"coding"} = "Y";
	}elsif($tmp[2] eq "exon"){
		$trsptCoordinateHref->{$trsptId}->{"exon"} .= $tmp[3] . "-" . $tmp[4] . ",";
		$trsptCoordinateHref->{$trsptId}->{"cDNALen"} += $tmp[4] - $tmp[3] +1;
	}elsif($tmp[2] eq "start_codon"){
		$trsptCoordinateHref->{$trsptId}->{"startCodonBegin"} .= $tmp[3] . ",";
		$trsptCoordinateHref->{$trsptId}->{"startCodonEnd"} .= $tmp[4] . ",";
	}elsif($tmp[2] eq "stop_codon"){
		$trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"} .= $tmp[3] . ",";
		$trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"} .= $tmp[4] . ",";
	}
}
close FF;

# 重新扫描每个trspt，处理存在orf但起始密码子和终止密码子不同时存在
# 处理方法，将缺失的startCodonBegin, startCodonEnd 或 stopCodonBegin, stopCodonEnd同时设置为-1
@trsptId = keys(%trsptCoordinate);
foreach $trsptId(@trsptId){
	# 不是编码基因则放弃！
	if(not exists($trsptCoordinateHref->{$trsptId}->{"coding"})){
		next;
	}
	if(not exists($trsptCoordinateHref->{$trsptId}->{"startCodonBegin"}) and exists($trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"})){
		$trsptCoordinateHref->{$trsptId}->{"startCodonBegin"} = -1;
		$trsptCoordinateHref->{$trsptId}->{"startCodonEnd"} = -1;
	}
	if(exists($trsptCoordinateHref->{$trsptId}->{"startCodonBegin"}) and not exists($trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"})){
		$trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"} = -1;
		$trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"} = -1;
	}
	if(not exists($trsptCoordinateHref->{$trsptId}->{"startCodonBegin"}) and not exists($trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"})){
		$trsptCoordinateHref->{$trsptId}->{"startCodonBegin"} = -1;
		$trsptCoordinateHref->{$trsptId}->{"startCodonEnd"} = -1;
		$trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"} = -1;
		$trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"} = -1;
	}
#	print join("\t", $trsptCoordinateHref->{$trsptId}->{"startCodonBegin"}, $trsptCoordinateHref->{$trsptId}->{"startCodonEnd"}, $trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"}, $trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"});
#	<STDIN>;
}

# 重新扫描每个trspt，然后计算起始密码子和终止密码子在cDNA中的位置
foreach $trsptId(@trsptId){

	# 对于不存在CDS的转录本，给予放弃
	next if($trsptCoordinateHref->{$trsptId}->{"coding"} ne "Y");
	# 求出起始密码子和终止密码子的开始和结束位置（防止出现密码子被split)
	if($trsptCoordinateHref->{$trsptId}->{"strand"} eq "+"){
		@tmp = ();
		@tmp = split(/,/, $trsptCoordinateHref->{$trsptId}->{"startCodonBegin"});
		$trsptCoordinateHref->{$trsptId}->{"startCodonBegin"} = $tmp[0];
		@tmp = ();
		@tmp = split(/,/, $trsptCoordinateHref->{$trsptId}->{"startCodonEnd"});
		$trsptCoordinateHref->{$trsptId}->{"startCodonEnd"} = $tmp[$#tmp];
	}else{
		my $tmpStartCodonBeginString = $trsptCoordinateHref->{$trsptId}->{"startCodonBegin"};
		my $tmpStartCodonEndString = $trsptCoordinateHref->{$trsptId}->{"startCodonEnd"};
		@tmp = ();
		@tmp = split(/,/, $tmpStartCodonEndString);
		$trsptCoordinateHref->{$trsptId}->{"startCodonBegin"} = $tmp[0];
		@tmp = ();
		@tmp = split(/,/, $tmpStartCodonBeginString);
		$trsptCoordinateHref->{$trsptId}->{"startCodonEnd"} = $tmp[$#tmp];

	}
	# 终止密码子被split时，种子密码子开始为开始坐标中的第1个数组元素，终止为结束坐标中的最后1个元素
	if($trsptCoordinateHref->{$trsptId}->{"strand"} eq "+"){
		@tmp = ();
		@tmp = split(/,/, $trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"});
		$trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"} = $tmp[0];
		@tmp = ();
		@tmp = split(/,/, $trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"});
		$trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"} = $tmp[$#tmp];
	}else{
		my $tmpStopCodonBeginString = $trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"};
		my $tmpStopCodonEndString = $trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"};
		@tmp = ();
		@tmp = split(/,/, $tmpStopCodonEndString);
		$trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"} = $tmp[0];
		@tmp = ();
		@tmp = split(/,/, $tmpStopCodonBeginString);
		$trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"} = $tmp[$#tmp];
	}


	# 将外显子坐标放到二维数组exon中
	@exon = ();
	@exonString = ();
	@exonString = split(/,/, $trsptCoordinateHref->{$trsptId}->{"exon"});
	$exonNum = 0;
	foreach $exonString(@exonString){
		@tmp = split(/\-/, $exonString);
		$exon[$exonNum][0] = $tmp[0];
		$exon[$exonNum][1] = $tmp[1];
		$exonNum = $exonNum+1;
	}

	
	# 开始计算起始密码子和终止密码子在cDNA中的开始和结束位置
	if($trsptCoordinateHref->{$trsptId}->{"startCodonBegin"}!=-1){
		$trsptCoordinateHref->{$trsptId}->{"startCodonBegin"} = &getCoordinateInCdna(\@exon, $trsptCoordinateHref->{$trsptId}->{"strand"}, $trsptCoordinateHref->{$trsptId}->{"startCodonBegin"});
	}
	# 计算起始密码子在cDNA中结束位置
	if($trsptCoordinateHref->{$trsptId}->{"startCodonEnd"}!=-1){
		$trsptCoordinateHref->{$trsptId}->{"startCodonEnd"} = &getCoordinateInCdna(\@exon, $trsptCoordinateHref->{$trsptId}->{"strand"}, $trsptCoordinateHref->{$trsptId}->{"startCodonEnd"});
	}

	# 计算终止密码子在cDNA中位置
	# 如果有终止密码子，那么orfEnd为终止码子前1位
	# 如果没有终止密码子，那么orfEnd为cDNA的长度
	if($trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"} !=-1){
		$trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"} = &getCoordinateInCdna(\@exon, $trsptCoordinateHref->{$trsptId}->{"strand"}, $trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"});
	}
	# 计算终止密码子的结束位置
	if($trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"} !=-1){
		$trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"} = &getCoordinateInCdna(\@exon, $trsptCoordinateHref->{$trsptId}->{"strand"}, $trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"});
	}


	# 计算orfBegin和orfEnd
	my ($cdsBeginCoord, $cdsEndCoord);
	&getCdsBeginAndEndCoord($trsptCoordinateHref->{$trsptId}->{"cdsExon"}, $trsptCoordinateHref->{$trsptId}->{"strand"}, \$cdsBeginCoord, \$cdsEndCoord);
#	print join("\t", $trsptId, $cdsBeginCoord, $cdsEndCoord);
#	<STDIN>;
	$trsptCoordinateHref->{$trsptId}->{"orfBegin"} = &getCoordinateInCdna(\@exon, $trsptCoordinateHref->{$trsptId}->{"strand"}, $cdsBeginCoord);
	$trsptCoordinateHref->{$trsptId}->{"orfEnd"} = &getCoordinateInCdna(\@exon, $trsptCoordinateHref->{$trsptId}->{"strand"}, $cdsEndCoord);

	# 输出orf在cDNA中的开始和结束位置，以及起始密码子和终止密码子在cDNA中开始和结束位置
	print WW join("\t", $trsptId, $trsptCoordinateHref->{$trsptId}->{"orfBegin"}, $trsptCoordinateHref->{$trsptId}->{"orfEnd"}, $trsptCoordinateHref->{$trsptId}->{"startCodonBegin"}, $trsptCoordinateHref->{$trsptId}->{"startCodonEnd"}, $trsptCoordinateHref->{$trsptId}->{"stopCodonBegin"}, $trsptCoordinateHref->{$trsptId}->{"stopCodonEnd"}, length($cDNAseq{$trsptId})) . "\n";
}

sub getTrsptIdInAttr{
	my ($attrString)=@_;
	# gene_id "Zm00001d027231"; transcript_id "Zm00001d027231_T002";
	my ($attr, @attr);
	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/transcript_id "(.*?)"/){
			return $1;
		}
	}
	return "";
}

# 计算一个基因组位置在由外显子坐标构成二维数组中的相对位置
sub getCoordinateInCdna{
	my ($exonInDna, $strand, $siteInDna) = @_;
	my ($exonNum, $position, @tt, $i);
	$exonNum = $#$exonInDna;
	if($strand eq "+"){
		for($i=0; $i<=$exonNum; $i++){
			if($siteInDna > $$exonInDna[$i][1]){
				$position+=$$exonInDna[$i][1] - $$exonInDna[$i][0] + 1;
			}else{
				$position+= $siteInDna - $$exonInDna[$i][0] + 1;
				return $position;
			}	
		}
	}else{
		for($i=0; $i<=$exonNum; $i++){
			if($siteInDna < $$exonInDna[$i][0]){
				$position+=$$exonInDna[$i][1] - $$exonInDna[$i][0] + 1;
			}else{
				$position+= $$exonInDna[$i][1] - $siteInDna + 1;
				return $position;
			}	
		}
	}
}

# 根据cdsExon串计算cds的begin和end的genome coordinate
sub getCdsBeginAndEndCoord{
	my ($cdsExonSeries, $strand, $cdsBeginCoord, $cdsEndCoord) = @_;
	my (@cdsExon, $cdsExon, @coor, $coor);
	$$cdsBeginCoord = 1000000000;
	$$cdsEndCoord = 0;
	@cdsExon = split(/,/, $cdsExonSeries);
	foreach $cdsExon(@cdsExon){
		@coor = split(/\-/, $cdsExon);
		foreach $coor(@coor){
			if($$cdsBeginCoord > $coor){
				$$cdsBeginCoord = $coor;
			}
			if($$cdsEndCoord < $coor ){
				$$cdsEndCoord = $coor;
			}
		}
	}
	if($strand eq "-"){
		($$cdsBeginCoord, $$cdsEndCoord) = ($$cdsEndCoord, $$cdsBeginCoord);
	}
}
