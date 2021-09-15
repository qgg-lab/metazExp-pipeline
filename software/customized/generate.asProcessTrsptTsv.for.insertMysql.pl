#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--alteredCdnafasta fina.cDNA.fa \\\n" .
                "--alteredPepFasta fina.pep.fa \\\n" .
                "--alteredOrfInCdnaTsv orf.in.cDNA.tsv \\\n" .
                "--alteredGoListTsv  go.uniq.term.list.tsv \\\n" .
                "--alteredPfamListTsv pfam.sorted.list.tsv \\\n" .
                "--origCdnaFasta fina.cDNA.fa \\\n" .
                "--origPepFasta fina.pep.fa \\\n" .
                "--origOrfInCdnaTsv orf.in.cDNA.tsv \\\n" .
                "--origGoListTsv  go.uniq.term.list.tsv \\\n" .
                "--origPfamListTsv pfam.sorted.list.tsv \\\n" .
		"--processTsv as.process.trspt.tsv \\\n" .
		"--asProcessTrsptTsv as.process.trspt.mysql.tsv \n";
	exit;
}

my ($alteredCdnaFasta, $alteredPepFasta, $alteredOrfInCdnaTsv, $alteredGoListTsv, $alteredPfamListTsv, $asToTrsptListTsv, $processTsv, $asProcessTrsptTsv, $origCdnaFasta, $origPepFasta, $origOrfInCdnaTsv, $origGoListTsv, $origPfamListTsv);

GetOptions(
        'alteredCdnaFasta=s'=>\$alteredCdnaFasta,
        'alteredPepFasta=s'=>\$alteredPepFasta,
        'alteredOrfInCdnaTsv=s'=>\$alteredOrfInCdnaTsv,
	'alteredGoListTsv=s'=>\$alteredGoListTsv,
	'alteredPfamListTsv=s'=>\$alteredPfamListTsv,
        'origCdnaFasta=s'=>\$origCdnaFasta,
        'origPepFasta=s'=>\$origPepFasta,
        'origOrfInCdnaTsv=s'=>\$origOrfInCdnaTsv,
	'origGoListTsv=s'=>\$origGoListTsv,
	'origPfamListTsv=s'=>\$origPfamListTsv,
	'processTsv=s'=>\$processTsv,
	'asProcessTrsptTsv=s'=>\$asProcessTrsptTsv,
);

my (%trsptAs, $trsptAsHref, $line);
my (%trspt, $trsptHref);
$trsptAsHref = \%trsptAs;
$trsptHref=\%trspt;

# 打开altered cDNAfasta文件，读取cDNA序列
my (@trsptAsId, $trsptAsId, $asId, $trsptId, @trsptId);
open FF, "<$alteredCdnaFasta";
# >SRX853408.3.6XXXATHAA3SS0000002088
# AGACCCGGACTCTAATTGCTCCGTATTCTTCTTCTCTTGAGAGAG 
while($line = <FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$trsptAsId = $1;
		@trsptAsId = ();
		@trsptAsId = split(/ /, $trsptAsId);
		$trsptAsId = $trsptAsId[0];
	}else{
		$trsptAsHref->{$trsptAsId}->{"cDNAseq"}.=$line;
	}
}
close FF;

# orig cdnaFasta
open FF, "<$origCdnaFasta";
# >SRX853408.3.6
# AGACCCGGACTCTAATTGCTCCGTATTCTTCTTCTCTTGAGAGAG 
while($line = <FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$trsptId = $1;
		@trsptId = ();
		@trsptId = split(/ /, $trsptId);
		$trsptId = $trsptId[0];
	}else{
		$trsptHref->{$trsptId}->{"cDNAseq"}.=$line;
	}
}
close FF;



############### pep fasta ###########################
#
# altered pep fata
# >AT5G51710.2XXXATHASE0000033759
# MVQVETVAQFGVVFLLFALGLEFSMTKLKVVGPVAVLG
open FF, "<$alteredPepFasta";
while($line = <FF>){
	chomp($line);
	if($line=~/>(.*)/){
                $trsptAsId = $1;
                @trsptAsId = ();
                @trsptAsId = split(/ /, $trsptAsId);
                $trsptAsId = $trsptAsId[0];
	}else{
		$trsptAsHref->{$trsptAsId}->{"pepSeq"}.=$line
	}
}
close FF;

# orig pepFasta
open FF, "<$origPepFasta";
while($line=<FF>){
	chomp($line);
        if($line=~/>(.*)/){
                $trsptId = $1;
                @trsptId = ();
                @trsptId = split(/ /, $trsptId);
                $trsptId = $trsptId[0];
        }else{
                $trsptHref->{$trsptId}->{"pepSeq"}.=$line
        }
}
close FF;



################## orf ##############################
my ($i, $j, @titleField, $titleField, @valueField, $valueField);

# altered orf
open FF, "<$alteredOrfInCdnaTsv";
# trsptId orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonEnd
# AT5G51710.2XXXATHASE0000033759  213     1907    213     215     1908    1910
$line = <FF>;
chomp($line);
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	$trsptAsId = $valueField[0];
	for($i=0; $i<=$#valueField; $i++){
		$trsptAsHref->{$trsptAsId}->{$titleField[$i]} = $valueField[$i];
	}
}
close FF;

# orig orf
open FF, "<$origOrfInCdnaTsv";
# trsptId orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonEnd
# AT3G53400.1     420     1817    420     422     1818    1820
# AT5G45110.2     515     2062    515     517     2063    2065
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	$trsptId = $valueField[0];
	for($i=0; $i<=$#valueField; $i++){
		$trsptHref->{$trsptId}->{$titleField[$i]} = $valueField[$i];
	}
}
close FF;




###################### GO term list  ############################
# altered go term list
my (@field);
open FF, "<$alteredGoListTsv";
# AT5G51710.2XXXATHASE0000033759  GO:0006812|GO:0015299|GO:0016021|GO:0055085
# SRX399565.17511.2XXXATHASE0000041372    GO:0004609|GO:0008654
while($line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);	
	$trsptAsId = $field[0];
	$trsptAsHref->{$trsptAsId}->{"goList"} = $field[1];
}
close FF;

# orig gotermlist
open FF, "<$origGoListTsv";
# AT5G51710.2  GO:0006812|GO:0015299|GO:0016021|GO:0055085
# SRX399565.17511.2    GO:0004609|GO:0008654
while($line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);	
	$trsptId = $field[0];
	$trsptHref->{$trsptId}->{"goList"} = $field[1];
}
close FF;



############### pfam list #############################

# altered pfam list
open FF, "<$alteredPfamListTsv";
# AT2G06005.1XXXATHASE0000043648  PF14802[491,991]
while($line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);	
	$trsptAsId = $field[0];
	$trsptAsHref->{$trsptAsId}->{"pfamList"} = $field[1];
}
close FF;

# orig pfam list
open FF, "<$origPfamListTsv";
# AT2G06005.1  PF14802[491,991]
# AT2G01690.1  PF12755[566,856]|PF11916[1658,2197]
while($line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);	
	$trsptId = $field[0];
	$trsptHref->{$trsptId}->{"pfamList"} = $field[1];
}
close FF;



### 读取as对trspt的加工过程 #######
open WW, ">$asProcessTrsptTsv";
# 打开process文件逐行输出
# chr     strand  asId    asType  trsptId residentType    editSite        cutSize insertSize      insertSeq
# 1       -       ATHAA3SS0000002088      A3SS    SRX853408.3.6   longAlt 305     152     0
my ($chr, $strand, $asId, $asType, $trsptId, $residentType, $editSite, $cutSize, $insertSize, $insertSeq, $trsptAsId, $goImpact, $pfamImpact, $hitPfamIdList, $frameImpact, $stopCodonImpact);
my ($fieldString, $valueString);
open FF, "<$processTsv";
<FF>;
@titleField = split(/\t/, $line);
while($line = <FF>){
	chomp($line);
	$insertSeq = "";
	($chr, $strand, $asId, $asType, $trsptId, $residentType, $editSite, $cutSize, $insertSize, $insertSeq) = split(/\t/, $line);

	if($insertSeq eq ""){
		$insertSeq = "-";
	}
	
	# 生成转录本和可变剪接的组合ID
	$trsptAsId = $trsptId . "XXX" . $asId;

	##### pep seq #####
	# altered pep
	if(not exists($trsptAsHref->{$trsptAsId}->{"pepSeq"})){
		$trsptAsHref->{$trsptAsId}->{"pepSeq"} = "-";
	}
	# orig pep 
	if(not exists($trsptHref->{$trsptId}->{"pepSeq"})){
		$trsptHref->{$trsptId}->{"pepSeq"} = "-";
	}

	##### go list #######
	# altered go list
	if(not exists($trsptAsHref->{$trsptAsId}->{"goList"}) or $trsptAsHref->{$trsptAsId}->{"goList"} eq ""){
		$trsptAsHref->{$trsptAsId}->{"goList"} = "NA";
	}
	# orig go list
	if(not exists($trsptHref->{$trsptId}->{"goList"}) or $trsptHref->{$trsptId}->{"goList"} eq ""){
		$trsptHref->{$trsptId}->{"goList"} = "NA";
	}
	# 判断AS编辑前后go是否发生改变
	if($trsptAsHref->{$trsptAsId}->{"goList"} eq "NA" and $trsptHref->{$trsptId}->{"goList"} eq "NA"){
		$goImpact = "N";
	}elsif($trsptAsHref->{$trsptAsId}->{"goList"} ne "NA" and $trsptHref->{$trsptId}->{"goList"} eq "NA"){
		$goImpact = "Y";
	}elsif($trsptAsHref->{$trsptAsId}->{"goList"} eq "NA" and $trsptHref->{$trsptId}->{"goList"} ne "NA"){
		$goImpact = "Y";
	}elsif($trsptAsHref->{$trsptAsId}->{"goList"} ne "NA" and $trsptHref->{$trsptId}->{"goList"} ne "NA"){
		# 对goList排序(不去冗余)，然后再比较两者是否相同
		$goImpact = &getGoImpact($trsptAsHref->{$trsptAsId}->{"goList"}, $trsptHref->{$trsptId}->{"goList"});
	}

	##### pfam list #####
	# altered pfam list
	if(not exists($trsptAsHref->{$trsptAsId}->{"pfamList"}) or $trsptAsHref->{$trsptAsId}->{"pfamList"} eq ""){
		$trsptAsHref->{$trsptAsId}->{"pfamList"} = "NA";
	}
	if(not exists($trsptHref->{$trsptId}->{"pfamList"}) or $trsptHref->{$trsptId}->{"pfamList"} eq ""){
		$trsptHref->{$trsptId}->{"pfamList"} = "NA";
	}
	# 判断AS编辑前后pfam是否发生改变
	if($trsptAsHref->{$trsptAsId}->{"pfamList"} eq "NA" and $trsptHref->{$trsptId}->{"pfamList"} eq "NA"){
		$pfamImpact = "N";
	}elsif($trsptAsHref->{$trsptAsId}->{"pfamList"} ne "NA" and $trsptHref->{$trsptId}->{"pfamList"} eq "NA"){
		$pfamImpact = "Y";
	}elsif($trsptAsHref->{$trsptAsId}->{"pfamList"} eq "NA" and $trsptHref->{$trsptId}->{"pfamList"} ne "NA"){
		$pfamImpact = "Y";
	}elsif($trsptAsHref->{$trsptAsId}->{"pfamList"} ne "NA" and $trsptHref->{$trsptId}->{"pfamList"} ne "NA"){
		$pfamImpact = &getPfamImpact($trsptAsHref->{$trsptAsId}->{"pfamList"}, $trsptHref->{$trsptId}->{"pfamList"});
	}

	# 判断AS对stopCodon的影响
	# premature:提前成熟，可能会导致NMD
	# postmature：退后成熟
	# constant：3UTR不变化
	# bothNo: AS作用前后都是无终止密码子
	# create：创建了终止密码子, lose：丢失了终止密码子
	# undetect：未检测到
	if($trsptHref->{$trsptId}->{"stopCodonBegin"}==-1){
		# AS作用前转录本无终止密码子
		if($trsptAsHref->{$trsptAsId}->{"stopCodonBegin"} == -1){
			# AS作用后转录本仍然无终止密码子
			$stopCodonImpact = "bothNo";
		}elsif($trsptAsHref->{$trsptAsId}->{"stopCodonBegin"} >= 0){
			# AS作用后创建了终止密码子
			$stopCodonImpact = "create";
		}
	}elsif($trsptAsHref->{$trsptAsId}->{"stopCodonBegin"}==-1){
		# AS作用前转录本有终止密码子，作用后无密码子
		$stopCodonImpact = "lose";
	}else{
		# AS作用前转录本有终止密码子，作用后也有密码子
		if(length($trsptHref->{$trsptId}->{"cDNAseq"}) - $trsptHref->{$trsptId}->{"stopCodonEnd"} == length($trsptAsHref->{$trsptAsId}->{"cDNAseq"}) - $trsptAsHref->{$trsptAsId}->{"stopCodonEnd"}){
			# 作用前后3UTR长度相同
			$stopCodonImpact = "constant";
		}elsif(length($trsptHref->{$trsptId}->{"cDNAseq"}) - $trsptHref->{$trsptId}->{"stopCodonEnd"} < length($trsptAsHref->{$trsptAsId}->{"cDNAseq"}) - $trsptAsHref->{$trsptAsId}->{"stopCodonEnd"}){
			# 相对于cDNA末端而言，AS作用后的3UTR变长，即终止密码子提前成熟
			$stopCodonImpact = "premature";
		}elsif(length($trsptHref->{$trsptId}->{"cDNAseq"}) - $trsptHref->{$trsptId}->{"stopCodonEnd"} > length($trsptAsHref->{$trsptAsId}->{"cDNAseq"}) - $trsptAsHref->{$trsptAsId}->{"stopCodonEnd"}){
			# 相对于cDNA末端而言，AS作用后的3UTR变短，即终止码子推迟成熟(因为相位改变导致原终止密码子失效)
			$stopCodonImpact = "postmature";
		}
	}


	# 转录本和AS关系
	if($residentType eq "longAlt"){
		$residentType = "inclusion";
		if($asType eq "MXE"){
			$residentType = "1st";
		}
	}else{
		$residentType = "exclusion";
		if($asType eq "MXE"){
			$residentType = "2nd";
		}
	}

	# 获得AS直接命中了哪些pfam：从编辑点开始 到 cutsize结束位置，这个范围涉及到哪些pfam
	if($trsptHref->{$trsptId}->{"pfamList"} eq "NA"){
		$hitPfamIdList = "NA";
	}else{
		$hitPfamIdList = &getHitPfamIdList($trsptHref->{$trsptId}->{"pfamList"}, $residentType, $editSite, $cutSize);
	}

	$frameImpact = "NA";
	# 判断AS是否对orf frame产生影响：cutSize + insertSize 是否为3的整数倍
	$frameImpact = &getFrameImpact($residentType, $editSite, $cutSize, $insertSize, $trsptHref->{$trsptId}->{"orfBegin"}, $trsptHref->{$trsptId}->{"orfEnd"});


##########  fieldString  ########
	$fieldString = join(", ", 
		"asId", 
		"trsptId", 
		"residentType", 
		"asType", 
		"editSiteInOrigCdna", 
		"cutSize", 
		"insertSize",
		"origOrfStart",
		"origOrfStop",
		"origStartCodonBegin",
		"origStartCodonEnd",
		"origStopCodonBegin",
		"origStopCodonEnd",
		"newOrfStart", 
		"newOrfStop", 
		"newStartCodonBegin", 
		"newStartCodonEnd", 
		"newStopCodonBegin", 
		"newStopCodonEnd", 
		"origGoTermList",
		"origPfamList",
		"newGoTermList", 
		"newPfamList",
		"goImpact",
		"pfamImpact",
		"hitPfamList",
		"frameImpact",
		"stopCodonImpact",
		"origCdnaSeq",
		"origPepSeq",
		"newCdnaSeq", 
		"newPepSeq",
		"insertSeq",
);

########  valueString ##############
	# orfBegin        orfEnd  startCodonBegin startCodonEnd   stopCodonBegin  stopCodonEnd
	$valueString = join(", ", 
		$asId, 
		$trsptId, 
		$residentType, 
		$asType, 
		$editSite, 
		$cutSize,
		$insertSize, 		
		$trsptHref->{$trsptId}->{"orfBegin"}, 
		$trsptHref->{$trsptId}->{"orfEnd"},
		$trsptHref->{$trsptId}->{"startCodonBegin"}, 
		$trsptHref->{$trsptId}->{"startCodonEnd"},
		$trsptHref->{$trsptId}->{"stopCodonBegin"}, 
		$trsptHref->{$trsptId}->{"stopCodonEnd"},
		$trsptAsHref->{$trsptAsId}->{"orfBegin"}, 
		$trsptAsHref->{$trsptAsId}->{"orfEnd"},
		$trsptAsHref->{$trsptAsId}->{"startCodonBegin"}, 
		$trsptAsHref->{$trsptAsId}->{"startCodonEnd"},
		$trsptAsHref->{$trsptAsId}->{"stopCodonBegin"}, 
		$trsptAsHref->{$trsptAsId}->{"stopCodonEnd"},
		$trsptHref->{$trsptId}->{"goList"}, 
		$trsptHref->{$trsptId}->{"pfamList"},
		$trsptAsHref->{$trsptAsId}->{"goList"}, 
		$trsptAsHref->{$trsptAsId}->{"pfamList"},
		$goImpact,
		$pfamImpact,
		$hitPfamIdList,
		$frameImpact,
		$stopCodonImpact,
		$trsptHref->{$trsptId}->{"cDNAseq"}, 
		$trsptHref->{$trsptId}->{"pepSeq"},
		$trsptAsHref->{$trsptAsId}->{"cDNAseq"}, 
		$trsptAsHref->{$trsptAsId}->{"pepSeq"},
		$insertSeq);
	print WW $fieldString . "___" . $valueString . "\n";
}
close WW;

# GO:0005509|GO:0009523|GO:0009654|GO:0015979|GO:0019898
sub getGoImpact{
	my ($origGoList, $alteredGoList) = @_;
	my (@origGo, @alteredGo);
	@origGo = split(/\|/, $origGoList);
	@alteredGo = split(/\|/, $alteredGoList);
	@origGo = sort(@origGo);
	@alteredGo = sort(@alteredGo);
	if($#origGo!=$#alteredGo){
		return "Y";
	}else{
		for(my $i=0; $i<=$#origGo; $i++){
			return "Y" if($origGo[$i] ne $alteredGo[$i]);
		}
	}
	return "N";
}

# PF15044[385,597]|PF12807[2344,2781]|PF13424[2989,3201]|PF13424[3241,3465]
sub getPfamImpact{
	my ($origPfamList, $alteredPfamList) = @_;
	my (@origPfam, @alteredPfam, $origPfamId, $alteredPfamId, $otherPart);
	@origPfam = split(/\|/, $origPfamList);
	@alteredPfam = split(/\|/, $alteredPfamList);
	if($#origPfam!=$#alteredPfam){
		return "Y";
	}else{
		for(my $i=0; $i<=$#origPfam; $i++){
			($origPfamId, $otherPart) = split(/\[/, $origPfam[$i]);
			($alteredPfamId, $otherPart) = split(/\[/, $alteredPfam[$i]);
			return "Y" if($origPfamId ne $alteredPfamId);
		}
	}
	return "N";
}


# 根据AS编辑切掉的位置，判断命中了那些pfamId
# PF04161[483,560]|PF04161[561,995], 
sub getHitPfamIdList{
	my($pfamPosList, $residentType, $editSite, $cutSize) = @_;
	my ($pfamIdList);
	my ($cutBeg, $cutEnd, $insertSite, $i);
	my (@pfamPos, $pfamPos, @pfamPosArr, $pfamNum);
	$pfamIdList = "";
	# 提取pfamId和位置到数组中
	@pfamPos=split(/\|/, $pfamPosList);
	foreach $pfamPos(@pfamPos){
		if($pfamPos=~/(PF\d+)\[(\d+),(\d+)\]/){
			($pfamPosArr[$pfamNum][0], $pfamPosArr[$pfamNum][1], $pfamPosArr[$pfamNum][2]) = ($1, $2, $3);
			$pfamNum++;
		}
	}

	if($residentType eq "inclusion" or $residentType eq "1st" or $residentType eq "2nd"){
		$cutBeg = $editSite + 1;
		$cutEnd = $editSite + $cutSize;
		for($i=0; $i<$pfamNum; $i++){
			if(not($cutBeg > $pfamPosArr[$i][2] or $cutEnd < $pfamPosArr[$i][1])){
				$pfamIdList .= $pfamPosArr[$i][0] . ",";
			}
		}
		if($pfamIdList ne ""){
			return substr($pfamIdList, 0, length($pfamIdList) - 1);
		}else{
			return "NA";
		}
	}elsif($residentType eq "exclusion"){
		$insertSite = $editSite + 1;
		for($i=0; $i<$pfamNum; $i++){
			if($insertSite<= $pfamPosArr[$i][2] and $insertSite>= $pfamPosArr[$i][1]){
				$pfamIdList .= $pfamPosArr[$i][0] . ",";
			}
		}
		if($pfamIdList ne ""){
			return substr($pfamIdList, 0, length($pfamIdList) - 1);
		}else{
			return "NA";
		}
	}
}



# 判断对frame的影响: Y:有影响；N:无影响；D:整个orf被删除
sub getFrameImpact{
	my ($residentType, $editSite, $cutSize, $insertSize, $orfBeg, $orfEnd) = @_;
	my ($frameImpact);
	my ($cutBeg, $cutEnd, $insertPos, $cutOrfLen, $changedSize);
	# 转录本和AS之间是inclusion关系时，需要从转录本中删除掉一段序列
	if($residentType eq "inclusion"){
		$cutBeg = $editSite + 1;
		$cutEnd = $editSite + $cutSize;
		#需要删除一段序列
		if($cutEnd < $orfBeg){
			# 删除序列在5UTR
			$cutOrfLen = 0;
			$frameImpact = "N_cut_5UTR";
		}elsif($cutBeg > $orfEnd){
			# 删除序列在3UTR
			$frameImpact = "N_cut_3UTR";
		}elsif($cutBeg >= $orfBeg and $cutEnd <= $orfEnd){
			# 删除序列完全在orf内部
			$cutOrfLen = $cutSize;
			if(int($cutOrfLen/3)*3 - $cutOrfLen == 0){
				# 删除碱基长度为3的整数倍，那么不会改变删除序列之外序列的frame
				$frameImpact = "N_cut_CDSInner";	
			}else{
				# 删除碱基长度不为3的整数倍，那么会改变被删除序列之后的frame，导致之后序列的编码发生大规模改变
				$frameImpact = "Y_cut_CDSInner";
			}
		}elsif($cutBeg <= $orfBeg and $cutEnd >= $orfEnd){
			# orf被完整删除掉，D表示orf被全部删除掉，frame不再存在
			$cutOrfLen = $orfEnd - $orfBeg + 1;
			$frameImpact = "Y_cut_CDSEntire";
		}elsif($cutBeg < $orfBeg and $cutEnd >= $orfBeg and $cutEnd<= $orfEnd){
			# 删除区域和orf左侧有重叠，也就是在orf左侧前导序列被删除
			# 由于CDS的起始密码子丢失，因此无法判断frame是否改变
			$frameImpact = "UN_cut_5UTRCDS";
		}elsif($cutBeg <= $orfEnd and $cutBeg >= $orfBeg and $cutEnd >= $orfEnd){
			# 删除区域和orf右侧有重叠，也就是orf右侧序列被删除，包含终止密码子被删除
			# 删除区域之后的orf区域不存在了，因此frame
			$frameImpact = "UN_cut_CDS3UTR";
		}
	}elsif($residentType eq "exclusion"){
		$insertPos = $editSite;
		#需要插入一段序列
		if($insertPos < $orfBeg){
			# 插入位置在orf的左侧，即5UTR内
			$frameImpact = "N_insert_5UTR";
		}elsif($insertPos > $orfEnd){
			# 插入位置orf的右侧，即在3UTR内
			$frameImpact = "N_insert_3UTR";
		}else{
			# 插入位置在orf内
			if(int($insertSize/3)*3-$insertSize == 0){
				# 插入的长度为3的整数倍，不改变frame
				$frameImpact = "N_insert_CDSInner";
			}else{
				# 插入长度不为3的整数倍，改变frame
				$frameImpact = "Y_insert_CDSInner";
			}
		}
	}elsif($residentType eq "1st" or $residentType eq "2nd"){
		#先删掉一段序列，然后再插入一段序列
                $cutBeg = $editSite + 1;
                $cutEnd = $editSite + $cutSize;
		if($cutEnd < $orfBeg){
			# 替换序列位置在5UTR
			$frameImpact = "N_replace_5UTR";
		}elsif($cutBeg > $orfEnd){
			# 替换序列位置在3UTR
			$frameImpact = "N_replace_3UTR";
		}elsif($cutBeg >= $orfBeg and $cutEnd <= $orfEnd){
			# 替换序列完全在CDS内部
			$changedSize = $insertSize - $cutSize;
			if($changedSize < 0){
				$changedSize = -1 * $changedSize;
			}
			if(int($changedSize/3)*3 - $changedSize == 0){
				# 替换长度为3的整数倍
				$frameImpact = "N_replace_CDSInner";
			}else{
				# 替换长度不为3的整数倍
				$frameImpact = "Y_replace_CDSInner";
			}
		}elsif($cutBeg <= $orfBeg and $cutEnd >= $orfEnd){
			# orf被完全替换
			$frameImpact = "Y_replace_CDSEntire";
		}elsif($cutBeg < $orfBeg and $cutEnd >= $orfBeg and $cutEnd<= $orfEnd){
			# 替换序列覆盖5UTR和CDS
			# 由于起始密码子丢失，因此无法确定frame是否改变
			$frameImpact = "UN_replace_5UTRCDS";
		}elsif($cutBeg <= $orfEnd and $cutBeg >= $orfBeg and $cutEnd >= $orfEnd){
			# 替换终止密码子被替换
			$frameImpact = "UN_replace_CDS3UTR";
		}
	}
	return $frameImpact;
}

# 判断对stopCodon影响: 提前，推后，消失
sub getStopCodonImpact{
	my ($residentType, $editSite, $cutSize, $insertSize, $origCdnaSeq, $orfBeg, $orfEnd, $stopCodonBegin, $stopCodonEnd);
	if($residentType eq "inclusion"){
		#需要删除一段序列
	}elsif($residentType eq "exclusion"){
		#需要插入一段序列
	}elsif($residentType eq "1st" or $residentType eq "2nd"){
		#先删除后插入
	}
}
