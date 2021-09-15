#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtf \\\n" .
                "--pepFasta \\\n" .
                "--uniqPepFasta \\\n" . 
		"--checkPepIdListTsv \n";
	exit;
}

my ($gtf, $pepFasta, $uniqPepFasta, $checkPepIdListTsv);

GetOptions(
        'gtf=s'=>\$gtf,
        'pepFasta=s'=>\$pepFasta,
        'uniqPepFasta=s'=>\$uniqPepFasta,
	'checkPepIdListTsv=s'=>\$checkPepIdListTsv,
);

my (%trsptIdToGeneId, %trsptIdToPepSeq, %pepSeqToGeneTrsptIdList, @trsptId, $hrefPepSeqToGeneTrsptIdList);
my ($trsptId, $geneId, @tt, $line, @field, $pepSeq);
my (%geneId, %trsptId);
$hrefPepSeqToGeneTrsptIdList = \%pepSeqToGeneTrsptIdList;

# 读取gtf获得trsptId 和  geneId 之间映射关系
open FF, "<$gtf";
while($line=<FF>){
	@field = ();
	@field = split(/\t/, $line);
	next if($field[2] ne "transcript");
#	if($field[8]=~/gene_id "(.*?)"; transcript_id "(.*?)"; /){
#		$trsptIdToGeneId{$2} = $1;
#	}
	($geneId, $trsptId) = ("", "");
	&getGeneIdAndTrsptId($field[8], \$geneId, \$trsptId);
	$trsptIdToGeneId{$trsptId} = $geneId;
}
close FF;

# 读取第pep文件，将geneId和trscriptId关联到pep序列上
# 原始转录本序列格式
# >SRX399565.4024.3
# AS作用后的转录本序列编号格式：
# >AT5G51710.2XXXATHASE0000033759
open FF, "<$pepFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$trsptId = $1;
		@trsptId = split(/ /, $trsptId);
		$trsptId = $trsptId[0];
	}else{
		$trsptIdToPepSeq{$trsptId}.=$line;
	}
}
close FF;

# 以pepSeq为关键字，将trsptId和geneId分别分配到以pepSeq为关键字的hash中
my ($origTrsptId);
@trsptId = keys(%trsptIdToPepSeq);
foreach $trsptId(@trsptId){
	$pepSeq = $trsptIdToPepSeq{$trsptId};
	$origTrsptId = "";
	# 当转录本对应的pep不存在：
	# 	将trsptId直接写入hash
	#	分析转录本是否为AS编辑后的编号，如果是那么用其原始转录本编号找到对应的基因编号
	if(not exists($hrefPepSeqToGeneTrsptIdList->{$pepSeq})){
		$hrefPepSeqToGeneTrsptIdList->{$pepSeq}->{"trsptIdList"} = $trsptId;
		# 找打原始trsptId，便于找对应的geneId
		if($trsptId=~/(.*)XXX(.*)/){
			$origTrsptId = $1;
		}else{
			$origTrsptId = $trsptId;
		}
		$hrefPepSeqToGeneTrsptIdList->{$pepSeq}->{"geneIdList"} = $trsptIdToGeneId{$origTrsptId};
	}else{
		$hrefPepSeqToGeneTrsptIdList->{$pepSeq}->{"trsptIdList"} .= "#" . $trsptId;
		if($trsptId=~/(.*)XXX(.*)/){
			$origTrsptId = $1;
		}else{
			$origTrsptId = $trsptId;
		}
		$hrefPepSeqToGeneTrsptIdList->{$pepSeq}->{"geneIdList"} .= "#" . $trsptIdToGeneId{$origTrsptId};
	}
}


# 以pepSeq为关键字输出trsptIdList和geneIdList
my ($geneIdListString, $trsptIdListString, @geneId);
my @pepSeq = keys(%pepSeqToGeneTrsptIdList);
open WW, ">$uniqPepFasta";
open CHECK, ">$checkPepIdListTsv";
foreach $pepSeq(@pepSeq){

	# 如果1条pep对应到n个geneId，那么重复产生n个相同pep。序列的id用geneId___trsptIdList的形式给出
	# 在混合的序列名称中，前面的geneIdList和后面的trsptIdlist完全一一对应
	%geneId = ();
	&getGeneIdListAndTrsptIdList($hrefPepSeqToGeneTrsptIdList->{$pepSeq}->{"geneIdList"}, $hrefPepSeqToGeneTrsptIdList->{$pepSeq}->{"trsptIdList"}, \%geneId);

	# 按照geneId重复输出pepSeq，确保每个geneId都单独有pep与之对应
	@geneId = ();
	@geneId = keys(%geneId);
	foreach $geneId(@geneId){
		print WW ">" . $geneId . "___" . $geneId{$geneId} . "\n";
		print WW $pepSeq . "\n";
	}
	# 对该pepSeq是否对应多个geneId做出判断，如果是，那么输出到check文件中便于人工检查
	if($#geneId>0){
		my $output = "";
		foreach $geneId(@geneId){
			if($output eq ""){
				$output = $geneId . "___" . $geneId{$geneId};
			}else{
				$output .= "===" . $geneId . "___" . $geneId{$geneId};
			}
		}
		print CHECK $output . "\n";
	}	

}
close WW,
close CHECK;

# 将uniqe geneId抽取出来构成一个字符串，用#分割
# 将每个geneId对应的trsptId也抽出来，构成字符串，用#分割
sub getGeneIdListAndTrsptIdList{
	my ($geneIdListString, $trsptIdListString, $geneIdHref) =@_;
	# input：g1#g1#g2#g1___g1.t1#g1.t2#g2.t1#g1.t3
	# output：g1#g2___g1.t1,g1.t2,g1.t3#g2.t1
	my (@geneId, @trsptId, $geneId, $trsptId, $i);
	@geneId = split(/#/, $geneIdListString);
	@trsptId = split(/#/, $trsptIdListString);
	for($i=0; $i<=$#geneId; $i++){
		if(not exists($geneIdHref->{$geneId[$i]})){
			$geneIdHref->{$geneId[$i]} = $trsptId[$i];
		}else{
			$geneIdHref->{$geneId[$i]} .= "#" . $trsptId[$i];
		}
	}
}

sub getGeneIdAndTrsptId{
	my ($attrString, $geneId, $trsptId) = @_;
	my (@attr, $attr);
	chomp($attrString);
	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$geneId = $1;
		}
		if($attr=~/transcript_id "(.*)"/){
			$$trsptId = $1;
		}
	}
}
