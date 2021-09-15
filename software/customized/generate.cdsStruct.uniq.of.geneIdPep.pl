#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--uniquePepFasta \\\n" .
                "--trsptCdsStrucGtf \\\n" .
                "--uniqPepCdsStrucGtf \n";
	exit;
}

my ($uniquePepFasta, $trsptCdsStrucGtf, $uniqPepCdsStrucGtf);

GetOptions(
        'uniquePepFasta=s'=>\$uniquePepFasta,
        'trsptCdsStrucGtf=s'=>\$trsptCdsStrucGtf,
        'uniqPepCdsStrucGtf=s'=>\$uniqPepCdsStrucGtf,
);

my (%trsptCdsStruc, $trsptCdsStrucHref, %pepToTrsptIdList);
$trsptCdsStrucHref=\%trsptCdsStruc;
# 将trspt的cds的外显子串结构以字符串形式读入到hash中
my ($cmd, $exonLineText, @exonLine, $exonLine, @field, $geneId, $trsptId);
$cmd = "cat $trsptCdsStrucGtf";
$exonLineText = `$cmd`;
@exonLine = split(/\n/, $exonLineText);
# print "finish read gtf into hash.\n";
foreach $exonLine(@exonLine){
#	print $exonLine;
#	<STDIN>;
	@field = ();
	@field = split(/\t/, $exonLine);
	if($field[8]=~/gene_id ".*?"; transcript_id "(.*?)";/){
		$trsptId = $1;
#		print $trsptId;
#		<STDIN>;
	}
	if($field[2] eq "transcript"){
		$trsptCdsStrucHref->{$trsptId} = $exonLine;
	}else{
		$trsptCdsStrucHref->{$trsptId} .= "\n" . $exonLine;
	}
}

# 将uniq pep fasta中序列名依次读取出来，只要输出一个trspt的外显子结构即可。
# >AT5G09410___SRX399568.18675.3XXXATHAA3SS0000000398#AT5G09410.2XXXATHASE0000000241#SRX1660582.15918.3XXXATHAA3SS0000000399#AT5G09410.1
# 注意：要将transcript_id改成如下：
# gene_id "AT5G09410"; transcript_id "SRX399568.18675.3XXXATHAA3SS0000000398#AT5G09410.2XXXATHASE0000000241#SRX1660582.15918.3XXXATHAA3SS0000000399#AT5G09410.1"
my ($geneIdAndTrsptIdListText, @geneIdAndTrsptIdList, $geneIdAndTrsptIdList);
my ($tsptIdList, @trsptId, $trsptId, $geneId, @exonLine, $exonLine, @field, @attr);
$cmd = "grep \">\" " . $uniquePepFasta . " |awk -F \'>\' \'{print \$2}\' ";
$geneIdAndTrsptIdListText = `$cmd`;
@geneIdAndTrsptIdList = split(/\n/, $geneIdAndTrsptIdListText);
open WW, ">$uniqPepCdsStrucGtf";
foreach $geneIdAndTrsptIdList(@geneIdAndTrsptIdList){
	($geneId, $tsptIdList) = split(/___/, $geneIdAndTrsptIdList);
	@trsptId = ();
	@trsptId = split(/#/, $tsptIdList);
	# 提取第1个trspt的exonStructure代表整个pep的结构
	$trsptId = $trsptId[0];
	@exonLine = ();
	@exonLine = split(/\n/, $trsptCdsStrucHref->{$trsptId});
	# 将每一行中transcript_id修改成 $geneIdAndTrsptIdList
	foreach $exonLine(@exonLine){
		@field = ();
		@field = split(/\t/, $exonLine);
		# field[8]: gene_id "AT3G53400"; transcript_id "AT3G53400.1"; exon_number "0";
		@attr = ();
		@attr = split(/; /, $field[8]);
		if($#attr==2){
			$field[8]= join("; ", $attr[0], "transcript_id \"" . $geneIdAndTrsptIdList . "\"", $attr[2]);
		}else{
			$field[8]= join("; ", $attr[0], "transcript_id \"" . $geneIdAndTrsptIdList . "\"");
		}
		print WW join("\t", @field) . "\n";
	}
	
}
close WW;
