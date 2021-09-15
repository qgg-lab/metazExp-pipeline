#!/usr/bin/perl
use strict;
my ($inputGenomeFile, $outputGenomeFile, $inpurtGtf, $outputGtf, $repBed) = @ARGV;
if($#ARGV<4){
	print "$0 orig.genome.fa replaced.genome.fa orig.gtf replaced.gtf reption.bed\n";
	exit;
}

# read genome sequence
my (%genomeSeq, $genomeSeqHref);
$genomeSeqHref=\%genomeSeq;
my ($line, @chr, $chr);

# 将基因组读入hash
open FF, "<$inputGenomeFile";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@chr = split(/ /, $1);
		$chr = $chr[0];
	}else{
		$genomeSeqHref->{$chr}.=$line;
	}
}
close FF;

#my $local_datestring = localtime();
#print "finish loading genome sequence: " . $local_datestring . "\n";

# 1、按照坐标从小到大对每个replace进行编号，然后将其读入到%replace
# 2、对replace按照染色体进行分组，然后按照坐标次序将replace编号归类到各个染色体中，存储在%chrToRepIdList
my (%replace, $replaceHref, $replaceNum, @field, %chrToRepIdList, $chrToRepIdListHref);
my (@repId, $repId, $i, $j);
$replaceHref=\%replace;
$chrToRepIdListHref = \%chrToRepIdList;
open FF, "<$repBed";
while($line=<FF>){
	# $replaceNum 为replace编号，格式为sprintf("%010d", $replaceNum)
	$replaceNum++;
	$repId=sprintf("%010d", $replaceNum);
	chomp($line);
	@field = split(/\t/, $line);
	$chr = $field[0];

	# 将replace操作信息保存到%replace，即$replaceHref引用的hash中
	$replaceHref->{$chr}->{$repId}->{"start"} = $field[1];
	$replaceHref->{$chr}->{$repId}->{"stop"} = $field[2];
	$replaceHref->{$chr}->{$repId}->{"origBases"} = $field[3];
	$replaceHref->{$chr}->{$repId}->{"repBases"} = $field[4];

	# 计算该replace操作执行后对其后续replace坐标的影响。
	#     如果替换碱基个数大于被替换碱基个数，那么会导致该replace之后的replace坐标都要增加；相反后replace坐标都要减小
	#     注意：如果是直接deletion型（即*），那么对后续的影响长度就是 负 原碱基串长度（即length($field[3])）
	if($field[4] eq "*"){
		$replaceHref->{$chr}->{$repId}->{"affectLen"} = 0 - length($field[3]);
	}else{
		$replaceHref->{$chr}->{$repId}->{"affectLen"} = length($field[4]) - length($field[3]);
	}

	# 将replace编号按照染色体归类
	if($chrToRepIdListHref->{$chr} eq ""){
		$chrToRepIdListHref->{$chr} = $repId;
	}else{
		$chrToRepIdListHref->{$chr} .= "," .$repId
	}

}
close FF;

#$local_datestring = localtime();
#print "finish loading replace: " . $local_datestring . "\n";


# 按照坐标从小到大对gtf中每一行进行编号，然后将每一行按照chr分组，保存到%chrToGtf
# 将gtf中行编号按照坐标从小到大登记到每个chr中，即保存到%chrToGtfLineIdList中
my (%chrToGtf, $chrToGtfHref, $gtfLineNum, %chrToGtfLineIdList, $chrToGtfLineIdListHref);
my ($gtfLineIdList, @gtfLineId, $gtfLineId);
$chrToGtfHref=\%chrToGtf;
$chrToGtfLineIdListHref=\%chrToGtfLineIdList;
open FF, "<$inpurtGtf";
while($line=<FF>){
	@field = split(/\t/, $line);
	next if($#field!=8);
	
	# 对gtf行统一编号，格式为sprintf("%010d", $gtfLineNum)
	$gtfLineNum++;
	$gtfLineId = sprintf("%010d", $gtfLineNum);
	$chr = $field[0];
	# gtf行坐标左侧部分登记进hash
	$chrToGtfHref->{$chr}->{$gtfLineId}->{"left"} = join("\t", $field[0], $field[1], $field[2]);
	# gtf行的两个坐标登记到hash
	$chrToGtfHref->{$chr}->{$gtfLineId}->{"start"} = $field[3];
	$chrToGtfHref->{$chr}->{$gtfLineId}->{"stop"} = $field[4];
	# gtf行坐标右侧部分登记进hash
	$chrToGtfHref->{$chr}->{$gtfLineId}->{"right"} = join("\t", $field[5], $field[6], $field[7], $field[8]);

	# 将gtf行编号按照染色体分组登记到hash中
	if($chrToGtfLineIdListHref->{$chr} eq ""){
		$chrToGtfLineIdListHref->{$chr} = $gtfLineId;
	}else{
		$chrToGtfLineIdListHref->{$chr} .= "," . $gtfLineId;
	}
}
close FF;


#$local_datestring = localtime();
#print "finish loading gtf: " . $local_datestring . "\n";
#print "Begin update genome sequence and cooradinates of gtf\n";

#
# 提取每个热染色体中的replaceId，按照坐标次序执行每个replace操作：
#     1、对genome产生序列替换的影响
#     2、对当前replace之后replace的坐标产生变大或者变小的影响
#     3、对当强replace之后gtf行的坐标产生变大或者变小的影响

my ($affectLen, $repIdList, $beforeSeq, $afterSeq);
my ($repStartBased1, $repStopBased1, $overlapBeginGtfLineJ, $overlapEndGtfLineJ);
my ($cutBeginBased1, $cutEndBased1, $insertPosBased1, $insertLen, $k);
@chr = keys(%replace);

foreach $chr(@chr){
	
	# 获得登记在该chr名下的所有replaceId
	$repIdList = $chrToRepIdListHref->{$chr};
	@repId = ();
	@repId = split(/,/, $repIdList);

	# 获得登记在该chr名下的所有的gtfLineId
	$gtfLineIdList = $chrToGtfLineIdListHref->{$chr};
	@gtfLineId = ();
	@gtfLineId = split(/,/, $gtfLineIdList);

	# 将当前chr所有repId输出查看
#	print "*******chr$chr" . "上, 所有的replace:\n";
#	for($j=0; $j<=$#repId; $j++){
#		print join("\t", $chr, $replaceHref->{$chr}->{$repId[$j]}->{"start"}, $replaceHref->{$chr}->{$repId[$j]}->{"stop"}, $replaceHref->{$chr}->{$repId[$j]}->{"origBases"}, $replaceHref->{$chr}->{$repId[$j]}->{"repBases"}) . "\n";
#	}
#	<STDIN>;
	# 将当前chr所有gtf输出查看
#	print "\n*******chr$chr" . "上,所有的gtfLine:\n";
#	for($j=0; $j<=$#gtfLineId; $j++){
#		print join("\t", $chr, $chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"start"}, $chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"stop"}) . "\n";
#	}
#	<STDIN>;

	# 依次执行每个replace操作
	for($i=0; $i<=$#repId; $i++){
	
		$repId = $repId[$i];

		# 将当前所有repId输出查看
#		print "\n===Begin execute:\"" . join("\t", $chr, $replaceHref->{$chr}->{$repId}->{"start"}, $replaceHref->{$chr}->{$repId}->{"stop"}, $replaceHref->{$chr}->{$repId}->{"origBases"}, $replaceHref->{$chr}->{$repId}->{"repBases"}) . "\",后续的replace如下：\n";
#		for($j=$i+1; $j<=$#repId; $j++){
#			print join("\t", $chr, $replaceHref->{$chr}->{$repId[$j]}->{"start"}, $replaceHref->{$chr}->{$repId[$j]}->{"stop"}, $replaceHref->{$chr}->{$repId[$j]}->{"origBases"}, $replaceHref->{$chr}->{$repId[$j]}->{"repBases"}) . "\n";
#		}
		#<STDIN>;


################# 1、执行replace操作后，更新基因组序列
		#     注意：要区分是否为直接deletion类型
		# 染色体上replace操作区域前后的序列
		$beforeSeq = substr($genomeSeqHref->{$chr}, 0, $replaceHref->{$chr}->{$repId}->{"start"});
		$afterSeq = substr($genomeSeqHref->{$chr}, $replaceHref->{$chr}->{$repId}->{"start"} + length($replaceHref->{$chr}->{$repId}->{"origBases"}));
		if($replaceHref->{$chr}->{$repId}->{"repBases"} ne "*"){
			$genomeSeqHref->{$chr} = $beforeSeq . $replaceHref->{$chr}->{$repId}->{"repBases"} . $afterSeq;
		}else{
			$genomeSeqHref->{$chr} = $beforeSeq . $afterSeq;
		}

################# 如果是长度一致的替换，那么不需要处理repId之后的坐标，以及gtf上的坐标
		if($replaceHref->{$chr}->{$repId}->{"affectLen"} == 0){
#			print "===After execute:\"" . join("\t", $chr, $replaceHref->{$chr}->{$repId}->{"start"}, $replaceHref->{$chr}->{$repId}->{"stop"}, $replaceHref->{$chr}->{$repId}->{"origBases"}, $replaceHref->{$chr}->{$repId}->{"repBases"}) . "\", 对后续replace无影响.\n";
#			print "===After execute:\"" . join("\t", $chr, $replaceHref->{$chr}->{$repId}->{"start"}, $replaceHref->{$chr}->{$repId}->{"stop"}, $replaceHref->{$chr}->{$repId}->{"origBases"}, $replaceHref->{$chr}->{$repId}->{"repBases"}) . "\", 对gtfLine无影响.\n";
			next
		}


################ 2、执行replace操作后，对当前replace操作后续所有操作的坐标都会产生影响
		#    修改当前replace之后的replace的坐标
		$affectLen = $replaceHref->{$chr}->{$repId}->{"affectLen"};
		for($j=$i+1; $j<=$#repId; $j++){
			$replaceHref->{$chr}->{$repId[$j]}->{"start"} += $affectLen;
			$replaceHref->{$chr}->{$repId[$j]}->{"stop"} += $affectLen;
		}
	
		# 将当前所有repId输出查看
#		print "\n===After execute \"" . $i . "\", 后续replace如下：\n";
#		for($j=$i+1; $j<=$#repId; $j++){
#			print join("\t", $chr, $replaceHref->{$chr}->{$repId[$j]}->{"start"}, $replaceHref->{$chr}->{$repId[$j]}->{"stop"}, $replaceHref->{$chr}->{$repId[$j]}->{"origBases"}, $replaceHref->{$chr}->{$repId[$j]}->{"repBases"}) . "\n";
#		}
#		<STDIN>;


################ 3、执行replace操作后，对当前replace操作位置后的gtfLine的坐标都会产生影响。
################    受影响的gtfLine被分为两部分：和replace重叠overlap部分; 和replace不重叠，落在replace后的部分

	######## (1) 找出和replace重叠的gtfLine，
		#     编号从$overlapBeginGtfLineJ到$overlapEndGtfLineJ
		$repStartBased1 = $replaceHref->{$chr}->{$repId}->{"start"} + 1;
		$repStopBased1 = $replaceHref->{$chr}->{$repId}->{"stop"};
		$overlapBeginGtfLineJ = 0;
		$overlapEndGtfLineJ = 0;

		# replace可能在gtfLine之间，可能在gtfLine内，可能覆盖部分gtfLine，可能完覆盖多个gtfLine
		# 覆盖区域为：
		#   replace和gtfLine有重叠
		for($j=0; $j<=$#gtfLineId; $j++){
			if($chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"stop"} >= $repStartBased1 and $chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"left"} ne "deleted"){
				$overlapBeginGtfLineJ = $j;
				last;
			}
		}
		#  找到和replace无重叠的第1个gtfLine
		for($j=$overlapBeginGtfLineJ; $j<=$#gtfLineId; $j++){
			if($chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"start"} > $repStopBased1 and $chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"left"} ne "deleted"){
				$overlapEndGtfLineJ = $j;
				last;
			}
		}
		
		# 将当前chr所有gtf输出查看
#		print "\n===Before execute:\"" . join("\t", $chr, $replaceHref->{$chr}->{$repId}->{"start"}, $replaceHref->{$chr}->{$repId}->{"stop"}, $replaceHref->{$chr}->{$repId}->{"origBases"}, $replaceHref->{$chr}->{$repId}->{"repBases"}) . "\"\n,将对如下gtfLine产生影响:\n";
#		for($j=$overlapBeginGtfLineJ; $j<=$#gtfLineId; $j++){
#			print join("\t", $chr, $chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"start"}, $chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"stop"}) . "\n";
#		}
		#<STDIN>;
#		print "affectLen=" . $affectLen . "\n";
	######## (2) 更新和replace重叠部分的gtfLine ######
		# 替换后的效果为剪掉一段序列
		if($affectLen < 0){

			$cutBeginBased1 = $replaceHref->{$chr}->{$repId}->{"start"} + length($replaceHref->{$chr}->{$repId}->{"repBases"}) + 1;
			$cutEndBased1 = $replaceHref->{$chr}->{$repId}->{"stop"};			
#			print "cutBeginBased1:$cutBeginBased1,cutEndBased1:$cutEndBased1\n";

			# 扫描重叠部分的gtfLine，判断replace和gtfLine关系
			for($k=$overlapBeginGtfLineJ; $k<$overlapEndGtfLineJ; $k++){

				if($cutBeginBased1 <= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} and $cutEndBased1 >= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"}){
					# gtfLine被完全剪切掉，那么直接将该gtfLine内容置空
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"left"} = "deleted";
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"right"} = "deleted";

				}elsif($cutEndBased1 >= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} and $cutEndBased1 <= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"}){ 
					# 剪掉gtfLine的左侧,由于左侧序列被剪掉，导致gtfLine剩余部分的坐标发生改变
					# a)将gtfLine的start坐标改到剪切点之后
					# b)gtfLine的start和stop的坐标都需要移位affectLen
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} = $cutEndBased1 + 1;
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} += $affectLen;
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"} += $affectLen;
				}elsif($cutBeginBased1 >= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} and $cutBeginBased1 <= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"}){
					# 剪掉的是gtfLine的右侧
					# 只需要将gtfLine的stop改小即可，不会对gtfLine的start和stop产生影响
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"} = $cutBeginBased1 - 1;
					
				}elsif($cutBeginBased1 >= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} and $cutEndBased1 <= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"}){
					# 剪掉的是gtfLine的中间部分，即cutBeginBased1落在gtfLine内部。
					# a)对gtfLine的stop执行移位即可，start坐标不受影响
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"} += $affectLen;

				}elsif($cutBeginBased1 < $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} and $cutEndBased1 < $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"}){
					# replace在gtfLine左侧，并没有落在gtfLine上，那么将start和stop执行移位操作
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} += $affectLen;
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"} += $affectLen;
				}
			}
		}elsif($affectLen > 0){	# 替换后的效果为插入一段序列
			$insertPosBased1 = $replaceHref->{$chr}->{$repId}->{"start"} + length($replaceHref->{$chr}->{$repId}->{"orgiBases"}) + 1;
#			print "insertPosBased1=" . $insertPosBased1 . "\n";
			for($k=$overlapBeginGtfLineJ; $k<$overlapEndGtfLineJ; $k++){
				if($insertPosBased1 <= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"} and $insertPosBased1 >= $chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"}){
					# 只需要将gtfLine的stop发生移位即可，对start无影响
					$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"} += $affectLen;
				}
			}
		}

	######## (3)更新replace操作位置$overlapEndGtfLineJ后gtfLine中的坐标
#		print "affect=$affectLen\n";
		for(my $k=$overlapEndGtfLineJ; $k<=$#gtfLineId; $k++){
			$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"start"} += $affectLen;
			$chrToGtfHref->{$chr}->{$gtfLineId[$k]}->{"stop"} += $affectLen;			
		}

		# 将当前chr所有gtf输出查看
#		print "\n===After execute:\"" . join("\t", $chr, $replaceHref->{$chr}->{$repId}->{"start"}, $replaceHref->{$chr}->{$repId}->{"stop"}, $replaceHref->{$chr}->{$repId}->{"origBases"}, $replaceHref->{$chr}->{$repId}->{"repBases"}) . "\"\n,将上述的gtfLine改为如下的gtfLine:\n";

#		for($j=$overlapBeginGtfLineJ; $j<=$#gtfLineId; $j++){
#			print join("\t", $chr, $chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"start"}, $chrToGtfHref->{$chr}->{$gtfLineId[$j]}->{"stop"}) . "\n";
#		}
#		<STDIN>;

	}
}

#$local_datestring = localtime();
#print "finish updating genome sequence and gtf: " . $local_datestring . "\n";


# 输出更新后的genome sequence
open WW, ">$outputGenomeFile";
@chr = ();
@chr = keys(%genomeSeq);
foreach $chr(@chr){
	print WW ">" .  $chr . "\n";
	print WW $genomeSeq{$chr} . "\n";
}
close WW;

#$local_datestring = localtime();
#print "finish writing new genome sequence into file: " . $local_datestring . "\n";


# 输出更新后的gtf
open WW, ">$outputGtf";
@chr = keys(%chrToGtfLineIdList);
print "chr num:" . $#chr . "\n";
foreach $chr(@chr){
	$gtfLineIdList = $chrToGtfLineIdListHref->{$chr};
	@gtfLineId = ();
	@gtfLineId = split(/,/, $gtfLineIdList);
#	print "$chr has " . $#gtfLineId . " gtfLine.\n";
	foreach $gtfLineId(@gtfLineId){
		if($chrToGtfHref->{$chr}->{$gtfLineId}->{"left"} ne "deleted"){
			print WW join("\t", $chrToGtfHref->{$chr}->{$gtfLineId}->{"left"}, $chrToGtfHref->{$chr}->{$gtfLineId}->{"start"}, $chrToGtfHref->{$chr}->{$gtfLineId}->{"stop"}, $chrToGtfHref->{$chr}->{$gtfLineId}->{"right"});
		}
	}
}
close WW;

#$local_datestring = localtime();
#print "finish write gtf into new gtf file: " . $local_datestring . "\n";

