#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw/sum/;
use Statistics::R;

if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--jobId 1234 \\\n" .
                "--sampleInfo sample.info.tsv \\\n" .
                "--taxonId1st 39947 \\\n" .
                "--taxonId2nd 3702 \\\n" .
                "--exptDir ./expts \\\n" .
		"--corPvalue 0.05 \\\n" .
		"--minReadCover 10 \\\n" .
		"--orthoAsGroupTsv orthoAsGroup.mysql.tsv\\\n" .
		"--outputDir ./outputDir \\\n" . 
		"--outputCoSplicingOrtho 123.39947.3702.cosplicing.orthoId.tsv\n";
	exit;
}

my ($sampleInfo, $taxonId1st, $taxonId2nd, $exptDir, $orthoAsGroupTsv);
my ($outputDir, $jobId, $corPvalue, $corCoef, $outputCoSplicingOrtho, $minReadCoverage);

GetOptions(
	'jobId=s'=>\$jobId,
        'sampleInfo=s'=>\$sampleInfo,
        'taxonId1st=s'=>\$taxonId1st,
        'taxonId2nd=s'=>\$taxonId2nd,
        'exptDir=s'=>\$exptDir,
	'corPvalue=s'=>\$corPvalue,
	'minReadCoverage=s'=>\$minReadCoverage,
	'orthoAsGroupTsv=s'=>\$orthoAsGroupTsv,
	'outputDir=s'=>\$outputDir,
	'outputCoSplicingOrtho=s'=>\$outputCoSplicingOrtho,
);

my ($line, @field, $taxonId, $exptId, $groupName);

system("mkdir -p $outputDir");

# +++ 1. 从expt到group和Taxon：%exptToGroup %exptToTaxon ++
# +++ 2. 从group到taxon和expt: %orthoToTaxonToAsIdList ++
# +++ 3. 登记所有exptId ++
my (%exptId, @exptId);
my (%exptToGroup, %exptToTaxon);
my (%groupToTaxonToExptIdList, $groupToTaxonToExptIdListHref);
$groupToTaxonToExptIdListHref = \%groupToTaxonToExptIdList;
open FF, "<$sampleInfo";
# taxonId exptId  	group
# 39947   SRX736223     root
<FF>;
while($line=<FF>){
	chomp($line);
	($taxonId, $exptId, $groupName) = split(/\t/, $line);
	$groupToTaxonToExptIdListHref->{$groupName}->{$taxonId} .= $exptId . ",";
	$exptToGroup{$exptId} = $groupName;
	$exptToTaxon{$exptId} = $taxonId;
	$exptId{$exptId} = 1;
}
close FF;

# 登记expt
@exptId = keys(%exptId);

# +++ 3. 所有group名称 +++
my @group = keys(%groupToTaxonToExptIdList);

# +++ 6. as在每个group中的psi列表：%asToGroupToPsi
my (%asToGroupToPsi, $asToGroupToPsiHref, %asId);
$asToGroupToPsiHref = \%asToGroupToPsi;
my ($asId, $IJC_SAMPLE_1, $SJC_SAMPLE_1, $IncFormLen, $SkipFormLen, $psi, $Experiment, $line);
foreach $exptId(@exptId){
	$groupName = $exptToGroup{$exptId};	
	open FF, "<$exptDir/$exptId";
	# ASID    		IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
	# OSJGA3SS0000000001    0               1               103             99              SRX864565
	<FF>;
	while($line = <FF>){
		chomp($line);
		($asId, $IJC_SAMPLE_1, $SJC_SAMPLE_1, $IncFormLen, $SkipFormLen) = split(/\t/, $line);

		# 检测read覆盖度是否达标
		next if($IJC_SAMPLE_1 + $SJC_SAMPLE_1 < $minReadCoverage);

		# print $line . "\n";
		$psi = ($IJC_SAMPLE_1/$IncFormLen)/(($IJC_SAMPLE_1/$IncFormLen) + ($SJC_SAMPLE_1/$SkipFormLen));
		$asToGroupToPsiHref->{$asId}->{$groupName} .= $psi . ",";
		$asId{$asId} = 1;
	}
	close FF;
}

# +++ 7. 计算as在每个group中psi均值，如果在某个group中没有检测到PSI，那么标志为-1 +++
# 如果在某个group没有检测到PSI，那么标志为-1
my (@psi, $withPsiInAllGroup, @asId);
@asId = keys(%asId);
print "总共搜集了: ";
print $#asId;
print "个asId\n";

foreach $asId(@asId){
	$withPsiInAllGroup = $#group+1;
	foreach $groupName(@group){
		if(exists($asToGroupToPsiHref->{$asId}->{$groupName})){
			@psi = ();
			@psi = split(/,/, $asToGroupToPsiHref->{$asId}->{$groupName});
			$asToGroupToPsiHref->{$asId}->{$groupName} = sum(@psi)/($#psi+1);
		}else{
			$withPsiInAllGroup--;
			$asToGroupToPsiHref->{$asId}->{$groupName} = -1;
		}
	}
	# 如果在任意一个group中没有检测到有效psi，那么放弃该asId
	if($withPsiInAllGroup < 3){
		delete($asId{$asId});
		delete($asToGroupToPsi{$asId});
	}
}

@asId = ();
@asId = keys(%asId);
print "具有3个以上group的AS数量：";
print $#asId;
print "\n";

# +++ 4. 从ortho到taxon和asId: %orthoToTaxonToAsIdList
# 读取orthoAs信息
my (%asToTaxon);
my (%orthoToTaxonToAsIdList, $orthoToTaxonToAsIdListHref);
$orthoToTaxonToAsIdListHref = \%orthoToTaxonToAsIdList;
my ($fieldString, $valueString);
my ($orthoId, $asId, $taxonId, $species, $abbr, $order, $speciesNum);
open FF, "<$orthoAsGroupTsv";
# asId, orthoAsGroupId, taxonId, species, speciesAbbr, orderName, speciesNum___TAESA3SS0000010572, orthA3SS00000017, 4565, Triticum aestivum, TAES, Poales, 1
<FF>;
while($line=<FF>){
	($fieldString, $valueString) = split(/___/, $line);
	($asId, $orthoId, $taxonId, $species, $abbr, $order, $speciesNum) = split(/, /, $valueString);
	# 只有在所有group中都能检测到有效PSI的asId才能被搜集到orthoToTaxonToAsIdList哈希中
	$orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId} .= $asId . "," if(exists($asId{$asId}));
	$asToTaxon{$asId} = $taxonId;
}
close FF;

# 保留至少有2个以上asId的orthoId
# (1) 2个以上as都在taxon1st中 
# (2) 2个以上as分布在taxon1st和taxon2nd中
# 换言之：
# (1) 在taxon1st中至少包含1个as
my (@orthoId, $orthoId);
my (%asId, $asId, @asId);
@orthoId = keys(%orthoToTaxonToAsIdList);
foreach $orthoId(@orthoId){	
	# taxon1st中不包含as的情况，将被剔除
	if(not(exists($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st}))){
		delete($orthoToTaxonToAsIdList{$orthoId});
		next;
	}
	# orthoId中只有1个as在taxon1st中，将被剔除
	if(not exists($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd})){
		@asId = ();
		@asId = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st});
		if($#asId == 0){
			delete($orthoToTaxonToAsIdList{$orthoId});
			next;
		}
	}

	# 登记asId到hash
	@asId = ();
	@asId = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st});
	foreach $asId(@asId){
		$asId{$asId} = 1;
	}
	@asId = ();
	@asId = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd});
	foreach $asId(@asId){
		$asId{$asId} = 1;
	}
}

### ----------------
# 此时向后：orthoToTaxonToAsIdList中每个orthoId中都至少包含1个TaxonId1st的as

# +++ 5. 所有orthoId和asId
@orthoId = ();
@orthoId = keys(%orthoToTaxonToAsIdList);
print join("\t", "最初orthogroup数量(至少2个as/orthogroup): ", $#orthoId) . "\n";
print "下面将剔除掉taxon1st中包含2个以上as，但是它们表达不一致的orthoId\n";

=beg
## ----过滤掉TaxonId1st内部as的表达模式不相关的orthoId -----
# 注意：会暂时保留只包含1个as的orthoId，将在以后计算它们是否和taxon2nd中as表达是否一致
my ($output1, $output2, $psiList);
my $innerTsv1 = "$outputDir/$jobId.$taxonId1st.ortho.inner.compare.group.mean.psi.list";
open WW, ">$innerTsv1";
print WW join("\t", "orthoId", "asId_1st", "meanPsiList1", "asId_2nd", "meanPsiList2") . "\n";
foreach $orthoId(@orthoId){
	# 
	@asId = ();
	@asId = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st});
	# 如果只有1个asId，那么暂时不关注
	next if($#asId == 0);
	# 列出该orthoId下所有两个asId之间的比较
	for(my $i=0; $i<=$#asId; $i++){
		for(my $j=$i+1; $j<=$#asId; $j++){
			# 列出asId[$i]在所有group中psi均值
			$output1 = "";
			foreach $groupName(@group){
				$output1 .= $asToGroupToPsiHref->{$asId[$i]}->{$groupName} . ",";
			}
			# 列出asId[$j]在所有group中psi均值
			$output2 = "";
			foreach $groupName(@group){
				$output2 .= $asToGroupToPsiHref->{$asId[$j]}->{$groupName} . ",";
			}
			# 输出整个比较行
			print WW join("\t", $orthoId, $asId[$i], substr($output1, 0, length($output1) - 1), 
					$asId[$j], substr($output2, 0, length($output2) - 1)
					) . "\n";

		}
	}
}
close WW;

# 利用R语言计算相关系数和pvalue值
# orthoId asId_1st      meanPsiList     asId_2nd      meanPsiList
my $R=Statistics::R->new(shared => 1);
my $cmd=<<EOF;
df <- read.table(file="$innerTsv1", sep="\\t", quote ="", header=T)
df\$corCoef <- apply(df, 1, function(x) cor(as.numeric(unlist(strsplit(x[3], ","))), as.numeric(unlist(strsplit(x[5], ","))), method = "pearson"))
df\$corPvalue <- apply(df, 1, function(x) cor.test(as.numeric(unlist(strsplit(x[3], ","))), as.numeric(unlist(strsplit(x[5], ","))), alternative="two.sided")\$p.value)
write.table(df, file='$innerTsv1', sep='\\t', quote=F, row.names = F)
EOF
#print $cmd;
$R->run($cmd);
$R->stop();

# 读取每个taxon1st内2个以上as的表达相关性结果，剔除掉那些表达模式不一致的orthoId
my (@field);
open FF, "<$innerTsv1";
# orthoId asId_1st meanPsiList1 asId_2nd meanPsiList2 corCoef corPvalue
<FF>;
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	delete($orthoToTaxonToAsIdList{$field[0]}) if($field[$#field] > $corPvalue or $field[$#field] eq "NA" or $field[$#field-1] <0);
}
close FF;
system("rm -rf $innerTsv1");

# 统计哈希中剩余的orthoId
@orthoId = ();
@orthoId = keys(%orthoToTaxonToAsIdList);
#print join("\t", "在taxonId 1st中剔除掉内部表达模式不一致的ortho，还剩下：", $#orthoId) . "\n";

# 目前哈希中剩余orthoId：
# (1) as只分布在taxon1st中，且它们表达模式一致
# (2) as分布在taxon1st和taxon2nd中
# 下面就是针对（2）进行相关性分析，剔除taxon1st和taxon2n中不一致的情况
### --- 计算交叉cross相关性 -----
my (@asId1st, @asId2nd);
my $crossTsv = "$outputDir/$jobId.$taxonId1st.$taxonId2nd.cross.compare.group.mean.psi.list";
open WW, ">$crossTsv";
print WW join("\t", "orthoId", "asId_taxonId1st", "meanPsiList1", "asId_taxonId2nd", "meanPsiList2") . "\n";
foreach $orthoId(@orthoId){

	# 如果taxon2nd中不包含as，那么不予检测
	next if(not(exists($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd})));

	# taxonId1st中的asId
	@asId1st = ();
	@asId1st = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st});

	# taxonId2nd中的asId
	@asId2nd = ();
	@asId2nd = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd});
	
	# 列出该orthoId下所有两个asId之间的比较
	for(my $i=0; $i<=$#asId1st; $i++){
		for(my $j=0; $j<=$#asId2nd; $j++){

			# 列出asId[$i]在所有group中psi均值
			$output1 = "";
			foreach $groupName(@group){
				$output1 .= $asToGroupToPsiHref->{$asId1st[$i]}->{$groupName} . ",";
			}

			# 列出asId[$j]在所有group中psi均值
			$output2 = "";
			foreach $groupName(@group){
				$output2 .= $asToGroupToPsiHref->{$asId2nd[$j]}->{$groupName} . ",";
			}

			# 输出整个比较行
			print WW join("\t", $orthoId, $asId1st[$i], substr($output1, 0, length($output1) - 1), 
					$asId2nd[$j], substr($output2, 0, length($output2) - 1)
					) . "\n";
		}
	}
}
close WW;

# 利用R语言计算相关系数和pvalue值
# orthoId asId_1st      meanPsiList     asId_2nd      meanPsiList
my $R=Statistics::R->new(shared => 1);
my $cmd=<<EOF;
df <- read.table(file="$crossTsv", sep="\\t", quote ="", header=T)
df\$corCoef <- apply(df, 1, function(x) cor(as.numeric(unlist(strsplit(x[3], ","))), as.numeric(unlist(strsplit(x[5], ","))), method = "pearson"))
df\$corPvalue <- apply(df, 1, function(x) cor.test(as.numeric(unlist(strsplit(x[3], ","))), as.numeric(unlist(strsplit(x[5], ","))), alternative="two.sided")\$p.value)
write.table(df, file='$crossTsv', sep='\\t', quote=F, row.names = F)
EOF
#print $cmd;
$R->run($cmd);
$R->stop();

# 筛选ortho，条件是：外部所有compare的相关系数和pvalue都满足条件
# 如果出现有相关性pvalue小于指定阈值，那么将该orthoId删除
open FF, "<$crossTsv";
# orthoId asId_1st meanPsiList1 asId_2nd meanPsiList2 corCoef corPvalue
<FF>;
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	delete($orthoToTaxonToAsIdList{$field[0]}) if($field[$#field] > $corPvalue or $field[$#field] eq "NA" or $field[$#field-1] <0);
}
close FF;
system("rm -rf $crossTsv");

# 剔除掉和taxon1st和taxon2nd中表达模式不一致的orthoId
@orthoId = ();
@orthoId = keys(%orthoToTaxonToAsIdList);
#print join("\t", "剔除as分布在2物种中，但表达模式不一致的ortho后，还剩下:", $#orthoId) . "\n";

my ($orthoIdNum_1st, $orthoIdNum_1_2);
foreach $orthoId(@orthoId){
	if(exists($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd}) and exists($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st})){
		$orthoIdNum_1_2++;
	}else{
		$orthoIdNum_1st++;
	}
}

#print "as只分布在taxonId1st中的orthoId数量：$orthoIdNum_1st\n";
#print "as分布在taxonId1st和taxonId2nd的orthoId数量：$orthoIdNum_1_2\n";



## 计算所有orthoId中的as表达模式的平均相关系数
my (@totalAsId);
my $finalTsv = "$outputDir/$jobId.$taxonId1st.$taxonId2nd.final.complete.cross.compare.group.mean.psi.list";
open WW, ">$finalTsv";
print WW join("\t", "orthoId", "asId_1st", "meanPsiList1", "asId_2nd", "meanPsiList2") . "\n";
foreach $orthoId(@orthoId){
	@totalAsId = ();
	# taxonId1st
	@asId = ();
	@asId = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st});
	push(@totalAsId, @asId);
	# taxonId2nd
	if(exists($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd})){
		@asId = ();
		@asId = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd});
		push(@totalAsId, @asId);
	}

	# 列出该orthoId下所有两个asId之间的比较
	for(my $i=0; $i<=$#totalAsId; $i++){
		for(my $j=$i+1; $j<=$#totalAsId; $j++){

			# 列出asId[$i]在所有group中psi均值
			$output1 = "";
			foreach $groupName(@group){
				$output1 .= $asToGroupToPsiHref->{$totalAsId[$i]}->{$groupName} . ",";
			}

			# 列出asId[$j]在所有group中psi均值
			$output2 = "";
			foreach $groupName(@group){
				$output2 .= $asToGroupToPsiHref->{$totalAsId[$j]}->{$groupName} . ",";
			}

			# 输出整个比较行
			print WW join("\t", $orthoId, $totalAsId[$i], substr($output1, 0, length($output1) - 1), 
					$totalAsId[$j], substr($output2, 0, length($output2) - 1)
					) . "\n";

		}
	}
}
close WW;

# 利用R语言计算相关系数和pvalue值
# orthoId asId_1st meanPsiList asId_2nd  meanPsiList corCoef corPvalue
my $R=Statistics::R->new(shared => 1);
my $cmd=<<EOF;
df <- read.table(file="$finalTsv", sep="\\t", quote ="", header=T)
df\$corCoef <- apply(df, 1, function(x) cor(as.numeric(unlist(strsplit(x[3], ","))), as.numeric(unlist(strsplit(x[5], ","))), method = "pearson"))
df\$corPvalue <- apply(df, 1, function(x) cor.test(as.numeric(unlist(strsplit(x[3], ","))), as.numeric(unlist(strsplit(x[5], ","))), alternative="two.sided")\$p.value)
write.table(df, file='$finalTsv', sep='\\t', quote=F, row.names = F)
EOF
#print $cmd;
$R->run($cmd);
$R->stop();

# 计算cross、taxon1st和taxon2nd中as表达模式平均相关系数
my (%orthoToTaxonCoef, $orthoToTaxonCoefHref, $asId1, $asId2, $tmpList1, $psiList2, $coef, $pvalue);
$orthoToTaxonCoefHref = \%orthoToTaxonCoef;
open FF, "<$finalTsv";
# orthoId asId_1st meanPsiList1 asId_2nd meanPsiList2 corCoef corPvalue
<FF>;
while($line=<FF>){
	chomp($line);
	($orthoId, $asId1, $tmpList1, $asId2, $psiList2, $coef, $pvalue) = split(/\t/, $line);
	# 属于taxonId1st内的比较
	if($asToTaxon{$asId1} eq $asToTaxon{$asId2} and $asToTaxon{$asId2} eq $taxonId1st){
		$orthoToTaxonCoefHref->{$orthoId}->{$taxonId1st} .= $coef . ",";
	# 属于taxonId2n内的比较
	}elsif($asToTaxon{$asId1} eq $asToTaxon{$asId2} and $asToTaxon{$asId2} eq $taxonId2nd){
		$orthoToTaxonCoefHref->{$orthoId}->{$taxonId2nd} .= $coef . ",";
	# 属于taxonId1st和taxonId2nd间的比较
	}elsif($asToTaxon{$asId1} ne $asToTaxon{$asId2}){
		$orthoToTaxonCoefHref->{$orthoId}->{"cross"} .= $coef . ",";
	}
}
close FF;
system("rm -rf $finalTsv");

# 重新输出每个orthoId，并计算平均值
my ($outputLine, @coef, $asNum, $avgCoef, $coefList);
#my $finalRlt = $outputDir . "/$jobId.$taxonId1st.$taxonId2nd.final.ortho.coexpress.result.tsv";
@orthoId = ();
@orthoId = keys(%orthoToTaxonCoef);
open WW, ">$outputDir/$outputCoSplicingOrtho";
print WW join("\t", "orthoId", $taxonId1st . "_asList", $taxonId1st . "_asNum", $taxonId1st . "_coef", $taxonId2nd . "_asList", $taxonId2nd . "_asNum", $taxonId2nd . "_coef", "cross_coef") . "\n";
foreach $orthoId(@orthoId){
	$outputLine = $orthoId;
	# taxonId1st
	$outputLine .= "\t" . substr($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st}, 0, length($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st})-1);
	@asId = ();
	@asId = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId1st});
	$asNum = $#asId + 1;
	$outputLine .= "\t" . $asNum;
	#  平均coef
	if(exists($orthoToTaxonCoef{$orthoId}->{$taxonId1st})){
		$coefList = $orthoToTaxonCoef{$orthoId}->{$taxonId1st};
		@coef = ();
		@coef = split(/,/, $orthoToTaxonCoef{$orthoId}->{$taxonId1st});
		$avgCoef = sum(@coef)/($#coef+1);
		$outputLine .= "\t" . $avgCoef;
	}else{
		$outputLine .= "\t-";
	}

	# taxonId2nd
	if(not(exists($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd}))){
		$outputLine .= "\t-\t-\t-";
	}else{
		$outputLine .= "\t" . substr($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd}, 0, length($orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd}) - 1);		
		@asId = ();
		@asId = split(/,/, $orthoToTaxonToAsIdListHref->{$orthoId}->{$taxonId2nd});
		$asNum = $#asId + 1;
		$outputLine .= "\t" . $asNum;
		if(exists($orthoToTaxonCoef{$orthoId}->{$taxonId2nd})){
			$coefList .= $orthoToTaxonCoef{$orthoId}->{$taxonId2nd};
			@coef = ();
			@coef = split(/,/, $orthoToTaxonCoef{$orthoId}->{$taxonId2nd});
			$avgCoef = sum(@coef)/($#coef+1);
			$outputLine .= "\t" . $avgCoef;
		}else{
			 $outputLine .= "\t-";
		}
	}

	# 最后计算交叉平均coef
	if(exists($orthoToTaxonCoef{$orthoId}->{"cross"})){
		$coefList .= $orthoToTaxonCoef{$orthoId}->{"cross"};
	}
	@coef = ();
	@coef = split(/,/, $coefList);
	$avgCoef = sum(@coef)/($#coef+1);
	$outputLine .= "\t" . $avgCoef . "\n";

	# 输出最后结果
	print WW $outputLine;
}
close WW;

=cut;
