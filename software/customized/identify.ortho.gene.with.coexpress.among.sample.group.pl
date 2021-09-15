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
		"--orthoGeneGroupTsv orthoGeneGroup.mysql.tsv\\\n" .
		"--outputDir ./outputDir \\\n" . 
		"--outputCoexpressOrtho 123.39947.3702.coexpress.orthoId.tsv\n";
	exit;
}

my ($sampleInfo, $taxonId1st, $taxonId2nd, $exptDir, $orthoGeneGroupTsv);
my ($outputDir, $jobId, $corPvalue, $corCoef, $outputCoexpressOrtho);

GetOptions(
	'jobId=s'=>\$jobId,
        'sampleInfo=s'=>\$sampleInfo,
        'taxonId1st=s'=>\$taxonId1st,
        'taxonId2nd=s'=>\$taxonId2nd,
        'exptDir=s'=>\$exptDir,
	'corPvalue=s'=>\$corPvalue,
	'orthoGeneGroupTsv=s'=>\$orthoGeneGroupTsv,
	'outputDir=s'=>\$outputDir,
	'outputCoexpressOrtho=s'=>\$outputCoexpressOrtho,
);

my ($line, @field, $taxonId, $exptId, $groupName);

system("mkdir -p $outputDir");

# +++ 1. 从expt到group和Taxon：%exptToGroup %exptToTaxon ++
# +++ 2. 从group到taxon和expt: %orthoToTaxonToGeneIdList ++
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

# +++ 4. 从ortho到taxon和geneId: %orthoToTaxonToGeneIdList
# 读取orthoGene信息
my (%geneToTaxon);
my (%orthoToTaxonToGeneIdList, $orthoToTaxonToGeneIdListHref);
$orthoToTaxonToGeneIdListHref = \%orthoToTaxonToGeneIdList;
my ($fieldString, $valueString);
my ($orthoId, $geneId, $taxonId, $species, $abbr, $order);
open FF, "<$orthoGeneGroupTsv";
# geneId, orthoGeneGroupId, taxonId, species, speciesAbbr, orderName___LOC110412305, OG0000000, 108875, Herrania umbratica, HUMB, Malvales
<FF>;
while($line=<FF>){
	($fieldString, $valueString) = split(/___/, $line);
	($geneId, $orthoId, $taxonId, $species, $abbr, $order) = split(/, /, $valueString);
	$orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId} .= $geneId . ",";
	$geneToTaxon{$geneId} = $taxonId;
}
close FF;

# 保留至少有2个以上geneId的orthoId
# (1) 2个以上gene都在taxon1st中 
# (2) 2个以上gene分布在taxon1st和taxon2nd中
# 换言之：
# (1) 在taxon1st中至少包含1个gene
my (@orthoId, $orthoId);
my (%geneId, $geneId, @geneId);
@orthoId = keys(%orthoToTaxonToGeneIdList);
foreach $orthoId(@orthoId){
	
	# taxon1st中不包含gene的情况，将被剔除
	if(not(exists($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st}))){
		delete($orthoToTaxonToGeneIdList{$orthoId});
		next;
	}

	# orthoId中只有1个gene在taxon1st中，将被剔除
	if(not exists($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd})){
		@geneId = ();
		@geneId = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st});
		if($#geneId == 0){
			delete($orthoToTaxonToGeneIdList{$orthoId});
			next;
		}
	}

	# 登记geneId到hash
	@geneId = ();
	@geneId = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st});
	foreach $geneId(@geneId){
		$geneId{$geneId} = 1;
	}
	@geneId = ();
	@geneId = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd});
	foreach $geneId(@geneId){
		$geneId{$geneId} = 1;
	}

}

### ----------------
# 此时向后：orthoToTaxonToGeneIdList中每个orthoId中都至少包含1个TaxonId1st的gene

# +++ 5. 所有orthoId和geneId
@orthoId = ();
@orthoId = keys(%orthoToTaxonToGeneIdList);
#print join("\t", "最初orthogroup数量(至少2个gene/orthogroup): ", $#orthoId) . "\n";
#print "下面将剔除掉taxon1st中包含2个以上gene，但是它们表达不一致的orthoId\n";

# +++ 6. gene在每个group中的tpm列表：%geneToGroupToTpm
my (%geneToGroupToTpm, $geneToGroupToTpmHref);
$geneToGroupToTpmHref = \%geneToGroupToTpm;
my ($geneId, $readCount, $coverage, $fpkm, $tpm, $line);
foreach $exptId(@exptId){
	$groupName = $exptToGroup{$exptId};	
	open FF, "<$exptDir/$exptId";
	# geneId        readCount  Coverage      FPKM        TPM
	# Os12g0542000  1784       136.099396    112.314980  161.494064	
	<FF>;
	while($line = <FF>){
		chomp($line);
		($geneId, $readCount, $coverage, $fpkm, $tpm) = split(/\t/, $line);
		$geneToGroupToTpmHref->{$geneId}->{$groupName} .= $tpm . ",";
	}
	close FF;
}


# +++ 7. 计算gene在每个group中tpm均值 +++
my (@tpm);
@geneId = keys(%geneId);
foreach $geneId(@geneId){
	foreach $groupName(@group){
		if(exists($geneToGroupToTpmHref->{$geneId}->{$groupName})){
			@tpm = ();
			@tpm = split(/,/, $geneToGroupToTpmHref->{$geneId}->{$groupName});
			$geneToGroupToTpmHref->{$geneId}->{$groupName} = sum(@tpm)/($#tpm+1);
		}
	}
}


## ----过滤掉TaxonId1st内部gene的表达模式不相关的orthoId -----
# 注意：会暂时保留只包含1个gene的orthoId，将在以后计算它们是否和taxon2nd中gene表达是否一致
my ($output1, $output2, $tpmList);
my $innerTsv1 = "$outputDir/$jobId.$taxonId1st.ortho.inner.compare.group.mean.tpm.list";
open WW, ">$innerTsv1";
print WW join("\t", "orthoId", "geneId_1st", "meanTpmList1", "geneId_2nd", "meanTpmList2") . "\n";
foreach $orthoId(@orthoId){
	# 
	@geneId = ();
	@geneId = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st});
	# 如果只有1个geneId，那么暂时不关注
	next if($#geneId == 0);
	# 列出该orthoId下所有两个geneId之间的比较
	for(my $i=0; $i<=$#geneId; $i++){
		for(my $j=$i+1; $j<=$#geneId; $j++){
			# 列出geneId[$i]在所有group中tpm均值
			$output1 = "";
			foreach $groupName(@group){
				$output1 .= $geneToGroupToTpmHref->{$geneId[$i]}->{$groupName} . ",";
			}
			# 列出geneId[$j]在所有group中tpm均值
			$output2 = "";
			foreach $groupName(@group){
				$output2 .= $geneToGroupToTpmHref->{$geneId[$j]}->{$groupName} . ",";
			}
			# 输出整个比较行
			print WW join("\t", $orthoId, $geneId[$i], substr($output1, 0, length($output1) - 1), 
					$geneId[$j], substr($output2, 0, length($output2) - 1)
					) . "\n";

		}
	}
}
close WW;

# 利用R语言计算相关系数和pvalue值
# orthoId geneId_1st      meanTpmList     geneId_2nd      meanTpmList
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

# 读取每个taxon1st内2个以上gene的表达相关性结果，剔除掉那些表达模式不一致的orthoId
my (@field);
open FF, "<$innerTsv1";
# orthoId geneId_1st meanTpmList1 geneId_2nd meanTpmList2 corCoef corPvalue
<FF>;
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	delete($orthoToTaxonToGeneIdList{$field[0]}) if($field[$#field] > $corPvalue or $field[$#field] eq "NA" or $field[$#field-1] <0);
}
close FF;
system("rm -rf $innerTsv1");

# 统计哈希中剩余的orthoId
@orthoId = ();
@orthoId = keys(%orthoToTaxonToGeneIdList);
#print join("\t", "在taxonId 1st中剔除掉内部表达模式不一致的ortho，还剩下：", $#orthoId) . "\n";

# 目前哈希中剩余orthoId：
# (1) gene只分布在taxon1st中，且它们表达模式一致
# (2) gene分布在taxon1st和taxon2nd中
# 下面就是针对（2）进行相关性分析，剔除taxon1st和taxon2n中不一致的情况
### --- 计算交叉cross相关性 -----
my (@geneId1st, @geneId2nd);
my $crossTsv = "$outputDir/$jobId.$taxonId1st.$taxonId2nd.cross.compare.group.mean.tpm.list";
open WW, ">$crossTsv";
print WW join("\t", "orthoId", "geneId_taxonId1st", "meanTpmList1", "geneId_taxonId2nd", "meanTpmList2") . "\n";
foreach $orthoId(@orthoId){

	# 如果taxon2nd中不包含gene，那么不予检测
	next if(not(exists($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd})));

	# taxonId1st中的geneId
	@geneId1st = ();
	@geneId1st = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st});

	# taxonId2nd中的geneId
	@geneId2nd = ();
	@geneId2nd = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd});
	
	# 列出该orthoId下所有两个geneId之间的比较
	for(my $i=0; $i<=$#geneId1st; $i++){
		for(my $j=0; $j<=$#geneId2nd; $j++){

			# 列出geneId[$i]在所有group中tpm均值
			$output1 = "";
			foreach $groupName(@group){
				$output1 .= $geneToGroupToTpmHref->{$geneId1st[$i]}->{$groupName} . ",";
			}

			# 列出geneId[$j]在所有group中tpm均值
			$output2 = "";
			foreach $groupName(@group){
				$output2 .= $geneToGroupToTpmHref->{$geneId2nd[$j]}->{$groupName} . ",";
			}

			# 输出整个比较行
			print WW join("\t", $orthoId, $geneId1st[$i], substr($output1, 0, length($output1) - 1), 
					$geneId2nd[$j], substr($output2, 0, length($output2) - 1)
					) . "\n";
		}
	}
}
close WW;

# 利用R语言计算相关系数和pvalue值
# orthoId geneId_1st      meanTpmList     geneId_2nd      meanTpmList
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
# orthoId geneId_1st meanTpmList1 geneId_2nd meanTpmList2 corCoef corPvalue
<FF>;
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	delete($orthoToTaxonToGeneIdList{$field[0]}) if($field[$#field] > $corPvalue or $field[$#field] eq "NA" or $field[$#field-1] <0);
}
close FF;
system("rm -rf $crossTsv");

# 剔除掉和taxon1st和taxon2nd中表达模式不一致的orthoId
@orthoId = ();
@orthoId = keys(%orthoToTaxonToGeneIdList);
#print join("\t", "剔除gene分布在2物种中，但表达模式不一致的ortho后，还剩下:", $#orthoId) . "\n";

my ($orthoIdNum_1st, $orthoIdNum_1_2);
foreach $orthoId(@orthoId){
	if(exists($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd}) and exists($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st})){
		$orthoIdNum_1_2++;
	}else{
		$orthoIdNum_1st++;
	}
}

#print "gene只分布在taxonId1st中的orthoId数量：$orthoIdNum_1st\n";
#print "gene分布在taxonId1st和taxonId2nd的orthoId数量：$orthoIdNum_1_2\n";



## 计算所有orthoId中的gene表达模式的平均相关系数
my (@totalGeneId);
my $finalTsv = "$outputDir/$jobId.$taxonId1st.$taxonId2nd.final.complete.cross.compare.group.mean.tpm.list";
open WW, ">$finalTsv";
print WW join("\t", "orthoId", "geneId_1st", "meanTpmList1", "geneId_2nd", "meanTpmList2") . "\n";
foreach $orthoId(@orthoId){
	@totalGeneId = ();
	# taxonId1st
	@geneId = ();
	@geneId = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st});
	push(@totalGeneId, @geneId);
	# taxonId2nd
	if(exists($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd})){
		@geneId = ();
		@geneId = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd});
		push(@totalGeneId, @geneId);
	}

	# 列出该orthoId下所有两个geneId之间的比较
	for(my $i=0; $i<=$#totalGeneId; $i++){
		for(my $j=$i+1; $j<=$#totalGeneId; $j++){

			# 列出geneId[$i]在所有group中tpm均值
			$output1 = "";
			foreach $groupName(@group){
				$output1 .= $geneToGroupToTpmHref->{$totalGeneId[$i]}->{$groupName} . ",";
			}

			# 列出geneId[$j]在所有group中tpm均值
			$output2 = "";
			foreach $groupName(@group){
				$output2 .= $geneToGroupToTpmHref->{$totalGeneId[$j]}->{$groupName} . ",";
			}

			# 输出整个比较行
			print WW join("\t", $orthoId, $totalGeneId[$i], substr($output1, 0, length($output1) - 1), 
					$totalGeneId[$j], substr($output2, 0, length($output2) - 1)
					) . "\n";

		}
	}
}
close WW;

# 利用R语言计算相关系数和pvalue值
# orthoId geneId_1st meanTpmList geneId_2nd  meanTpmList corCoef corPvalue
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

# 计算cross、taxon1st和taxon2nd中gene表达模式平均相关系数
my (%orthoToTaxonCoef, $orthoToTaxonCoefHref, $geneId1, $geneId2, $tmpList1, $tpmList2, $coef, $pvalue);
$orthoToTaxonCoefHref = \%orthoToTaxonCoef;
open FF, "<$finalTsv";
# orthoId geneId_1st meanTpmList1 geneId_2nd meanTpmList2 corCoef corPvalue
<FF>;
while($line=<FF>){
	chomp($line);
	($orthoId, $geneId1, $tmpList1, $geneId2, $tpmList2, $coef, $pvalue) = split(/\t/, $line);
	# 属于taxonId1st内的比较
	if($geneToTaxon{$geneId1} eq $geneToTaxon{$geneId2} and $geneToTaxon{$geneId2} eq $taxonId1st){
		$orthoToTaxonCoefHref->{$orthoId}->{$taxonId1st} .= $coef . ",";
	# 属于taxonId2n内的比较
	}elsif($geneToTaxon{$geneId1} eq $geneToTaxon{$geneId2} and $geneToTaxon{$geneId2} eq $taxonId2nd){
		$orthoToTaxonCoefHref->{$orthoId}->{$taxonId2nd} .= $coef . ",";
	# 属于taxonId1st和taxonId2nd间的比较
	}elsif($geneToTaxon{$geneId1} ne $geneToTaxon{$geneId2}){
		$orthoToTaxonCoefHref->{$orthoId}->{"cross"} .= $coef . ",";
	}
}
close FF;
system("rm -rf $finalTsv");

# 重新输出每个orthoId，并计算平均值
my ($outputLine, @coef, $geneNum, $avgCoef, $coefList);
#my $finalRlt = $outputDir . "/$jobId.$taxonId1st.$taxonId2nd.final.ortho.coexpress.result.tsv";
@orthoId = ();
@orthoId = keys(%orthoToTaxonCoef);
open WW, ">$outputDir/$outputCoexpressOrtho";
print WW join("\t", "orthoId", $taxonId1st . "_geneList", $taxonId1st . "_geneNum", $taxonId1st . "_coef", $taxonId2nd . "_geneList", $taxonId2nd . "_geneNum", $taxonId2nd . "_coef", "cross_coef") . "\n";
foreach $orthoId(@orthoId){
	$outputLine = $orthoId;
	# taxonId1st
	$outputLine .= "\t" . substr($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st}, 0, length($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st})-1);
	@geneId = ();
	@geneId = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId1st});
	$geneNum = $#geneId + 1;
	$outputLine .= "\t" . $geneNum;
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
	if(not(exists($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd}))){
		$outputLine .= "\t-\t-\t-";
	}else{
		$outputLine .= "\t" . substr($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd}, 0, length($orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd}) - 1);		
		@geneId = ();
		@geneId = split(/,/, $orthoToTaxonToGeneIdListHref->{$orthoId}->{$taxonId2nd});
		$geneNum = $#geneId + 1;
		$outputLine .= "\t" . $geneNum;
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


=beg
# 为echart生成数据文件
my ($orthoId, $geneList1st, $geneNum1st, $coef1st, $geneList2nd, $geneNum2nd, $coef2nd, $coefCross);
my ($taxon1stOrthoList, $taxon2ndOrthoList, $y);
my $finalEchart=$outputDir . "/$jobId.$taxonId1st.$taxonId2nd.echart.orth.gene.coexpress.tsv";
open WW, ">$finalEchart";
print WW "[\n";
open FF, "<$finalRlt";
# orthoId 	1stGeneList  		   1stGeneNum 39947_coef        2ndGeneList 2ndGeneNum  2nd_coef  cross_coef
# OG0014369     Os01g0674400,Os05g0560500  2          0.823835913645117 -           -           -         0.823835913645117 
<FF>;
while($line=<FF>){
	chomp($line);
	($orthoId, $geneList1st, $geneNum1st, $coef1st, $geneList2nd, $geneNum2nd, $coef2nd, $coefCross) = split(/\t/, $line);

	next if($coefCross <= 0.5);

	$coefCross = int($coefCross * 10000) - 5000;
	
	if($coef1st eq "-"){
		$coef1st = 0.9;
	}
	$taxon1stOrthoList .= "[$coefCross, $coef1st, $geneNum1st, \'$orthoId\', $taxonId1st],";

	if($geneList2nd ne "-"){
		if($coef2nd eq "-"){
			$coef2nd = 0.8;
		}
		$taxon2ndOrthoList .= "[$coefCross, $coef2nd, $geneNum2nd, \'$orthoId\', $taxonId2nd],";
	}
}
close FF;

print WW "[" . substr($taxon1stOrthoList, 0, length($taxon1stOrthoList) -1) . "],\n";
print WW "[" . substr($taxon2ndOrthoList, 0, length($taxon2ndOrthoList) -1) . "],\n";
print WW "]";
close WW;
=cut;
