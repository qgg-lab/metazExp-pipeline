#!/usr/bin/perl
#use strict;
use Getopt::Long;
use List::Util qw/sum/;
use Statistics::R;

if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonId 39947 \\\n" .
                "--jobId 1234 \\\n" .
                "--exptDir ../013/expts/gene \\\n" .
                "--exptGroupInfo sample.info.tsv \\\n" .
                "--statTest Yes \\\n" .
                "--foldChange 2\\\n" .
                "--pvalue 0.05 \\\n" .
                "--minReplicateNum 3\\\n" .
                "--minGroupNum 5 \\\n" .
		"--outputDir specificityByTpmDir\n";
	exit;
}

my ($taxonId, $jobId, $exptDir, $exptGroupInfo, 
$statTest, $foldChange, $pvalue, 
$minReplicateNum, $minGroupNum, $outputDir);

GetOptions(
        'taxonId=s'=>\$taxonId,
        'jobId=s'=>\$jobId,
        'exptDir=s'=>\$exptDir,
        'exptGroupInfo=s'=>\$exptGroupInfo,
	'statTest=s'=>\$statTest,
	'foldChange=s'=>\$foldChange,
	'pvalue=s'=>\$pvalue,
	'minReplicateNum=s'=>\$minReplicateNum,
	'minGroupNum=s'=>\$minGroupNum,
	'outputDir=s'=>\$outputDir,
);

#print join("\t", $taxonId, $jobId, $exptDir, $exptGroupInfo, $foldChange, $pvalue, $minReplicateNum, $minGroupNum, $outputDir) . "\n";
system("mkdir -p $outputDir");
# 把exptId读入到group对应的hash中
# root->SRX1203042,SRX1203043,SRX1203045
# leaf->SRX1203049,SRX1203050,SRX1203055
my (%groupToExptIdList, %tmpHash, @fieldName, @fieldValue, $line);
open FF, "<$exptGroupInfo";
my $line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@fieldValue = split(/\t/, $line);
	for(my $i=0; $i<=$#fieldValue; $i++){
		$tmpHash{$fieldName[$i]} = $fieldValue[$i];
	}
	$groupToExptIdList{$tmpHash{"groupName"}} .= $tmpHash{"exptId"} . ",";
}
close FF;

# 获得所有group存入数组，同时去掉exptId列表中末尾的逗号
my (@group, $group);
@group = keys(%groupToExptIdList);
for(my $i=0; $i<=$#group; $i++){
	$groupToExptIdList{$group[$i]} = substr($groupToExptIdList{$group[$i]}, 0, length($groupToExptIdList{$group[$i]}) - 1);
}

# 将gene在每个experiment中表达水平TPM按照group进行分类，登记到每个gene的hash中
my (@exptId, $exptId, @expr, $geneId, $tpm);
my ($exprFile);
my (%geneExprListInGroup, $geneExprListInGroupHref);
$geneExprListInGroupHref = \%geneExprListInGroup;
foreach $group(@group){
	@exptId = ();
	@exptId = split(/,/, $groupToExptIdList{$group});
	foreach $exptId(@exptId){
		$exprFile = $exptDir . "/" . $exptId;
		open FF, "<$exprFile";
		# geneId  	readCount       Coverage        FPKM    	TPM
		# Os02g0307300  0       	0.005889        0.001047        0.002106
		<FF>;
		while($line = <FF>){
			chomp($line);
			@expr = split(/\t/, $line);
			$geneId = $expr[0];
			$tpm = $expr[4];
			$geneExprListInGroupHref->{$geneId}->{$group} .= $tpm . ",";
		}
		close FF;
	}
}

# 将gene的表达量列表按照group分组输出到文件
my (@geneId, $outputLine);
@geneId = keys(%geneExprListInGroup);
my $exprByGroupFile = $outputDir . "/$jobId.$taxonId.tpm.list.in.group.tsv";
open WW, ">$exprByGroupFile";
print WW join("\t", "geneId", @group) . "\n";
foreach $geneId(@geneId){

	$outputLine = $geneId;

	for(my $i=0; $i<=$#group; $i++){

		$geneExprListInGroupHref->{$geneId}->{$group[$i]} = substr($geneExprListInGroupHref->{$geneId}->{$group[$i]}, 0, length($geneExprListInGroupHref->{$geneId}->{$group[$i]}) - 1);

		$outputLine .= "\t" . $geneExprListInGroupHref->{$geneId}->{$group[$i]};
	}
	
	print WW $outputLine . "\n";
}
close WW;

my $tmpRlt = $outputDir . "/$jobId.$taxonId.diff.expr.by.tpm.tsv";
my $R=Statistics::R->new(shared => 1);
my $cmd=<<EOF;
df<-read.table(file="$exprByGroupFile", sep="\\t", quote ="", header=T)
geneId <- df\$geneId
expr<-df[,-1]
groupNum = dim(expr)[2]
groupName=colnames(expr)
rlt <- data.frame(a=c(1:length(geneId)))
for(i in c(1:groupNum )){
  if(i==groupNum)
    break
  for(j in c((i+1):groupNum)){
    compareName <- paste("FC", groupName[i], groupName[j],sep="/")
    rlt[,compareName] = apply(expr,1,function(x) mean(as.numeric(unlist(strsplit(x[i], ","))))/mean(as.numeric(unlist(strsplit(x[j], ",")))))
  }
}
rlt<-rlt[,-1]

for(i in c(1:groupNum )){
   if(i==groupNum)
      break
   for(j in c((i+1):groupNum)){
	compareName <- paste("adjPvalue", groupName[i],groupName[j],sep="/")
        rlt[,compareName] = p.adjust(apply(expr,1,function(x) t.test(as.numeric(unlist(strsplit(x[i], ","))),as.numeric(unlist(strsplit(x[j], ","))), alternative='two.sided', paired = F)\$p.value), method="fdr")
   }
}
rownames(rlt)=geneId
write.table(rlt, file='$tmpRlt', sep='\\t', quote=F, row.names = T)
EOF
$R->run($cmd);
$R->stop();

# 开始筛选
# geneId spGroup   spDirection   meanTpm   minPvalue maxPvalue minFC maxFC tpmListInAllGroup
# Os01   root      high          12.3456   0.01      0.04      2     4     root:10.3,12.9,3.8;leaf:12.4,23.8,48.5;
my (@tpm, $meanTpm, $tpmListInAllGroups, $tmpMinFC, $tmpMaxFC);
my ($line, @nameField, @valueField, $geneId);
my (%foldchange, %pvalue, $groupX, $groupY);
my ($minPvalue, $maxPvalue, $minFC, $maxFC, $spDirection, $giveUp);
open WW, ">$outputDir/$jobId.$taxonId.final.specific.expr.gene.list.by.tpm.with.t-test.tsv";
print WW join("\t", "geneId", "spGroup", "spDirection", "meanTpm", "minPvalue", "maxPvalue", "minFC", "maxFC", "tpmListInAllGroup") . "\n";
open FF, "<$tmpRlt";
$line = <FF>;
# "geneId无" FC/callus/panicle  FC/callus/leaf ... adjPvalue/callus/endosperm adjPvalue/panicle/leaf ..
chomp($line);
my @nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	$geneId = shift(@valueField);

	%foldchange = ();
	%pvalue = ();
	# 为每个group生成和其它所有group比的结果：
	#    (1)FC: 每个FC/groupY/groupX转换成2个groupX_groupY和groupY_groupX存放到FC
	#    (2)PVALUE: 每个adjPvalue/groupY/groupX转换成2个groupX_groupY和groupY_groupX存放到PVALUE
	for(my $i=0; $i<=$#valueField; $i++){
		if($nameField[$i]=~/FC\/(.*)\/(.*)/){
			$groupX = $1;
			$groupY = $2;
			$foldchange{$groupX . "/" . $groupY} = $valueField[$i];
			if($valueField[$i] ne "NA"){
				if($valueField[$i] != 0){
					# 即使$valueField[$i] == Inf，也可以执行1/$valueField[$i]运算，结果为0
					$foldchange{$groupY . "/" . $groupX} = 1/$valueField[$i];
				}else{
					$foldchange{$groupY . "/" . $groupX} = Inf;
				}
			}else{
				$foldchange{$groupY . "/" . $groupX} = "NA";
			}
		}elsif($nameField[$i]=~/adjPvalue\/(.*)\/(.*)/){
			$groupX = $1;
			$groupY = $2;
			$pvalue{$groupX . "/" . $groupY} = $valueField[$i];
			$pvalue{$groupY . "/" . $groupX} = $valueField[$i];
		}
	}

	# 对每个group检查:
	#    (1) 以它为开头的FC是否都大于$foldchange 或者 都小于1/$foldchange
	#    (2) 以它为开头的PVALUE是否都小于等于$pvalue
	foreach $groupX(@group){

		($minPvalue, $maxPvalue, $minFC, $maxFC, $spDirection, $giveUp, $tpmListInAllGroups) = (Inf, -Inf, Inf, -Inf, "-", "no", "");

		foreach $groupY(@group){

			next if($groupX eq $groupY);

			# 只要发现不能比较现象，那么放弃该group特异性检查: Inf, NA, 数值
			if($foldchange{$groupX . "/" . $groupY} eq "NA"){
				$giveUp = "yes";
				last;
			}
			# 只要发现不能比较现象，那么放弃该group特异性检查
			if($pvalue{$groupX . "/" . $groupY} eq "NA"){ # 只有NA和数值
				$giveUp = "yes";
				last;
			}

			# 更新$minFC, $maxFC: 注意: Inf,-Inf是数值，可以进行大小比较
			if($minFC > $foldchange{$groupX . "/" . $groupY}){
				$minFC = $foldchange{$groupX . "/" . $groupY};
			}
			if($maxFC < $foldchange{$groupX . "/" . $groupY}){
				$maxFC = $foldchange{$groupX . "/" . $groupY};
			}
			# 更新minPvalue, maxPvalue
			if($minPvalue > $pvalue{$groupX . "/" . $groupY}){
				$minPvalue = $pvalue{$groupX . "/" . $groupY};
			}
			if($maxPvalue < $pvalue{$groupX . "/" . $groupY}){
				$maxPvalue = $pvalue{$groupX . "/" . $groupY};
			}
		}

		# 如果检测到groupX和任意一个groupX不可以比（NA），那么认为该groupX肯定不是特异的
		next if($giveUp eq "yes");

		# $minPvalue和$maxPvalue都小于$pvalue 且 $minFC和$maxFC都大于$foldchange或者都小于1/$foldchange
		# 那么此时的groupX为特异的
		if(($minFC >= $foldChange or $maxFC <= 1/$foldChange) and $maxPvalue <= $pvalue){
			if($minFC > 1){
				$spDirection = "high";
			}else{
				$spDirection = "low";
				# $maxFC = 1/$minFC
				# $minFC = 1/$maxFC
				if($minFC ==0 ){
					$tmpMaxFC = Inf;
				}else{
					$tmpMaxFC = 1/$minFC;
				}
				if($maxFC ==0 ){
					$tmpMinFC = Inf;
				}else{
					$tmpMinFC = 1/$maxFC;
				}
				$maxFC = $tmpMaxFC;
				$minFC = $tmpMinFC;
			}

			# 求当前 groupX的tpm均值			
			@tpm = ();
			@tpm = split(/,/,$geneExprListInGroupHref->{$geneId}->{$groupX});
			$meanTpm = sum(@tpm)/($#tpm+1);

			# 将所有group的tpm都列出来
			$tpmListInAllGroups = "";
			foreach $groupY(@group){
				$tpmListInAllGroups .= $groupY . ":" . $geneExprListInGroupHref->{$geneId}->{$groupY} . ";";
			}
			print WW join("\t", $geneId, $groupX, $spDirection, sprintf("%.4f", $meanTpm), $minPvalue, $maxPvalue, $minFC, $maxFC, substr($tpmListInAllGroups, 0, length($tpmListInAllGroups)-1) ) . "\n";
		}
	}
}
close FF;
close WW;

