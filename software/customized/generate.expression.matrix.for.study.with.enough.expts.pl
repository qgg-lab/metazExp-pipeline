#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw/sum/;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--studyExptListTsv 3702.study.exptIdList.tsv\\\n" .
		"--exptDir ../008-pickup-psi-of-ASs-in-all-expers/psiOutputDir/ \\\n" .
                "--tailVarianceGenePercentage 20 \\\n" .
                "--minNonZeroExprSamplePercentage 50 \\\n" .
		"--outputDir coExpression \n";
	exit;
}

my ($taxonId, $studyExptListTsv, $tailVarianceGenePercentage, $minNonZeroExprSamplePercentage, $outputDir, $exptDir, $geneIdList);

GetOptions(
        'studyExptListTsv=s'=>\$studyExptListTsv,
	'exptDir=s'=>\$exptDir,
        'tailVarianceGenePercentage=s'=>\$tailVarianceGenePercentage,
        'minNonZeroExprSamplePercentage=s'=>\$minNonZeroExprSamplePercentage,
	'outputDir=s'=>\$outputDir,
);

my ($geneLine, $line, $exptNum, $exptIdList, @exptId, $exptId, $study, $expressionFile, $geneId, $tpm, %geneId, @geneId);
my (@fieldName, @fieldValue, $fieldName, $fieldValue);
my (%studyExptGeneExpress, $studyExptGeneExpressHref, @study);
$studyExptGeneExpressHref = \%studyExptGeneExpress;

# 搜集表达数据整理进入hash
open FF, "<$studyExptListTsv";
# 3702    SRP082532       35      SRX2039036,SRX2039037,SRX2039038,SRX2039039,SRX2039040,SRX2039041,SRX2039042
while($line=<FF>){
	chomp($line);
	($taxonId, $study, $exptNum, $exptIdList) =  split(/\t/, $line);

	@exptId = ();
	@exptId = split(/,/, $exptIdList);
	foreach $exptId(@exptId){
		$expressionFile = $exptDir . "/" . $exptId;
		open EE, "<$expressionFile";
		# geneId  	readCount   Coverage        FPKM       TPM
		# Os04g0492900  8           0.188577        2.221850   4.092893
		$geneLine = <EE>;
		while($geneLine=<EE>){
			chomp($geneLine);
			@fieldValue = ();
			@fieldValue = split(/\t/, $geneLine);
			$geneId = $fieldValue[0];
			$geneId{$geneId} = 1;
			$tpm = $fieldValue[$#fieldValue];
			$studyExptGeneExpressHref->{$study}->{$exptId}->{$geneId} = $tpm;
		}
		close EE;
	}
}
close FF;

# 按照study输出表达数据
my ($exptNumCutoff, $exptNum, $geneExpressList, @tmpTpm, $tpm, $tpmList, $nonZeroTpmNumTag, $variance, @tpmArr, $avg);
@geneId = keys(%geneId);
@study = keys(%studyExptGeneExpress);

foreach $study(@study){

	# 获得当前study下的所有exptId，存放在数组中
	@exptId = ();
	@exptId = keys(%{$studyExptGeneExpressHref->{$study}});
	$exptNumCutoff = int(($#exptId+1)*$minNonZeroExprSamplePercentage/100);

	@tpmArr = ();
	# 收集geneID、表达值列表、tpm个数是否达标、表达量方差存放在二维数组@tpmArr
	for(my $i=0; $i<=$#geneId; $i++){
		
		# 将gene在所有expt中的表达值搜集起来，并统计大于0表达值的个数
		$exptNum = 0;
		$tpmList = "";
		$nonZeroTpmNumTag = 0;
		$variance = -1;
		for(my $j=0; $j<=$#exptId; $j++){
			if(exists($studyExptGeneExpressHref->{$study}->{$exptId[$j]}->{$geneId[$i]}) and $studyExptGeneExpressHref->{$study}->{$exptId[$j]}->{$geneId[$i]}>0){
				$tpmList .= $studyExptGeneExpressHref->{$study}->{$exptId[$j]}->{$geneId[$i]} . ",";
				$exptNum++;
			}else{
				$tpmList .= "0,";
			}
		}
		
		# 放弃在所有experiment中都没有检测到gene
		next if($exptNum==0);

		# 判断当前大于0表达值个数是否达标
		if($exptNum >= $exptNumCutoff){
			$nonZeroTpmNumTag = "Y";
		}else{
			$nonZeroTpmNumTag = "N";
		}

		# 计算当前基因表达值的方差
		@tmpTpm = ();
		$tpmList = substr($tpmList, 0, length($tpmList) -1);
		@tmpTpm = split(/,/, $tpmList);
		$variance = 0;
		$avg = sum(@tmpTpm)/($#tmpTpm+1);
		foreach $tpm(@tmpTpm){
			$variance += ($tpm-$avg)*($tpm-$avg);
		}
		$variance=$variance/($#tmpTpm+1);

		# 将当前gene、tpmList、nonZeroTpmNumTag、variance存储到二维数组中便于后续排序
		($tpmArr[$i][0], $tpmArr[$i][1], $tpmArr[$i][2], $tpmArr[$i][3]) = ($geneId[$i], $tpmList, $nonZeroTpmNumTag, $variance);		
	}

	# 按照方差升序排序二维数组
	@tpmArr = sort{$a->[3]<=> $b->[3]} @tpmArr;

	# 计算删除gene的终止编号
	my $k = int(($#tpmArr+1)*$tailVarianceGenePercentage/100);

	open WW, ">$outputDir" . "/" . $study . ".geneExpr.matrix.xls";
	#print WW join("\t", "GeneID", "NonZeroTpm", "variance", "PASS", @exptId) . "\n";
	print WW join("\t", "GeneID", @exptId) . "\n";
	for(my $t=0; $t<$#tpmArr; $t++){
		my @tt = ();
		@tt = split(/,/, $tpmArr[$t][1]);
		#my $ss = join("\t", @tt);
		if($t<$k or $tpmArr[$t][2] eq "N"){
			#print WW join("\t", $tpmArr[$t][0], $tpmArr[$t][2], $tpmArr[$t][3], "delete", @tt) . "\n";
			#print WW join("\t", $tpmArr[$t][0], @tt) . "\n";
		}else{
			#print WW join("\t", $tpmArr[$t][0], $tpmArr[$t][2], $tpmArr[$t][3], "retain", @tt) . "\n";
			print WW join("\t", $tpmArr[$t][0], @tt) . "\n";
		}
	}
	close WW;
}

