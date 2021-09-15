#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonIdList \\\n" .
                "--baseDir \\\n" .
                "--outputTsv\n";
	exit;
}

my ($taxonIdList, $baseDir, $outputTsv);

GetOptions(
        'taxonIdList=s'=>\$taxonIdList,
        'baseDir=s'=>\$baseDir,
        'outputTsv=s'=>\$outputTsv,
);



my (@taxonId, $taxonId, @name, @value, %tmpHash, $line, %study, $experimentNum, $baseVolume, $runNum, @run, $tissueNum, %tissue);
my (@study, @tissue);

open FF, "<$taxonIdList";
@taxonId = <FF>;
close FF;
open WW, ">$outputTsv";
print WW join("\t", "分类号", "研究数", "实验数", "测序RUN", "组织数", "总测序容量", "优化前剪接位点数", "优化后剪接位点数", "优化剪接位点增量", "优化前外显子数", "优化后外显子数", "优化外显子增量", "单外显子基因数", "多外显子基因数", "优化前多外显子基因上转录本总数", "优化前多外显子基因上转录本平均数", "优化后多外显子基因上转录本总数", "优化后多外显子基因上转录本平均数", "优化前具有多转录本的多外显子基因数", "优化前具有多转录本的多外显子基因占比", "优化后有多转录本的多外显子基因数", "优化后具有多转录本的多外显子基因占比") . "\n";
foreach $taxonId(@taxonId){
	chomp($taxonId);


	my $origSpliceNum = `cat $baseDir/$taxonId/001-prepare-local-datasource/splicesites.txt |wc -l`;
	chomp($origSpliceNum);
	my $improvedSpliceNum = `cat $baseDir/$taxonId/007-build-genomeDb-on-final-anno/splicesites.txt |wc -l`;
	chomp($improvedSpliceNum);
	my $spliceIncrementPercent=sprintf("%.2f", ($improvedSpliceNum-$origSpliceNum)/$origSpliceNum*100);


####  == 1 sra 统计信息 == #########

	# 采集用于构建数据库的sra数据情况
	my $filterExperimentTsv = $baseDir . "/" . $taxonId . "/010-gather-alignment-info-of-all-expts/filtered.alignment.info.of.assembled.experiment.tsv";
	# 清空studyhash
	%study = ();
	# 清空experiment数量
	$experimentNum = 0;
	# 清空测序数据容量
	$baseVolume = 0;
	# 清空tissue
	%tissue = ();
	# 清空run个数
	$runNum = 0;
	
	@name = ();
	@value = ();
	open FF, "<$filterExperimentTsv";
	$line = <FF>;
	chomp($line);
	@name = split(/\t/, $line);
	while($line=<FF>){
		@value = split(/\t/, $line);
		for(my $i=0; $i<=$#value; $i++){
			$tmpHash{$name[$i]}=$value[$i];
		}
		$study{$tmpHash{"Study"}} = 1;
		$experimentNum++;
		$baseVolume+=$tmpHash{"Base"};
		@run = ();
		@run = split(/,/, $tmpHash{"RunList"});
		$runNum += $#run+1;
		$tissue{$tmpHash{"Tissue"}} = 1;
	}
	close FF;
	
	@study = keys(%study);
	@tissue = keys(%tissue);



##### == 2 exon 优化统计信息 == #######
	my $origExonNum=`grep -P \"\\texon\\t\" $baseDir/$taxonId/001-prepare-local-datasource/ensembl.gtf | awk -F \'\\t\' \'{print \$1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$7}\' |sort -u |wc -l`;
	chomp($origExonNum);
	my $improvedExonNum=`grep -P \"\\texon\\t\" $baseDir/$taxonId/006-form-final-trspt-annotation/final.complete.trspt.anno.gtf |awk -F \'\\t\' \'{print \$1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$7}\' |sort -u |wc -l`;
	chomp($improvedExonNum);
	my $improvedPercentage = sprintf("%.2f",($improvedExonNum-$origExonNum)/$origExonNum*100);


##### == 3 多外显子基因上平均转录本数量 ######
	# 基因总数
	my $origGeneNum = `grep -cP \"\\tgene\\t\" $baseDir/$taxonId/001-prepare-local-datasource/ensembl.gtf`;
	chomp($origGeneNum);

	# 原始注释单外显子基因数量
	my $origSingexonGeneNum = `grep -P \"\\texon\\t\" $baseDir/$taxonId/001-prepare-local-datasource/ensembl.gtf |awk -F \'\\t\' \'{print \$9}\' | awk -F \';\' \'{print \$1}\' |sort |uniq -c |grep -c \" 1 \"`;
	chomp($origSingexonGeneNum);

	# 原始注释中，多外显子基因数量
	my $multiExonGeneIdList = "./multi-exon-geneId.list";
	`grep -P \"\\texon\\t\" $baseDir/$taxonId/001-prepare-local-datasource/ensembl.gtf |awk -F \'\\t\' \'{print \$9}\' | awk -F \';\' \'{print \$1}\' |sort |uniq -c |grep -v \" 1 \" |awk -F \'gene_id\' \'{print \$2}\' > $multiExonGeneIdList`;
	my $origMultiExonGeneNum = `cat $multiExonGeneIdList |wc -l`;
	chomp($origMultiExonGeneNum);

	# 原始注释中多外显子基因上transcript数量
	my $origTrsptNumInMultExonGene=`grep -Ff $multiExonGeneIdList $baseDir/$taxonId/001-prepare-local-datasource/ensembl.gtf |grep -cP \"\\ttranscript\\t\"`;
	chomp($origTrsptNumInMultExonGene);

	# 原始注释中多外显子基因上平均transcript数量
	my $origAvgTrsptNumInPerMultExonGene = sprintf("%.2f", $origTrsptNumInMultExonGene/$origMultiExonGeneNum);

	# 原始注释中，具有多个trspt的多外显子基因数量
	my $origMultiExonGeneNumWithMultTrspt = `grep -Ff $multiExonGeneIdList $baseDir/$taxonId/001-prepare-local-datasource/ensembl.gtf |grep -P \"\\ttranscript\\t\" |awk -F \'\\t\' \'{print \$9}\' |awk -F \';\' \'{print \$1}\' | sort |uniq -c |grep -cv \" 1 \"`;
	chomp($origMultiExonGeneNumWithMultTrspt);

	# 原始注释中，具有多个trspt的多外显子基因占比
	my $origMultiExonGenePercentWithMultTrspt = sprintf("%.2f", $origMultiExonGeneNumWithMultTrspt/$origMultiExonGeneNum*100);

	# 优化后注释中多外显子基因上的trspt总数量
	my $improvedTrsptNumInMultExonGene = `grep -Ff $multiExonGeneIdList $baseDir/$taxonId/006-form-final-trspt-annotation/final.complete.trspt.anno.gtf |grep -cP \"\\ttranscript\\t\"`;
	chomp($improvedTrsptNumInMultExonGene);

	# 优化后注释中多外显子基因上的平均trspt数量
	my $improvedAvgTrsptNumInPerMultExonGene = sprintf("%.2f", $improvedTrsptNumInMultExonGene/$origMultiExonGeneNum);

	# 优化后，具有多个trspt的多外显子基因数量
	my $improvedMultiExonGeneNumWithMultTrspt=`grep -Ff $multiExonGeneIdList $baseDir/$taxonId/006-form-final-trspt-annotation/final.complete.trspt.anno.gtf |grep -P \"\\ttranscript\\t\" |awk -F \'\\t\' \'{print \$9}\' |awk -F \';\' \'{print \$1}\' | sort |uniq -c |grep -cv \" 1 \"`;
	chomp($improvedMultiExonGeneNumWithMultTrspt);
	
	# 优化后，具有多个trspt的多外显子基因占比
	my $improvedMultiExonGenePercentWithMultTrspt = sprintf("%.2f", $improvedMultiExonGeneNumWithMultTrspt/$origMultiExonGeneNum*100);
	
	


	print WW join("\t", $taxonId, $#study+1, $experimentNum, $runNum, $#tissue+1, sprintf("%.2f", $baseVolume), $origSpliceNum, $improvedSpliceNum, $spliceIncrementPercent, $origExonNum, $improvedExonNum, $improvedPercentage, $origSingexonGeneNum, $origMultiExonGeneNum, $origTrsptNumInMultExonGene, $origAvgTrsptNumInPerMultExonGene, $improvedTrsptNumInMultExonGene, $improvedAvgTrsptNumInPerMultExonGene, $origMultiExonGeneNumWithMultTrspt, $origMultiExonGenePercentWithMultTrspt, $improvedMultiExonGeneNumWithMultTrspt, $improvedMultiExonGenePercentWithMultTrspt) . "\n";

}
