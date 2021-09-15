#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "\t\t--inclusionTagFile total.inclusion.with.tissue.expt.tsv\\\n" .
                "\t\t--skipingTagFile total.skiping.with.tissue.expt.tsv\\\n" .
                "\t\t--filterFile 00000.filter.by.sample.cfg \\\n" .
                "\t\t--validatedAsListFile total.as.with.validated.flag.tsv \n";
	exit;
}

my ($inclusionTagFile, $skipingTagFile, $filterFile, $validatedAsListFile);

GetOptions(
        'inclusionTagFile=s'=>\$inclusionTagFile,
        'skipingTagFile=s'=>\$skipingTagFile,
        'filterFile=s'=>\$filterFile,
        'validatedAsListFile=s'=>\$validatedAsListFile,
);

my ($line, @tmp, $i);
# 将过滤标准文件读入到hash中
my (%tissueToExptNum);
open FF, "<$filterFile";
# tissue  tissueCutoff
# ear     2
# embryo  3
# total       210
<FF>;
while($line=<FF>){
	chomp($line);
	@tmp = ();
	@tmp = split(/\t/, $line);	
	$tissueToExptNum{$tmp[0]}=$tmp[1];
}
close FF;


my (%asValidation);
my ($asId, $tissue, $exptList, @exptId, $flag, $tissueId,  $tissueExptList, $totalExptNum, @tt);

# 取inclusion中每个AS，查看其是否符合条件。如果不符合条件，那么标识为1，否则标识为0
open FF, "<$inclusionTagFile";
#ZMAYA3SS0000004719      total:10 seedling:SRX5794451,SRX5794453,SRX5794454       stem:SRX5387719,SRX5387724
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);
	$asId = $tmp[0];

	$flag = 0;
	# 检查总experiment数量
	@tt = ();
	@tt = split(/:/, $tmp[1]);
	if($tt[1]>=$tissueToExptNum{"total"}){
		$flag = 1;
		goto ASSIGN;
	}

	# 检查组织中expt数量是否通过
	for($tissueId=2; $tissueId<=$#tmp; $tissueId++){
		$tissueExptList = $tmp[$tissueId];
		@tt = ();
		@tt = split(/:/, $tissueExptList);
		$tissue = $tt[0];
		$exptList = $tt[1];
		@exptId = ();
		@exptId = split(/,/, $exptList);
		if(exists($tissueToExptNum{$tissue}) and $#exptId+1 >= $tissueToExptNum{$tissue}){
			$flag = 1;
			last;
		}
	}	
ASSIGN:
	${$asValidation{$asId}}{"inclusion"} = $flag;
}
close FF;

# 读取skiping中每个AS，查看其是否符合条件。如果不符合条件，那么标识为1，否则标识为0
open FF, "<$skipingTagFile";
#ZMAYA3SS0000004719      total:10 seedling:SRX5794451,SRX5794453,SRX5794454       stem:SRX5387719,SRX5387724
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);
	$asId = $tmp[0];

	$flag = 0;
	# 检查总experiment数量
	@tt = ();
	@tt = split(/:/, $tmp[1]);
	if($tt[1]>=$tissueToExptNum{"total"}){
		$flag = 1;
		goto ASSIGN;
	}

	# 检查组织中expt数量是否通过
	for($tissueId=2; $tissueId<=$#tmp; $tissueId++){
		$tissueExptList = $tmp[$tissueId];
		@tt = ();
		@tt = split(/:/, $tissueExptList);
		$tissue = $tt[0];
		$exptList = $tt[1];
		@exptId = ();
		@exptId = split(/,/, $exptList);
		if(exists($tissueToExptNum{$tissue}) and $#exptId+1 >= $tissueToExptNum{$tissue}){
			$flag = 1;
			last;
		}
	}	
ASSIGN:
	${$asValidation{$asId}}{"skiping"} = $flag;
}
close FF;

# 输出AS编号及其验证标识
open WW, ">$validatedAsListFile";
print WW "ASID\tinclusion\tskiping\tbothInclAndSkip\n";
my (@asId);
@asId = keys(%asValidation);
foreach $asId(@asId){
	print WW join("\t", $asId,  ${$asValidation{$asId}}{"inclusion"}, ${$asValidation{$asId}}{"skiping"}, ${$asValidation{$asId}}{"inclusion"}*${$asValidation{$asId}}{"skiping"}) . "\n";
}
close WW;
