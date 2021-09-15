#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--trFasta \\\n" .
		"--trRlt \\\n" .
		"--spFasta \\\n" .
		"--spRlt \\\n" .
                "--geneToTrspt \\\n" .
                "--outputGeneAnno \n";
	exit;
}

my ($trRlt, $spRlt, $trFasta, $spFasta, $geneToTrspt, $outputGeneAnno);

GetOptions(
	'trFasta=s'=>\$trFasta,
        'trRlt=s'=>\$trRlt,
	'spFasta=s'=>\$spFasta,
        'spRlt=s'=>\$spRlt,
        'geneToTrspt=s'=>\$geneToTrspt,
        'outputGeneAnno=s'=>\$outputGeneAnno,
);

my (%trsptToGene, %geneToSwissId, $line, @field, $field, @tt, $swissId, $geneId, $firstSpacePos);

# 获得tr中symbol和fulllName
my (%tr, $trHref, $trId, $trSymbol, $trFullName);
$trHref=\%tr;
open FF, "<$trFasta";
# >tr|K4CJ08|K4CJ08_SOLLC Nodulin26-like intrinsic protein 51 OS=Solanum lycopersicum OX=4081 GN=101244210 PE=2 SV=1
while($line=<FF>){
	chomp($line);
	next if(not($line=~/>/));
	@field = ();
	@field = split(/ /, $line);

	# 提取trId
	@tt = ();
	@tt = split(/\|/, $field[0]);
	$trId = $tt[2];

	# 获得tr的symbol
	$trSymbol = "NA";
	foreach $field(@field){
		if($field=~/GN=(.*)/){
			$trSymbol = $1;
		}
	}
	
	# 获得tr的全名
	@field = ();
	@field = split(/OS=/, $line);
	$firstSpacePos = index($field[0], " ");
	$trFullName = substr($field[0], $firstSpacePos+1);

	$trHref->{$trId}->{"symbol"} = $trSymbol;
	$trHref->{$trId}->{"fullName"} = $trFullName;
}
close FF;

# 读取sp序列中fullName和symbol
my (%sp, $spHref, $spId, $spSymbol, $spFullName);
$spHref=\%sp;
open FF, "<$spFasta";
# >sp|Q9FVN0|AMT13_SOLLC Ammonium transporter 1 member 3 OS=Solanum lycopersicum OX=4081 GN=AMT1-3 PE=2 SV=1
while($line=<FF>){
	chomp($line);
	next if(not($line=~/>/));
	@field = ();
	@field = split(/ /, $line);

	# 提取spId
	@tt = ();
	@tt = split(/\|/, $field[0]);
	$spId = $tt[2];

	# 获得sp的symbol
	$spSymbol = "NA";
	foreach $field(@field){
		if($field=~/GN=(.*)/){
			$spSymbol = $1;
		}
	}
	
	# 获得sp的全名
	@field = ();
	@field = split(/OS=/, $line);
	$firstSpacePos = index($field[0], " ");
	$spFullName = substr($field[0], $firstSpacePos+1);

	$spHref->{$spId}->{"symbol"} = $spSymbol;
	$spHref->{$spId}->{"fullName"} = $spFullName;

}
close FF;

# 读取基因组中transcript对弈的symbol和geneId
my (%annoTrsptToGeneId, %annoGeneToSymbol);
open FF, "<$geneToTrspt";
# >Solyc04g072140.3.1 transcript_name:NA gene_id:Solyc04g072140.3 gene_name:NA
while($line=<FF>){
	chomp($line);
	if($line=~/(.*)\t(.*)/){
		#print join("\t", $2, $1);
		#<STDIN>;
		$annoTrsptToGeneId{$2} = $1;
		$annoGeneToSymbol{$1} = "NA";
	}
}
close FF;


# 读取tr注释的结果
# 为每个tr蛋白找出对应geneId列表和evalue值： Q94FU1_SOLLC -> {geneId} = evalue
# 如果匹配到多个geneId片段，那么只保留最小evalue的那个
my (%geneToTr, $geneToTrHref, $trId, $trsptId, $geneId);
$geneToTrHref=\%geneToTr;
#tr|Q94FU1|Q94FU1_SOLLC  Solyc06g062920.3.1      100.000 464     0       0       1       464     1       464     0.0     966
#tr|K4CPZ0|K4CPZ0_SOLLC  Solyc09g005000.1.1      100.000 674     0       0       1       674     1       674     0.0     1392
#sp|Q6F473|RS29_PLUXY    XM_011552265.1  100.000 56      0       0       1       56      1       56      5.68e-38        119
#sp|Q8I7U0|RS13_PLUXY    XM_011559990.1  99.338  151     1       0       1       151     1       151     1.07e-108       305
open FF, "<$trRlt";
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);

	# 获得trId
	@tt = ();
	@tt = split(/\|/, $field[0]);
	$trId = $tt[2];

	# 获得trsptId和geneId
	$trsptId = $field[1];
	$geneId = $annoTrsptToGeneId{$trsptId};

	# 登记tr匹配的geneId的最小evalue
	if(not exists($geneToTrHref->{$geneId}->{$trId})){
		$geneToTrHref->{$geneId}->{$trId} = $field[10];
	}elsif($field[10]<$geneToTrHref->{$geneId}->{$trId}){
		$geneToTrHref->{$geneId}->{$trId} = $field[10];
	}
#	print join("\t", $geneId, $trId) . "\n";
#	<STDIN>;
}
close FF;

# 上述过程是以蛋白组为数据库，为tr序列做注释
#      这样可能出现一个多个tr对应到1个蛋白组基因的情况
# 我们最终是想为蛋白组中的序列做注释，因此需要反置处理：
# 	为每个蛋白组基因只提供一个最佳tr序列
my (%finalGeneToTr, @geneId, @trId, $minEvalue, $minTrId);
my @geneId = keys(%geneToTr);
foreach $geneId(@geneId){
	@trId = ();
	next if(not exists($geneToTrHref->{$geneId}));
	@trId = keys(%{$geneToTrHref->{$geneId}});
	# 寻找最小evalue的trId
	$minEvalue = 10;
	foreach $trId(@trId){
		if($geneToTrHref->{$geneId}->{$trId} < $minEvalue){
			$minEvalue = $geneToTrHref->{$geneId}->{$trId};
			$minTrId = $trId;
		}
	}
	$finalGeneToTr{$geneId} = $minTrId;
#	print join("\t", $geneId, $trId) . "\n";
#	<STDIN>;
}


# 读取sp注释的结果
# 为每个sp蛋白找出对应geneId列表和evalue值： Q94FU1_SOLLC -> {geneId} = evalue
# 如果匹配到多个geneId片段，那么只保留最小evalue的那个
my (%geneToSp, $geneToSpHref, $spId, $trsptId);
$geneToSpHref=\%geneToSp;
open FF, "<$spRlt";
# sp|P29535|1A14_SOLLC    Solyc05g050010.3.1      99.580  476     2       0       1       476     1       476     0.0     987
# sp|Q9FVN0|AMT13_SOLLC   Solyc03g045070.1.1      100.000 460     0       0       1       460     1       460     0.0     924 
# tr|A0A1L8D6F4|A0A1L8D6F4_PLUXY  XM_011549556.1  98.643  516     7       0       21      536     1       516     0.0     1059
# tr|A0A1L8D6H7|A0A1L8D6H7_PLUXY  XM_011566673.1  100.000 333     0       0       1       333     1       333     0.0     699
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);

	# 获得spId
	@tt = ();
	@tt = split(/\|/, $field[0]);
	$spId = $tt[2];

	# 获得geneId
	$trsptId = $field[1];
	$geneId = $annoTrsptToGeneId{$trsptId};

	# 登记tr匹配的geneId的最小evalue
	if(not exists($geneToSpHref->{$geneId}->{$spId})){
		$geneToSpHref->{$geneId}->{$spId} = $field[10];
	}elsif($field[10]<$geneToSpHref->{$geneId}->{$spId}){
		$geneToSpHref->{$geneId}->{$spId} = $field[10];
	}
#	print join("\t", $geneId, $spId);
#	<STDIN>;
}
close FF;

# 上述过程是以蛋白组为数据库，为sp序列做注释
#      这样可能出现一个多个sp对应到1个蛋白组基因的情况
# 我们最终是想为蛋白组中的序列做注释，因此需要反置处理：
#       为每个蛋白组基因只提供一个最佳sp序列
my (%finalGeneToSp, @geneId, @spId, $minEvalue, $minSpId);
my @geneId = keys(%geneToSp);
foreach $geneId(@geneId){
#	print $geneId;
#	<STDIN>;
	next if(not exists($geneToSpHref->{$geneId}));
	@spId = ();
	@spId = keys(%{$geneToSpHref->{$geneId}});
	# 寻找最小evalue的spId
	$minEvalue = 10;
	$minSpId = "";
	foreach $spId(@spId){
#		print join("\t", "geneId:".$geneId, "spId:" . $spId, "evalue:" . $geneToSpHref->{$geneId}->{$spId});
#		<STDIN>;
		if($geneToSpHref->{$geneId}->{$spId} < $minEvalue){
			$minEvalue = $geneToSpHref->{$geneId}->{$spId};
			$minSpId = $spId;
		}
	}
	$finalGeneToSp{$geneId} = $minSpId;
#	print join("\t", "geneId:".$geneId, "minSpId:" . $minSpId) . "\n";
#	<STDIN>;
}


# 输出所有蛋白组中的基因注释
my ($trSymbol, $trFullName, $spSymbol, $spFullName);
my @geneId = keys(%annoGeneToSymbol);
open WW, ">$outputGeneAnno";
print WW join("\t", "geneId", "annoSymbol", "trId", "trSymbol", "trFullName", "spId", "spSymbol", "spFullName") . "\n";
foreach $geneId(@geneId){
	if(exists($finalGeneToSp{$geneId})){
		$spId = $finalGeneToSp{$geneId};
		$spSymbol = $spHref->{$spId}->{"symbol"};
		$spFullName = $spHref->{$spId}->{"fullName"};
	}else{
		$spId = "NA";
		$spSymbol = "NA";
		$spFullName = "NA";
	}
	if(exists($finalGeneToTr{$geneId})){
		$trId = $finalGeneToTr{$geneId};
		$trSymbol = $trHref->{$trId}->{"symbol"};
		$trFullName = $trHref->{$trId}->{"fullName"};
	}else{
		$trId = "NA";
		$trSymbol = "NA";
		$trFullName = "NA";
	}

	print WW join("\t", $geneId, $annoGeneToSymbol{$geneId}, $trId, $trSymbol, $trFullName, $spId, $spSymbol, $spFullName) . "\n";
}
close WW;
