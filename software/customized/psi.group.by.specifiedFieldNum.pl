#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--specifiedField \\\n" .
                "--groupPctgNum \\\n" .
                "--groupSize \\\n" .
                "--psiMysqlTsv \\\n" .
                "--filterValid \\\n" .
		"--outputPsiTsvForViolin \\\n" .
		"--outputAsNumInGroup \n";
	exit;
}

my ($specifiedField, $groupPctgNum, $groupSize, $psiMysqlTsv, $filterValid, $outputPsiTsvForViolin, $outputAsNumInGroup);

GetOptions(
        'specifiedField=s'=>\$specifiedField,
        'groupPctgNum=s'=>\$groupPctgNum,
        'groupSize=s'=>\$groupSize,
	'psiMysqlTsv=s'=>\$psiMysqlTsv,
	'filterValid=s'=>\$filterValid,
	'outputPsiTsvForViolin=s'=>\$outputPsiTsvForViolin,
	'outputAsNumInGroup=s'=>\$outputAsNumInGroup,
);


my (@name, @value, $nameList, $valueList, $name, $value, %tmpPsi);
my (%psi, $psiHref, %totalSpecifiedField);
$psiHref=\%psi;

# 扫描所有PSI值，以ASID为关键字搜集AS对应的PSI和Tissue列表
open FF, "<$psiMysqlTsv";
while(my $line=<FF>){
	#asId, experimentId, tissue, treatment, jcecIncl, jcecSkip, jcecInclFormLen, jcecSkipFormLen, jcecPsi, jcecPsiValid, jcIncl, jcSkip, jcFormLen, jcSkipFormLen, jcPsi, jcPsiValid___ATHAA3SS0000003578, SRX4020057, leaf, NA, 0, 74, 138, 125, 0, Y, 0, 74, 137, 125, 0, Y
	chomp($line);
	($nameList, $valueList) = split(/___/, $line);
	@name = split(/, /, $nameList);
	@value = split(/, /, $valueList);
	for(my $i=0; $i<=$#name; $i++){
		$tmpPsi{$name[$i]} = $value[$i];
	}
	next if(uc($filterValid) eq "Y" and $tmpPsi{"jcecPsiValid"} ne "Y");

	# 登记specifiedField入hash
	$totalSpecifiedField{$tmpPsi{$specifiedField}} = 1;

	# 登记psi值
	$psiHref->{$tmpPsi{"asId"}}->{"psiList"} .= $tmpPsi{"jcecPsi"} . ",";
	
	# 检测指定字段值是否存在
	if(index($psiHref->{$tmpPsi{"asId"}}->{$specifiedField}, $tmpPsi{$specifiedField})<0){
		$psiHref->{$tmpPsi{"asId"}}->{$specifiedField} = $psiHref->{$tmpPsi{"asId"}}->{$specifiedField} . $tmpPsi{$specifiedField} . ",";
	}
}
close FF;

# 根据指定Field的总数，生成各个组开始、结束范围以及标签
my @specifiedField = keys(%totalSpecifiedField);
my $specifiedFieldNum = $#specifiedField+1;
my (%groupHash, $groupHashHref, $finalGroupSize, $groupBegin, $groupEnd, $groupName);
$groupHashHref=\%groupHash;
if($groupPctgNum eq "Pctg"){
	$finalGroupSize = int($groupSize/100*$specifiedFieldNum) + 1;
}else{
	$finalGroupSize = $groupSize;
}

#print "specifiedFieldNum=" . $specifiedFieldNum . "\n";
#print "finalGroupSize:" . $finalGroupSize . "\n";
# 开始生成groupHash
for(my $i=0; $i<=10000; $i++){
	# 计算group的开始和结束
	$groupBegin = $i * $finalGroupSize + 1;
	$groupEnd = ($i+1) * $finalGroupSize;

	# 生成group的name
	if($groupPctgNum ne "Pctg"){
		$groupName = $groupBegin . "-" . $groupEnd;
	}else{
		my $start = $i * $groupSize + 1;
		my $stop = ($i+1) * $groupSize;
		$groupName = $start . "-" . $stop;
	}

	# 登记group信息
	$groupHashHref->{$groupName}->{"start"} = $groupBegin;
	$groupHashHref->{$groupName}->{"stop"} = $groupEnd;
	$groupHashHref->{$groupName}->{"count"} = 0;

	print join("\t", $groupName, $groupHashHref->{$groupName}->{"start"}, $groupHashHref->{$groupName}->{"stop"}) . "\n";

	last if($groupEnd >= $specifiedFieldNum);
}


open WW, ">$outputPsiTsvForViolin";
print WW join("\t", "group", "psi") . "\n";
# 重新扫描AS，计算每个AS的specifiedField的数量
my $specifiedFieldList="";
my @specifiedField = keys(%psi);
my @asId = keys(%psi);
foreach my $asId(@asId){

	# 获得当前AS在多少specifiedField中出现
	@specifiedField = ();
	$specifiedFieldList = substr($psiHref->{$asId}->{$specifiedField}, 0, length($psiHref->{$asId}->{$specifiedField}) -1);
	@specifiedField = split(/,/, $specifiedFieldList);
	my $num = $#specifiedField + 1;

	# 找出当前AS对应的groupName
	my $groupName = &getGroupName($groupHashHref, $num);
	$psiHref->{$asId}->{"groupName"} = $groupName;
	$groupHashHref->{$groupName}->{"count"}++;
	
	# 输出所有的psi值，格式：
	# groupName PSI
	# 1-5  0.345
	# 1-5  0.123
	# 6-10 0.987
	# 6-10 0.901
	my $psiList = substr($psiHref->{$asId}->{"psiList"}, 0, length($psiHref->{$asId}->{"psiList"})-1);
	my @psi = ();
	@psi = split(/,/, $psiList);
	foreach my $psi(@psi){
		print WW join("\t", $psiHref->{$asId}->{"groupName"}, $psi) . "\n";
	}
}
close WW;

# 输出每个group中的AS的数量
open WW, ">$outputAsNumInGroup";
print WW join("\t", "group", "asNum") . "\n";
my @groupName = keys(%groupHash);
foreach $groupName(@groupName){
	print WW join("\t", $groupName, $groupHashHref->{$groupName}->{"count"}) . "\n";
}
close WW;

# 计算groupName
sub getGroupName{
	my ($href, $num) = @_;
	my (@groupName, $groupName);
	@groupName = keys(%{$href});
	foreach $groupName(@groupName){
		if($num<=$href->{$groupName}->{"stop"} and $num>=$href->{$groupName}->{"start"}){
			return $groupName;
		}
	}
}
