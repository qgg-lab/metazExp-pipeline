#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--psiTsvFile psi.tsv\\\n" .
                "--experimentTsvFile experiment.tsv \\\n" .
                "--minMappedSpotNum 20 \\\n" .
                "--minReadLength 40 \\\n" .
		"--libraryList RF,FR,UN,R,F,U \\\n" .
		"--layoutList SINGLE,PAIRED \\\n" .
                "--tissueExpNum 6 \\\n" .
                "--inclusionReadNum 1 \\\n" .		
                "--skippingReadNum 1 \\\n" .
		"--outputAsTissuePsiFile final.selected.as.with.tissue.psi.tsv \n";
	exit;
}

my ($psiTsvFile, $experimentTsvFile, $minMappedSpotNum, $minReadLength, $tissueExpNum, $inclusionReadNum, $skippingReadNum, $outputAsTissuePsiFile, $layoutList, $libraryList);

GetOptions(
        'psiTsvFile=s'=>\$psiTsvFile,
        'experimentTsvFile=s'=>\$experimentTsvFile,
        'minMappedSpotNum=s'=>\$minMappedSpotNum,
        'minReadLength=s'=>\$minReadLength,
	'librarylist=s'=>\$libraryList,
	'layoutList=s'=>\$layoutList,
        'tissueExpNum=s'=>\$tissueExpNum,
        'inclusionReadNum=s'=>\$inclusionReadNum,
        'skippingReadNum=s'=>\$skippingReadNum,
        'outputAsTissuePsiFile=s'=>\$outputAsTissuePsiFile,
);

my (%expInfo, $expId, @tt, $line, @field);
my (@field);
my ($fieldNameString, $valueString);

# 读取experiment的相关信息，包括组织类型、比对情况等等
# 将这些experiment的相关信息读入到哈希中，便于访问
open FF, "<$experimentTsvFile";
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);
        @field = split(/, /, $valueString);
        $expId = substr($field[2], 1, length($field[2]) - 2);
	${$expInfo{$expId}}{"mappedSpotNum"} = $field[1];
	${$expInfo{$expId}}{"readLength"} = substr($field[6], 1, length($field[6]) - 2);
	@tt = split(/,/, ${$expInfo{$expId}}{"readLength"});
	${$expInfo{$expId}}{"readLength"} = $tt[0];
	${$expInfo{$expId}}{"tissue"} = substr($field[$#field], 1, length($field[$#field]) - 2);
	${$expInfo{$expId}}{"library"} = substr($field[4], 1, length($field[4]) - 2);
	${$expInfo{$expId}}{"layout"} = substr($field[5], 1, length($field[5]) - 2);
}
close FF;
close WW;

# 从总psi表格中读取每个pis值及其对应的experiment的编号（同时可以从上面哈希中获得experment的相关信息）
# 这里不该用inclusion和skipping read num对psi值执行了相关过滤。这里传递过来的inclusion和skipping参数为0，因此相当于不再过滤
# 我们在生成AS目录时，其实已经用inclusion和skipping 的支持read数量进行过滤了。
my ($asId, $inReadNum, $skipReadNum, $psi, %asInfo, @asId);
open FF, "<$psiTsvFile";
while($line=<FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);
	@field = split(/, /, $valueString);
	$asId = substr($field[0], 1, length($field[0]) - 2);
	$inReadNum = $field[4];
	$skipReadNum = $field[5];
	$psi = $field[8];
	$expId = substr($field[2], 1, length($field[2]) - 2);

	next if(${$expInfo{$expId}}{"mappedSpotNum"} < $minMappedSpotNum);
	next if(${$expInfo{$expId}}{"readLength"} < $minReadLength);
#	next if($inReadNum <= $inclusionReadNum);
#	next if($skipReadNum <= $skippingReadNum);
	next if($inReadNum <= $inclusionReadNum and $skipReadNum <= $skippingReadNum);
	next if(&checkLibrary(${$expInfo{$expId}}{"library"}, $libraryList) == 0);
	next if(&checkLayout(${$expInfo{$expId}}{"layout"}, $layoutList) == 0);

	${$asInfo{$asId}}{${$expInfo{$expId}}{"tissue"}} .= $psi . ",";
}
close FF;

print "Finish load psi into hash\n";

my (@tissue, @psi, $tissue);
@asId = keys(%asInfo);
open WW, ">$outputAsTissuePsiFile";
foreach $asId(@asId){
	print WW $asId;
	@tissue = ();
	@tissue = keys(%{$asInfo{$asId}});
	foreach $tissue(@tissue){
		@psi = ();
		@psi = split(/,/, ${$asInfo{$asId}}{$tissue});
		if($#psi >= $tissueExpNum - 1){
			print WW "\t" . $tissue;
			foreach $psi(@psi){
				print WW "," . $psi;
			}
		}
	}
	print WW "\n";
}
close WW;


sub checkLibrary{
	my ($expLibraryList, $conditionLibraryList) = @_;
	my (@expLibrary, $expLibrary, @conditionLibrary, $conditionLibrary);

	@expLibrary = split(/,/, $expLibraryList);
	@conditionLibrary = split(/,/, $conditionLibraryList);
	foreach $expLibrary(@expLibrary){
		foreach $conditionLibrary(@conditionLibrary){
			return 1 if($expLibrary eq $conditionLibrary);
		}
	}
	return 0;
}

sub checkLayout{
	my ($expLayoutList, $conditionLayoutList) = @_;
	my (@expLayout, $expLayout, @conditionLayout, $conditionLayout);

	@expLayout = split(/,/, $expLayoutList);
	@conditionLayout = split(/,/, $conditionLayoutList);

	foreach $expLayout(@expLayout){
		foreach $conditionLayout(@conditionLayout){
			return 1 if($expLayout eq $conditionLayout);
		}
	}
	return 0;
}

