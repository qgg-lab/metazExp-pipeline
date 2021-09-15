#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;

if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--process \\\n" .
                "--cDNAfasta \\\n" .
                "--altedCdnaFasta \n";
	exit;
}

my ($process, $cDNAfasta, $altedCdnaFasta);

GetOptions(
        'process=s'=>\$process,
        'cDNAfasta=s'=>\$cDNAfasta,
        'altedCdnaFasta=s'=>\$altedCdnaFasta,
);

# 将原始cDNA序列读入到hash
my (%origTrsptSeq, $id, @tt, $line);
open FF, "<$cDNAfasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$id = $1;
		@tt = ();
		@tt = split(/ /, $id);
		$id = $tt[0];
	}else{
		$origTrsptSeq{$id}.=$line;
	}
}
close FF;


# 读取process文件对原来转录本序列进行操作改变，然后输出
my (@titleField, $titleField, @valueField, $valueField, $i, %process, $processHref);
my ($preSeq, $tailSeq, $cutSeq, $insertSeq, $cutBegin, $cutEnd, $tailBegin);
$processHref = \%process;
open WW, ">$altedCdnaFasta";
open FF, "<$process";
# chr     strand  asId    asType  trsptId residentType    editSite        cutSize insertSize      insertSeq
$line = <FF>;
chomp($line);
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line);
	%process = ();
	for($i=0; $i<=$#valueField; $i++){
		$processHref->{$titleField[$i]}=$valueField[$i];
	}
	
	# |--- preSeq ----||----- cutSeq -----|| ------------ tailSeq ---------|
	# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	# 计算editSite之前的序列（包括编辑位点在内）
#	if($processHref->{"trsptId"} eq "AT3G15390.1" and $processHref->{"asId"} eq "ATHASE0000035094"){
#		print $line . "\n";
#		print join("\t", "editSite:", $processHref->{"editSite"});
#		<STDIN>;
#	}
	
	$preSeq = substr($origTrsptSeq{$processHref->{"trsptId"}}, 0, $processHref->{"editSite"});
	# 计算被切掉的序列
	# based on 1
	$cutBegin = $processHref->{"editSite"} + 1;
	# based on 0
	$cutBegin = $cutBegin - 1;
	$cutSeq = substr($origTrsptSeq{$processHref->{"trsptId"}}, $cutBegin, $processHref->{"cutSize"});

	# 计算tailSeq(剩余)
	# based on 1 
	$tailBegin = $processHref->{"editSite"} + $processHref->{"cutSize"} + 1;
	# based on 0
	$tailBegin = $tailBegin - 1;
	$tailSeq = substr($origTrsptSeq{$processHref->{"trsptId"}}, $tailBegin);

#	if($processHref->{"trsptId"} eq "AT3G15390.1" and $processHref->{"asId"} eq "ATHASE0000035094"){
#		print join("\t", "preSeqLen:", length($preSeq), "cutSeqLen:", length($cutSeq), "InsertSeqLen:", length($processHref->{"insertSeq"}), "tailSeqLen:", length($tailSeq)) . "\n";
#		<STDIN>;
#	}
	print WW ">" . $processHref->{"trsptId"} . "XXX" . $processHref->{"asId"} . "\n";
	print WW $preSeq . $processHref->{"insertSeq"} . $tailSeq . "\n";
}
close FF;
