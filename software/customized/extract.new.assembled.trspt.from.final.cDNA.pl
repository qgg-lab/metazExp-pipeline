#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--finalCdnaFasta final.cdna.fa \\\n" .
                "--finalGtf finalGtf.gtf \\\n" .
		"--methodList StringTie \\\n" .
                "--cDNAOfNewAssembledTrsptFasta cDNA.of.new.assembled.trspt.fa \n";
	exit;
}

my ($finalCdnaFasta, $finalGtf, $methodList, $newassembledtrsptcdnafasta);

GetOptions(
        'finalCdnaFasta=s'=>\$finalCdnaFasta,
        'finalGtf=s'=>\$finalGtf,
	'methodList=s'=>\$methodList,
        'newassembledtrsptcdnafasta=s'=>\$newassembledtrsptcdnafasta,
);

# 将所有转录本序列读入到hash
my (%allCdna, @tmp, $line, $id, @id);
open FF, "<$finalCdnaFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$id = $1;
		@tmp = ();
		@tmp = split(/ /, $id);
		$id = $tmp[0];
		#print $id;
		#<STDIN>;
	}else{
		$allCdna{$id}.=$line;
	}
}
close FF;

# 将标注为StringTie的转录本挑选出来放在hash中
my (%newAssembledTrspt, @field, $trsptId);
open FF, "<$finalGtf";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "transcript" or index($methodList, $field[1])<0);
	$trsptId = &getTrsptIdInAttr($field[8]);
	$newAssembledTrspt{$trsptId}=1;
#	print $trsptId;
#	<STDIN>;
}
close FF;

# 输出新组装的转录本
my (@trsptId);
@trsptId=keys(%newAssembledTrspt);
open WW, ">$newassembledtrsptcdnafasta";
foreach $trsptId(@trsptId){
	print WW ">$trsptId\n";
	print WW $allCdna{$trsptId} . "\n";
}
close WW;

sub getTrsptIdInAttr{
	my ($attrString) = @_;
	my (@attr, $attr);
	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/transcript_id "(.*)"/){
			return $1;
		}
	}
}
