#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "(1)将pep序列的Id转换成对应的cDNA序列编号\n";
	print "(2)将多样序列转换成单行序列\n";
	print "\nperl $0 \\\n" . 
                "--gtfFile \\\n" .
                "--oldCdnaSeqFile \\\n" .
                "--oldPepSeqFile \\\n" .
		"--newCdnaSeqFile \\\n" .
		"--newPepSeqFile \n";
	exit;
}

my ($gtfFile, $oldCdnaSeqFile, $oldPepSeqFile, $newPepSeqFile, $newCdnaSeqFile);

GetOptions(
        'gtfFile=s'=>\$gtfFile,
        'oldCdnaSeqFile=s'=>\$oldCdnaSeqFile,
        'oldPepSeqFile=s'=>\$oldPepSeqFile,
        'newPepSeqFile=s'=>\$newPepSeqFile,
	'newCdnaSeqFile=s'=>\$newCdnaSeqFile,
);

my (%pepSeq, %trsptSeq);
my (%pepIdToTrsptId, $trsptId, $pepId);
my ($line, @field);
open FF, "<$gtfFile";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "CDS");
	($trsptId, $pepId) = ("", "");
	&getId($field[8], \$trsptId, \$pepId);
	if($trsptId ne ""){
		$pepIdToTrsptId{$pepId} = $trsptId;
	}
}
close FF;

my ($line, $id, @tt);
# 将trspt序列读入
open FF, "<$oldCdnaSeqFile";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$id = $1;
		@tt = ();
		@tt = split(/ /, $id);
		$id = $tt[0];
	}else{
		$trsptSeq{$id} .= $line;
	}
}
close FF;

# 将pep序列读入
open FF, "<$oldPepSeqFile";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$id = $1;
		@tt = ();
		@tt = split(/ /, $id);
		$id = $tt[0];
	}else{
		$pepSeq{$id} .= $line;
	}
}
close FF;

# 将新trspt序列输出
my (@id);
@id = keys(%trsptSeq);
open WW, ">$newCdnaSeqFile";
foreach $id(@id){
	print WW ">$id\n";
	print WW $trsptSeq{$id} . "\n";
}
close WW;

# 将新pep序列输出
@id = ();
@id = keys(%pepSeq);
open WW, ">$newPepSeqFile";
foreach $id(@id){
	if(exists($pepIdToTrsptId{$id})){
		print WW ">" . $pepIdToTrsptId{$id} . "\n";
	}else{
		print WW ">$id\n";
	}
	print WW $pepSeq{$id} . "\n";
}
close WW;

sub getId{
	my ($attr, $trsptId, $pepId) = @_;
	my (@tmp, $tmp);
	@tmp = split(/;/, $attr);
	foreach $tmp(@tmp){
		if($tmp=~/transcript_id "(.*?)"/){
			$$trsptId = $1;
		}
		if($tmp=~/protein_id "(.*?)"/){
			$$pepId = $1;
		}
	}
}
