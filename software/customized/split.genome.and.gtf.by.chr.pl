#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--fastaFile \\\n" .
                "--gtfFile \\\n" .
                "--outputDir \n";
	exit;
}

my ($fastaFile, $gtfFile, $outputDir);

GetOptions(
        'fastaFile=s'=>\$fastaFile,
        'gtfFile=s'=>\$gtfFile,
        'outputDir=s'=>\$outputDir,
);

my (%seqAndGtf, $seqAndGtfHref, $line, @seqId, $seqId, @field);
$seqAndGtfHref=\%seqAndGtf;
# 将sequence读入hash
open FF, "<$fastaFile";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@field = ();
		@field = split(/ /, $1);
		$seqId = $field[0];
	}else{
		$seqAndGtfHref->{$seqId}->{"seq"}.=$line;
	}
}
close FF;

# 将gtf读入
open FF, "<$gtfFile";
while($line=<FF>){
	@field = ();
	@field = split(/\t/, $line);
	next if($#field !=8);
	$seqAndGtfHref->{$seqId}->{"gtf"}.=$line;
}
close FF;

# 按照chr编号和other输出文件
@seqId = ();
@seqId = keys(%seqAndGtf);
system("mkdir -p $outputDir");
open NONCHRSEQ, ">$outputDir" . "/nonChr.fa";
open NONCHRGTF, ">$outputDir" . "/nonChr.gtf";
foreach $seqId(@seqId){
	if($seqId=~/\d+/){
		open CHRSEQ, ">$outputDir" . "/" . $seqId . ".fa";
		print CHRSEQ ">" . $seqId . "\n";
		print CHRSEQ $seqAndGtfHref->{$seqId}->{"seq"} . "\n";
		close CHRSEQ;

		open CHRGTF, ">$outputDir" . "/" . $seqId . ".gtf";
		print CHRGTF $seqAndGtfHref->{$seqId}->{"gtf"};
		close CHRGTF;
	}else{
		print NONCHRSEQ ">" . $seqId . "\n";
		print NONCHRSEQ $seqAndGtfHref->{$seqId}->{"seq"};

		print NONCHRGTF $seqAndGtfHref->{$seqId}->{"gtf"};
	}
}
close NONCHRSEQ;
close NONCHRGTF;
