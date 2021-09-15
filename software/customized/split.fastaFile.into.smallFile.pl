#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--largeFastaFile \\\n" .
                "--seqNum 5000 \\\n" .
		"--smallFastaFilePrefix ./4081\n";
	exit;
}

my ($largeFastaFile, $seqNum, $smallFastaFilePrefix);

GetOptions(
        'largeFastaFile=s'=>\$largeFastaFile,
        'seqNum=s'=>\$seqNum,
        'smallFastaFilePrefix=s'=>\$smallFastaFilePrefix,
);

my (%seq, @id, $id, $line);
open FF, "<$largeFastaFile";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		@id = split(/ / , $1);
		$id = $id[0];
	}else{
		$seq{$id}.=$line;
	}
}
close FF;

my ($smallFileNum);
@id = ();
@id = keys(%seq);
if(int(($#id + 1)/$seqNum)*$seqNum == $#id + 1){
	$smallFileNum = ($#id + 1)/$seqNum;
}else{
	$smallFileNum = int(($#id + 1)/$seqNum) + 1;
}

for(my $i=0; $i<$smallFileNum; $i++){
	open WW, ">$smallFastaFilePrefix" . "_" . $i . ".fa";
	for(my $j=$i*$seqNum; $j<($i+1)*$seqNum and $j<=$#id; $j++){
		print WW ">" . $id[$j] . "\n";
		print WW $seq{$id[$j]} . "\n";
	}
	close WW;
}
