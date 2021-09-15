#!/usr/bin/perl
use strict;
if($#ARGV<0){
	print "\nperl $0 final.complete.trspt.pep.fa \n\n\n";
	exit;
}
my $inputPepSeqFile = $ARGV[0];
my (%pepSeq, $id, $line, @id);
open FF, "<$inputPepSeqFile";
while(my $line =<FF>){
        chomp($line);
        if($line=~/>(.*)/){
                $id = $1;
        }else{
                $pepSeq{$id}.=$line;
        }
}
close FF;

open WW, ">$inputPepSeqFile";
@id=keys(%pepSeq);
foreach $id(@id){
	$pepSeq{$id}=~s/\*/X/g;
	print WW ">" . $id . "\n";
	print WW $pepSeq{$id} . "\n";
}
close WW;
