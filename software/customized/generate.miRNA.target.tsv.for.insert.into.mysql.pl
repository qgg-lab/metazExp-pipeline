#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--targetTsv \\\n" .
		"--outputMysqlTsv \n";
	exit;
}

my ($targetTsv, $outputMysqlTsv);

GetOptions(
        'targetTsv=s'=>\$targetTsv,
        'outputMysqlTsv=s'=>\$outputMysqlTsv,
);

# miRNA_Acc.      Target_Acc.     Expectation     UPE$    miRNA_start     miRNA_end       Target_start    Target_end      miRNA_aligned_fragment alignment       Target_aligned_fragment Inhibition      Target_Desc.    Multiplicity
# ath-miR156i     SRX399568.3733.4        0.0     -1.0    1       20      643     662     UGACAGAAGAGAGAGAGCAG    ::::::::::::::::::::    CUGCUCUCUCUCUUCUGUCA    Cleavage        transcript_name:NA gene_id:AT1G53160 gene_name:NA       1 
my ($fieldString, $valueString, $line, @nameField, @valueField, $i);
open WW, ">$outputMysqlTsv";
open FF, "<$targetTsv";
<FF>;
$line=<FF>;
chomp($line);
@nameField = split(/\t/, $line);
for(my $i=0; $i<=$#nameField; $i++){
	$nameField[$i]=~s/\.//g;
}
$nameField[3]=substr($nameField[3], 0, length($nameField[3])-1);
while($line=<FF>){
	chomp($line);
	@valueField=split(/\t/, $line);

	# fieldString
	$fieldString = join(", ", @nameField);

	# valueString
	$valueString = join(", ", @valueField);

	print WW join("___", $fieldString, $valueString) . "\n";
}
close FF;
close WW;
