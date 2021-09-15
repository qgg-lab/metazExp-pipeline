#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--species \"Equus caballus\" \\\n" .
		"--genomeFile  genome.fa \\\n" . 
		"--inputAsTsvFile 89462.as.tsv\\\n" .
		"--inputSpliceFile splice.tsv \\\n" .
		"--inputVcfFile genome.vcf \\\n" .
		"--variantTsvFile 89462.variant.Tsv \n\n";
	exit;
}

my ($genomeFile, $inputVcfFile, $outputVcfFile, $species, $variantTsvFile, $inputAsTsvFile, $inputSpliceFile);
my ($dbName, $dbUser, $dbPWD);
GetOptions(
        'genomeFile=s'=>\$genomeFile,
	'species=s'=>\$species,
        'inputVcfFile=s'=>\$inputVcfFile,
	'inputAsTsvFile=s'=>\$inputAsTsvFile,
	'inputSpliceFile=s'=>\$inputSpliceFile,
	'variantTsvFile=s'=>\$variantTsvFile,
);

# read genome sequence into hash
my (%genomeSeq, $line, @tt, $id);
open FF, "<$genomeFile";
while($line = <FF>){
        chomp($line);
        if($line=~/>/){
                @tt = ();
                @tt = split(/ /, $line);
                if($tt[0]=~/>ref\|(.*?)\|/ or $tt[0]=~/>gi.*\|ref\|(.*?)\|/ or  $tt[0]=~/>(.*)/){
                        $id = $1;
                }
        }else{
                $genomeSeq{$id} .=uc($line);
        }
}
close FF;

print "Finish load genome sequence into hash.\n";

my ($line, $sblen, $i, %as);
my (@tmp, $fieldNameString, @fieldName, $fieldName, $valueString, @value, $value);
open FF, "<$inputAsTsvFile";
while($line=<FF>){
	chomp($line);
        @tmp = ();
        @tmp = split(/_____/, $line);
        $fieldNameString = $tmp[0];
        @fieldName = ();
        @fieldName = split(/, /, $fieldNameString);

	$valueString = $tmp[1];
        @value = ();
        @value = split(/, /, $valueString);

        for($i=0; $i<=$#value; $i++){
                if($value[$i]=~/"(.*)"/){
                        $value[$i] = $1;
                }
                $as{$fieldName[$i]}=$value[$i];
        }

	$sblen = $as{"end"} - $as{"start"} + 1;
	substr($genomeSeq{$as{"chr"}}, $as{"start"}-1, $sblen) = "0"x$sblen;
}
close FF;

# 读取内含子，提取donor和acceptor位点
my ($donor1, $donor2, $acceptor1, $acceptor2, $chr, $base);
open FF, "<$inputSpliceFile";
#        base0
#1       339311  342546  -
#1       342720  346601  -
#1       346888  346924  -
#1       346923  350266  -
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);
	$chr = $tmp[0];
	if($tmp[3] eq "+"){
		$donor1 = $tmp[1] + 1 + 1;
		$donor2 = $tmp[1] + 1 + 2;
		$acceptor1 = $tmp[2] + 1 - 1;
		$acceptor2 = $tmp[2] + 1 - 1 - 1;
	}else{
		$acceptor1 = $tmp[1] + 1 + 1;
		$acceptor2 = $tmp[1] + 1 + 2;
		$donor1 = $tmp[2];
		$donor2 = $tmp[2] - 1;
	}
	# acceptor site
	$base = substr($genomeSeq{$chr}, $acceptor1-1, 1);
	if($base eq "0" or $base eq "1"){
		substr($genomeSeq{$chr}, $acceptor1-1, 1)+=2;
	}
	$base = substr($genomeSeq{$chr}, $acceptor2, 1);
	if($base eq "0" or $base eq "1"){
		substr($genomeSeq{$chr}, $acceptor2, 1) +=2;
	}
	# donor site
	$base = substr($genomeSeq{$chr}, $donor1-1, 1);
	if($base eq "0" or $base eq "2"){
		substr($genomeSeq{$chr}, $donor1-1, 1)+=1;
	}
	$base = substr($genomeSeq{$chr}, $donor2-1, 1);
	if($base eq "0" or $base eq "2"){
		substr($genomeSeq{$chr}, $donor2-1, 1) +=1;
	}

}
close FF;

print "Finish mask as region.\n";

# filter vcf position
open WW, ">$variantTsvFile";
my ($chr, $pos, $ref, $alt, $accId, $varType);
open FF, "<" . $inputVcfFile;
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#1       207381  rs719849645     C       T       .       .       dbSNP_150;TSA=SNV
while($line=<FF>){
	next if($line=~/#/);
	@tt = ();
	@tt = split(/\t/, $line);
	($chr, $pos, $accId, $ref, $alt, $varType) = ($tt[0], $tt[1], $tt[2], $tt[3], $tt[4], $tt[7]);
	if($varType=~/=SNV/){
		$varType="SNV";
	}elsif($varType=~/=deletion/){
		$varType="deletion";
	}elsif($varType=~/=insertion/){
		$varType="insertion";
	}else{
		next;	
	}

	if(substr($genomeSeq{$tt[0]}, $tt[1]-1, 1) eq "0" or substr($genomeSeq{$tt[0]}, $tt[1]-1, 1) eq "1" or substr($genomeSeq{$tt[0]}, $tt[1]-1, 1) eq "2" or substr($genomeSeq{$tt[0]}, $tt[1]-1, 1) eq "3"){

		print WW "species, chr, pos, ref, alt, accId, varType, hitElement" . "_____" . "\"" . $species . "\", \"" . $chr . "\", $pos, \"$ref\", \"$alt\", \"$accId\", \"$varType\", " . substr($genomeSeq{$tt[0]}, $tt[1]-1, 1) . "\n";
	}
}
close FF;
close WW;
