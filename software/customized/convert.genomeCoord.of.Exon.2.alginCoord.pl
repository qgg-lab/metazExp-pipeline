#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
                "--exonFile 9031.exon.tsv \\\n" .
                "--genomeAlignCoordFile total.genome.align.coord.tsv \\\n" .
                "--outputExonFile  9031.exon.align.coord\n";
        exit;
}

my ($exonFile, $genomeAlignCoordFile, $outputExonFile);

GetOptions(
        'exonFile=s'=>\$exonFile,
        'genomeAlignCoordFile=s'=>\$genomeAlignCoordFile,
        'outputExonFile=s'=>\$outputExonFile,
);

my (%genomeToAlignCoord);
my ($line, @fields, $blockId, @blockId);
open FF, "<$genomeAlignCoordFile";
#gallus_gallus   2       +       9937535 9937535 15460000008917
#7366
#gallus_gallus   2       +       9937536 9988269 15460000008919
#467245
#467246
print "Begin load genome coord time:\n";
print &getTime() . "\n";
#print scalar(@{${$hs{"001"}}{"coord"}});
while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	if($#fields == 5){
		$blockId = $fields[5];
		${$genomeToAlignCoord{$blockId}}{"chr"} = $fields[1];
		${$genomeToAlignCoord{$blockId}}{"chain"} = $fields[2];
		${$genomeToAlignCoord{$blockId}}{"start"} = $fields[3];
		${$genomeToAlignCoord{$blockId}}{"stop"} = $fields[4];
		${$genomeToAlignCoord{$blockId}}{"coord"} = ();
	}else{
		${$genomeToAlignCoord{$blockId}}{"coord"}[scalar(@{${$genomeToAlignCoord{$blockId}}{"coord"}})] = $fields[0];	
	}
}
close FF;
print "Finish load genome align coord time:\n";
print &getTime() . "\n";

# obtain block Id into array
@blockId = keys(%genomeToAlignCoord);
#print "@blockId";
#<STDIN>;

my ($chr, $i, $chain, $start, $stop, $exonId);
open WW, ">$outputExonFile";
open FF, "<$exonFile";
#EXON1   1       5273    5524    -
#EXON2   1       5536    5952    -
#EXON3   1       6755    6829    -
while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);

	$exonId = $fields[0];
	$chr = $fields[1];
	$start = $fields[2];
	$stop = $fields[3];
	$chain = $fields[4];
	
	print WW $exonId;
	print WW "\t" . &obtainAlignCoord(\%genomeToAlignCoord, \@blockId, $chr, $start, $chain);
	print WW "\t" . &obtainAlignCoord(\%genomeToAlignCoord, \@blockId, $chr, $stop, $chain);
	print WW "\n";
#	<STDIN>;
}
close FF;
close WW;
print "Finish convert time:\n";
print &getTime() . "\n";

sub obtainAlignCoord{
	my ($block, $blockIdArr, $chr, $pos, $chain) = @_;
	my ($blockId, $findFlag, $hitBlockId);


	$findFlag = 0;
	foreach $blockId(@{$blockIdArr}){
		if( ${$$block{$blockId}}{"chr"} eq $chr and ${$$block{$blockId}}{"chain"} eq $chain and $pos >= ${$$block{$blockId}}{"start"} and $pos <= ${$$block{$blockId}}{"stop"}){
			$findFlag = 1;
			$hitBlockId = $blockId;
			last;
		}
	}

	if($findFlag == 0){
		return "-";
	}else{
		return $hitBlockId . "_" . ${$$block{$hitBlockId}}{"coord"}[$pos - ${$$block{$hitBlockId}}{"start"}];
	}
}

sub getTime(){
        my @months = qw( 01 02 03 04 05 06 07 08 09 10 11 12 );
        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
        return $months[$mon] . "/" . $mday . " " . $hour . ":" . $min . ":" . $sec;
}

