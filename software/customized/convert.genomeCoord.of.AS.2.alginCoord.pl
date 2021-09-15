#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
                "--asFile SE.catalog \\\n" .
                "--genomeAlignCoordFile total.genome.align.coord.tsv \\\n" .
                "--outputAsFile  SE.catalog.align.coord\n";
        exit;
}

my ($asFile, $genomeAlignCoordFile, $outputAsFile);

GetOptions(
        'asFile=s'=>\$asFile,
        'genomeAlignCoordFile=s'=>\$genomeAlignCoordFile,
        'outputAsFile=s'=>\$outputAsFile,
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

my ($chr, $i, $chain);
open WW, ">$outputAsFile";
open FF, "<$asFile";
#ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
#GGALA3SS0000004913      "ENSGALG00000000003"    "PANX2" chr1    +       20336702        20336832        20336707        20336832        20329069        20329191
while($line=<FF>){
	if($line=~/^ASID/){
		print WW $line;
		next;
	}
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);

	$chr = $fields[3];
	$chr = $1 if($fields[3]=~/chr0*(\d+)/);
	$chain = $fields[4];
	
	print WW $fields[0];
	for($i=1; $i<=4; $i++){
		print WW "\t" . $fields[$i];
	}
	for($i=5; $i<=$#fields; $i++){
		print WW "\t" . &obtainAlignCoord(\%genomeToAlignCoord, \@blockId, $chr, $fields[$i], $chain);
	}
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

