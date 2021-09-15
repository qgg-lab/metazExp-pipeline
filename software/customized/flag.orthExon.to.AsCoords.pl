#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputOriginExonCoord orth.exon.coord.in.self.genome.tsv \\\n" .
                "--asCatalogFile A5SS.catalog \\\n" .
                "--asType A5SS \\\n" .
		"--outputAsCatalogWithOrthTag A5SS.catalog.with.orthTag \n";
	exit;
}

my ($inputOriginExonCoord, $asCatalogFile, $asType, $outputAsCatalogWithOrthTag);
GetOptions(
	'inputOriginExonCoord=s'=>\$inputOriginExonCoord, 
	'asCatalogFile=s'=>\$asCatalogFile, 
	'asType=s'=>\$asType, 
	'outputAsCatalogWithOrthTag=s'=>\$outputAsCatalogWithOrthTag
);

my ($line, @fields, $field);
my (%chrToOrthCoord, @chr, $chr, @orthId, $orthId);
open FF, "<$inputOriginExonCoord";
#chr28   2372812 2372887 Exon116411      0       +       ORTH0000000002
while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	${$chrToOrthCoord{$fields[0]}}{"chainList"} .= $fields[5] . "#";
	${$chrToOrthCoord{$fields[0]}}{"startList"} .= $fields[1] . "#";
	${$chrToOrthCoord{$fields[0]}}{"stopList"} .= $fields[2] . "#";
	${$chrToOrthCoord{$fields[0]}}{"orthIdList"} .= $fields[6] . "#";
}
close FF;

print "finish load orth exon\n";
# read AS
my (@fieldTitle, $i,  %as, $orthId);
open WW, ">$outputAsCatalogWithOrthTag";
open FF, "<$asCatalogFile";
$line=<FF>;
print WW $line;

chomp($line);
@fieldTitle = ();
@fieldTitle = split(/\t/, $line);
#print "@fieldTitle";
#<STDIN>;
while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	%as = ();
	
	for($i=0; $i<=$#fields; $i++){
		$as{$fieldTitle[$i]}=$fields[$i];
	}		
	
	print WW join("\t", $fields[0], $fields[1], $fields[2], $fields[3], $fields[4]); 

	if(uc($asType) eq "A5SS"){
		#ASID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"longExonStart_0base"} + 1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"longExonEnd"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"shortES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"shortEE"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"flankingES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"flankingEE"}) . "\n";
	}elsif(uc($asType) eq "A3SS"){
		#ASID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"longExonStart_0base"} + 1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"longExonEnd"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"shortES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"shortEE"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"flankingES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"flankingEE"}) . "\n";
	}elsif(uc($asType) eq "SE"){
		#ASID GeneID geneSymbol chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"exonStart_0base"} + 1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"exonEnd"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"upstreamES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"upstreamEE"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"downstreamES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"downstreamEE"}) . "\n";

	}elsif(uc($asType) eq "RI"){
		#ASID GeneID geneSymbol chr strand riExonStart_0base riExonEnd upstreamES upstreamEE downstreamES downstreamEE
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"riExonStart_0base"} + 1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"riExonEnd"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"upstreamES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"upstreamEE"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"downstreamES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"downstreamEE"}) . "\n";

	}elsif(uc($asType) eq "MXE"){
		#ASID GeneID geneSymbol chr strand 1stExonStart_0base 1stExonEnd 2ndExonStart_0base 2ndExonEnd upstreamES upstreamEE downstreamES downstreamEE
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"1stExonStart_0base"} + 1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"1stExonEnd"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"2ndExonStart_0base"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"2ndExonEnd"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"upstreamES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"upstreamEE"});
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"downstreamES"}+1);
		print WW "\t" . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"downstreamEE"}) . "\n";
	
	}
}
close FF;

sub getOrthId{
	my ($hs, $chr, $chain, $pos) = @_;
	my (@chain, @start, @stop, @orthId, $i);
	
	@chain = split(/#/, ${$hs->{$chr}}{"chainList"});
	@start = split(/#/, ${$hs->{$chr}}{"startList"});
	@stop = split(/#/, ${$hs->{$chr}}{"stopList"});
	@orthId = split(/#/, ${$hs->{$chr}}{"orthIdList"});
	for($i=0; $i<=$#chain; $i++){
		if($chain[$i] eq $chain and $pos <= $stop[$i] and $pos >= $start[$i]){
			return $orthId[$i];
		}
	}
	return "-";
}
