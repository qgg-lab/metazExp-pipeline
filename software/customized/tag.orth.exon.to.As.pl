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
print WW "orthIdList\t";

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
#	print $line;
	@fields = ();
	@fields = split(/\t/, $line);
	%as = ();
	
	for($i=0; $i<=$#fields; $i++){
		$as{$fieldTitle[$i]}=$fields[$i];
	}		
#	my @ff = keys(%as);
#	print "@ff";
#<STDIN>;
	if(uc($asType) eq "A5SS"){

		if($as{"chr"} eq "+"){
			$orthId = &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"longExonEnd"});
			$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"shortEE"});
		}else{
			$orthId = &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"longExonStart_0base"} + 1);
			$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"shortES"} + 1);
		}

	}elsif(uc($asType) eq "A3SS"){

		if($as{"chr"} eq "+"){
			$orthId = &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"longExonStart_0base"}+1);
			$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"shortES"}+1);
		}else{
			$orthId = &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"longExonEnd"});
			$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"shortEE"});
		}

	}elsif(uc($asType) eq "SE"){

		$orthId = &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"exonStart_0base"}+1);
		$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"exonEnd"});				

	}elsif(uc($asType) eq "RI"){

		$orthId = &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"upstreamEE"});
		$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"downstreamES"}+1);

	}elsif(uc($asType) eq "MXE"){

		$orthId = &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"1stExonStart_0base"}+1);
		$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"1stExonEnd"});
		$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"2ndExonStart_0base"}+1);
		$orthId .= "," . &getOrthId(\%chrToOrthCoord, $as{"chr"}, $as{"strand"}, $as{"2ndExonEnd"});

	}
	print WW $orthId . "\t" . $line . "\n";
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
