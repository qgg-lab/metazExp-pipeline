#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV<0){
	print "\tperl $0 \\\n" . 
		"\t\t--asFile total.A5SS.tsv \\\n" .
		"\t\t--asType A5SS \\\n" .
		"\t\t--outputAsFile total.A5SS.with.exon\n\n";
	print "This script is used to append concat exon coordinates for each AS.\n\n";
	exit;
}

my ($asFile, $asType, $outputAsFile);

GetOptions(
	'asFile=s'=>\$asFile,
	'asType=s'=>\$asType,
	'outputAsFile=s'=>\$outputAsFile
);

if(uc($asType) ne "A5SS" and uc($asType) ne "A3SS" and uc($asType) ne "RI" and uc($asType) ne "SE" and uc($asType) ne "MXE"){
	print "asType must be specified as A5SS, A3SS, RI, SE or MXE\n";
	exit;
}


my ($line, @field, @tmpArr, %asHash, @fieldName);
my ($chain, $chrId,  $i,  $trspt1, $trspt2, $exon1, $exon2, $exon3, $exon4);
open FF, "<$asFile";
$line=<FF>;
chomp($line);
open WW, ">$outputAsFile";
print WW $line . "\tfirstLocalConcatExons\tsecondLocalConcatExons\n";

# construct %asHash and register fieldName into @fieldName
@field = ();
@field = split(/\t/, $line);
for($i=0; $i<=$#field; $i++){
	$asHash{$field[$i]}=1;
	$fieldName[$#fieldName+1]=$field[$i];
}

# 16___43110298___43110580___-
# read each field into
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$chrId = $field[3];
	$chain = $field[4];
	
	# read each field into asHash according to fieldName
	for($i=0; $i<=$#fieldName; $i++){
		$asHash{$fieldName[$i]}=$field[$i];
	}
	
	# output AS position
	if(uc($asType) eq "A5SS"){ # outputs are same with each other when on +/- chain 

		# obtain trspt1: chr1____12345____45678___+!chr1____567890___890123___+
		$trspt1 = "";
		$exon1 = join("___", $chrId, $asHash{"longExonStart_0base"}+1, $asHash{"longExonEnd"}, $chain);
		$exon2 = join("___", $chrId, $asHash{"flankingES"}+1, $asHash{"flankingEE"}, $chain);
		$trspt1 = join("!", $exon1, $exon2);
		
		# obtain trspt2: 
		$trspt2 = "";
		$exon1 = join("___", $chrId, $asHash{"shortES"}+1, $asHash{"shortEE"}, $chain);
		$exon1 = join("___", $chrId, $asHash{"flankingES"}+1, $asHash{"flankingEE"}, $chain);
		$trspt2 = join("!", $exon1, $exon2);

		# output as with concat exon coordinates
		print WW join("\t", $line, $trspt1, $trspt2) . "\n";
						
	}
	if(uc($asType) eq "A3SS"){ # outputs are same with each other when on +/- chain

		# obtain trspt1
		$trspt1 = "";
		$exon1 = join("___", $chrId, $asHash{"flankingES"}+1, $asHash{"flankingEE"}, $chain);
		$exon2 = join("___", $chrId, $asHash{"longExonStart_0base"}+1, $asHash{"longExonEnd"}, $chain);
		$trspt1 = join("!", $exon1, $exon2);

		# obtain trspt2
		$trspt2 = "";
		$exon1 = join("___", $chrId, $asHash{"flankingES"}+1, $asHash{"flankingEE"}, $chain);
		$exon2 = join("___", $chrId, $asHash{"shortES"}+1, $asHash{"shortEE"}, $chain);
		$trspt2 = join("!", $exon1, $exon2);

		# output as with concat exon coordinates
		print WW join("\t", $line, $trspt1, $trspt2) . "\n";

	}

	if(uc($asType) eq "RI"){

		if($chain eq "+"){

			# obtain trspt1: upstreamES+1…upstreamEE|downstreamES+1 … downstreamEE
			$trspt1 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$trspt1 = join("!", $exon1, $exon2);

			# obtain trspt2: riExonStart_0base + 1 … riExonEnd
			$trspt2 = "";
			$exon1 = join("___", $chrId, $asHash{"riExonStart_0base"}+1, $asHash{"riExonEnd"}, $chain);
			$trspt2 = join("!", $exon1);

			# output as with concat exon coordinates
			print WW join("\t", $line, $trspt1, $trspt2) . "\n";

		}elsif($chain eq "-"){
			
			# obtain trspt1:downstreamES+1 … downstreamEE | upstreamES+1…upstreamEE
			$trspt1 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$trspt1 = join("!", $exon1, $exon2);

			# obtain trspt2: riExonStart_0base … riExonEnd
			$trspt2 = "";
			$exon1 = join("___", $chrId, $asHash{"riExonStart_0base"}+1, $asHash{"riExonEnd"}, $chain);
			$trspt2 = $exon1;

			# output as with concat exon coordinates
			print WW join("\t", $line, $trspt1, $trspt2) . "\n";
		}
	}

	if(uc($asType) eq "SE"){
		
		if($chain eq "+"){

			# obtain trspt1: upstreamES+1 … upstreamEE|downstreamES+1 … downstreamEE
			$trspt1 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);	
			$exon2 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);	
			$trspt1 = join("!", $exon1, $exon2);

			# obtain trspt2: upstreamES+1…upstreamEE|exonStart_0base+1…exonEnd|downstreamES+1…downstreamEE
			$trspt2 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"exonStart_0base"}+1, $asHash{"exonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$trspt2 = join("!", $exon1, $exon2, $exon3);

			# output as with concat exon coordinates
			print WW join("\t", $line, $trspt1, $trspt2) . "\n";

		}elsif($chain eq "-"){

			# obtain trspt1: downstreamES+1…downStreamEE|upstreamES+1…upstreamEE
			$trspt1 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$trspt1 = join("!", $exon1, $exon2);

			# obtain trspt2: downstreamES+1…downstreamEE|exonStart_0base+1…exonEnd|upstreamES+1…upstreamEE
			$trspt2 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"exonStart_0base"}+1, $asHash{"exonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$trspt2 = join("!", $exon1, $exon2, $exon3);

			# output as with concat exon coordinates
			print WW join("\t", $line, $trspt1, $trspt2) . "\n";
		
		}
	}

	if(uc($asType) eq "MXE"){
		
		if($chain eq "+"){

			# obtain trspt1: upstreamES+1…upstreamEE|1stExonStart_0base+1…1stExonEnd|downstreamES+1…downstreamEE
			$trspt1 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"1stExonStart_0base"}+1, $asHash{"1stExonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$trspt1 = join("!", $exon1, $exon2, $exon3);
	
			# obtain trspt2: upstreamES+1…upstreamEE|2stExonStart_0base+1…2stExonEnd|downstreamES+1…downstreamEE
			$trspt2 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"2stExonStart_0base"}+1, $asHash{"2stExonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$trspt2 = join("!", $exon1, $exon2, $exon3);	

			# output as with concat exon coordinates
			print WW join("\t", $line, $trspt1, $trspt2) . "\n";

		}elsif($chain eq "-"){

			# obtain trspt1: downstreamES+1…downstreamEE|1stExonStart_0base+1…1stExonEnd|upstreamES+1…upstreamEE 
			$trspt1 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"1stExonStart_0base"}+1, $asHash{"1stExonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$trspt1 = join("!", $exon1, $exon2, $exon3);

			# obtain trspt2: downstreamES+1…downstreamEE|2stExonStart_0base+1…2stExonEnd|upstreamES+1…upstreamEE
			$trspt2 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"2stExonStart_0base"}+1, $asHash{"2stExonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$trspt2 = join("!", $exon1, $exon2, $exon3);

			# output as with concat exon coordinates
			print WW join("\t", $line, $trspt1, $trspt2) . "\n";
		}
	}
}
close WW;
close FF;
