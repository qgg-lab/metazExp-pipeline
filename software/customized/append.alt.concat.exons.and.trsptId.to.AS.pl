#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV<0){
	print "\tperl $0 \\\n" . 
		"\t\t--asFile total.A5SS.tsv \\\n" .
		"\t\t--asType A5SS \\\n" .
		"\t\t--gtfFile final.complete.trspt.anno.gtf \\\n" .
		"\t\t--outputAsFile total.A5SS.with.alt.concat.exons.and.trsptId.tsv\n\n";
	print "This script is used to append concat exons and trspts for each AS.\n\n";
	exit;
}

my ($asFile, $asType, $gtfFile, $outputAsFile);

GetOptions(
	'asFile=s'=>\$asFile,
	'asType=s'=>\$asType,
	'gtfFile=s'=>\$gtfFile,
	'outputAsFile=s'=>\$outputAsFile
);

if(uc($asType) ne "A5SS" and uc($asType) ne "A3SS" and uc($asType) ne "RI" and uc($asType) ne "SE" and uc($asType) ne "MXE"){
	print "asType must be specified as A5SS, A3SS, RI, SE or MXE\n";
	exit;
}


# obtain exons of each trspt
my (%trsptToExonMultiLineText);
my ($featureLine, @col, $trsptId, @twoDimExon, @srtTwoDimExon, $line, $geneId);
open FF, "<$gtfFile";
while($featureLine=<FF>){
	@col = ();
	@col = split(/\t/, $featureLine);
	
	next if($col[2] ne "exon");

	$trsptId = &getTranscriptIdInAttrs($col[8]);
	$geneId = &getGeneIdInAttrs($col[8]);

	${$trsptToExonMultiLineText{$trsptId}}{"geneId"} = $geneId;
	${$trsptToExonMultiLineText{$trsptId}}{"chrId"} = $col[0];
	${$trsptToExonMultiLineText{$trsptId}}{"chain"} = $col[6];
	${$trsptToExonMultiLineText{$trsptId}}{"exons"} .= join("\t", $col[3], $col[4]) . "\n";
}
close FF;

# obtain concat exon coordinates for each trspt and assign these trspts to cooresponding gene
my (%geneToTrspt); 
# $geneToTrspt{"gene1"} = "rna1:chr1___13751531___13751661___+;chr1___13787497___13788466___+;"
#			"rna2:chr1___12132424___34234534___+;chr1___8634232____98612334___+;"
#			"rna3:chr1___42199424___54239814___+;chr1___9983423____98612334___+;"
my (@trsptId, @exonLineText, $exonLineText, $chrId, $chain);
@trsptId = ();
@trsptId = keys(%trsptToExonMultiLineText);
foreach $trsptId(@trsptId){

	$geneId = ${$trsptToExonMultiLineText{$trsptId}}{"geneId"};
	$chrId = ${$trsptToExonMultiLineText{$trsptId}}{"chrId"};
	$chain = ${$trsptToExonMultiLineText{$trsptId}}{"chain"};
	@exonLineText = ();
	@exonLineText = split(/\n/, ${$trsptToExonMultiLineText{$trsptId}}{"exons"});

	# load exon into @twoDimExon
	@twoDimExon = ();
	foreach $exonLineText(@exonLineText){
		push @twoDimExon, [split '\t', $exonLineText];
	}
	
	# sort exon into @srtTwoDimExon
	if($chain eq "+"){
		@srtTwoDimExon = sort{$a->[0]<=>$b->[0]}@twoDimExon;
	}else{
		@srtTwoDimExon = sort{$b->[0]<=>$a->[0]}@twoDimExon;
	}

	# generate concat exon coordinates for each trspt: rna1:chr1___123___234___+;
	my $trspt_concat_exon_text = "";
	$trspt_concat_exon_text = $trsptId . ":";
	for(my $i=0; $i<=$#srtTwoDimExon; $i++){
		$trspt_concat_exon_text .= join("___", $chrId, $srtTwoDimExon[$i][0], $srtTwoDimExon[$i][1], $chain) . ";";
	}
	
	$geneToTrspt{$geneId} .= $trspt_concat_exon_text . "\n";
}

my @geneId = keys(%geneToTrspt);
open GFFTEXT, ">trspt.coordinate.txt";
foreach $geneId(@geneId){
	print GFFTEXT $geneId . ":\n" . $geneToTrspt{$geneId} . "\n";
}
close GFFTEXT;

# generate local alternative exon coordinates
my (@field, @tmpArr, %asHash, @fieldName);
my ($chain, $chrId,  $i,  $altConcatExons1, $altConcatExons2, $trspt1Ids, $trspt2Ids, $exon1, $exon2, $exon3, $exon4);
open FF, "<$asFile";
$line=<FF>;
chomp($line);
open WW, ">$outputAsFile";
print WW $line . "\t1stAltExons\t1stTrsptIds\t2ndAltExons\t2ndTrsptIds\n";

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

	$geneId = $field[1];
	$geneId=~s/"//g;	

	$chrId = $field[3];
	$chrId = substr($chrId, 3);
	$chain = $field[4];
	
	# read each field into asHash according to fieldName
	for($i=0; $i<=$#fieldName; $i++){
		$asHash{$fieldName[$i]}=$field[$i];
	}

	# output AS position
	if(uc($asType) eq "A5SS"){ # outputs are same with each other when on +/- chain 

		# obtain trspt1: chr1____12345____45678___+;chr1____567890___890123___+
		$altConcatExons1 = "";
		$exon1 = join("___", $chrId, $asHash{"longExonStart_0base"}+1, $asHash{"longExonEnd"}, $chain);
		$exon2 = join("___", $chrId, $asHash{"flankingES"}+1, $asHash{"flankingEE"}, $chain);
		$altConcatExons1 = join(";", $exon1, $exon2);
	
		# obtain trspt2: 
		$altConcatExons2 = "";
		$exon1 = join("___", $chrId, $asHash{"shortES"}+1, $asHash{"shortEE"}, $chain);
		$exon2 = join("___", $chrId, $asHash{"flankingES"}+1, $asHash{"flankingEE"}, $chain);
		$altConcatExons2 = join(";", $exon1, $exon2);
						
	}
	if(uc($asType) eq "A3SS"){ # outputs are same with each other when on +/- chain

		# obtain trspt1
		$altConcatExons1 = "";
		$exon1 = join("___", $chrId, $asHash{"flankingES"}+1, $asHash{"flankingEE"}, $chain);
		$exon2 = join("___", $chrId, $asHash{"longExonStart_0base"}+1, $asHash{"longExonEnd"}, $chain);
		$altConcatExons1 = join(";", $exon1, $exon2);

		# obtain trspt2
		$altConcatExons2 = "";
		$exon1 = join("___", $chrId, $asHash{"flankingES"}+1, $asHash{"flankingEE"}, $chain);
		$exon2 = join("___", $chrId, $asHash{"shortES"}+1, $asHash{"shortEE"}, $chain);
		$altConcatExons2 = join(";", $exon1, $exon2);

	}

	if(uc($asType) eq "RI"){

		if($chain eq "+"){

			# obtain trspt1: upstreamES+1…upstreamEE|downstreamES+1 … downstreamEE
			$altConcatExons1 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$altConcatExons1 = join(";", $exon1, $exon2);

			# obtain trspt2: riExonStart_0base + 1 … riExonEnd
			$altConcatExons2 = "";
			$exon1 = join("___", $chrId, $asHash{"riExonStart_0base"}+1, $asHash{"riExonEnd"}, $chain);
			$altConcatExons2 = join(";", $exon1);

		}elsif($chain eq "-"){
			
			# obtain trspt1:downstreamES+1 … downstreamEE | upstreamES+1…upstreamEE
			$altConcatExons1 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$altConcatExons1 = join(";", $exon1, $exon2);

			# obtain trspt2: riExonStart_0base … riExonEnd
			$altConcatExons2 = "";
			$exon1 = join("___", $chrId, $asHash{"riExonStart_0base"}+1, $asHash{"riExonEnd"}, $chain);
			$altConcatExons2 = $exon1;
		}
	}

	if(uc($asType) eq "SE"){
		
		if($chain eq "+"){

			# obtain trspt1: upstreamES+1 … upstreamEE|downstreamES+1 … downstreamEE
			$altConcatExons1 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);	
			$exon2 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);	
			$altConcatExons1 = join(";", $exon1, $exon2);

			# obtain trspt2: upstreamES+1…upstreamEE|exonStart_0base+1…exonEnd|downstreamES+1…downstreamEE
			$altConcatExons2 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"exonStart_0base"}+1, $asHash{"exonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$altConcatExons2 = join(";", $exon1, $exon2, $exon3);

		}elsif($chain eq "-"){

			# obtain trspt1: downstreamES+1…downStreamEE|upstreamES+1…upstreamEE
			$altConcatExons1 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$altConcatExons1 = join(";", $exon1, $exon2);

			# obtain trspt2: downstreamES+1…downstreamEE|exonStart_0base+1…exonEnd|upstreamES+1…upstreamEE
			$altConcatExons2 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"exonStart_0base"}+1, $asHash{"exonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$altConcatExons2 = join(";", $exon1, $exon2, $exon3);
		
		}
	}

	if(uc($asType) eq "MXE"){
		
		if($chain eq "+"){

			# obtain trspt1: upstreamES+1…upstreamEE|1stExonStart_0base+1…1stExonEnd|downstreamES+1…downstreamEE
			$altConcatExons1 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"1stExonStart_0base"}+1, $asHash{"1stExonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$altConcatExons1 = join(";", $exon1, $exon2, $exon3);
	
			# obtain trspt2: upstreamES+1…upstreamEE|2stExonStart_0base+1…2stExonEnd|downstreamES+1…downstreamEE
			$altConcatExons2 = "";
			$exon1 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"2stExonStart_0base"}+1, $asHash{"2stExonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$altConcatExons2 = join(";", $exon1, $exon2, $exon3);	

		}elsif($chain eq "-"){

			# obtain trspt1: downstreamES+1…downstreamEE|1stExonStart_0base+1…1stExonEnd|upstreamES+1…upstreamEE 
			$altConcatExons1 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"1stExonStart_0base"}+1, $asHash{"1stExonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$altConcatExons1 = join(";", $exon1, $exon2, $exon3);

			# obtain trspt2: downstreamES+1…downstreamEE|2stExonStart_0base+1…2stExonEnd|upstreamES+1…upstreamEE
			$altConcatExons2 = "";
			$exon1 = join("___", $chrId, $asHash{"downstreamES"}+1, $asHash{"downstreamEE"}, $chain);
			$exon2 = join("___", $chrId, $asHash{"2stExonStart_0base"}+1, $asHash{"2stExonEnd"}, $chain);
			$exon3 = join("___", $chrId, $asHash{"upstreamES"}+1, $asHash{"upstreamEE"}, $chain);
			$altConcatExons2 = join(";", $exon1, $exon2, $exon3);

		}
	}
	#print "geneId:$geneId";
	#<STDIN>;
	#print $geneToTrspt{$geneId};
	#<STDIN>;
	#print $altConcatExons1;
	#<STDIN>;
	#print $altConcatExons2;
	#<STDIN>;
	# obtain trsptIds which contain alter concat exons
	$trspt1Ids = &getTrsptIdsByAltConcatExons($geneToTrspt{$geneId}, $altConcatExons1);
	$trspt2Ids = &getTrsptIdsByAltConcatExons($geneToTrspt{$geneId}, $altConcatExons2);
	# output as with concat exon coordinates and corresponding trsptIds
	print WW join("\t", $line, $altConcatExons1, $trspt1Ids, $altConcatExons2, $trspt2Ids) . "\n";
}
close WW;
close FF;

sub getTranscriptIdInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/transcript_id "(.*)"/){
                        return $1;
                }
        }
        return "";
}

sub getGeneIdInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/gene_id "(.*)"/){
                        return $1;
                }
        }
        return "";
}


# obtain trsptIds if they contain alt exons in their complete concate exon coordinates
# $geneToTrspt{"gene1"} = "rna1:chr1___13751531___13751661___+;chr1___13787497___13788466___+;"
#                       "rna2:chr1___12132424___34234534___+;chr1___8634232____98612334___+;"
#                       "rna3:chr1___42199424___54239814___+;chr1___9983423____98612334___+;"
sub getTrsptIdsByAltConcatExons{
	my ($multiTrsptText, $altConcatExons) = @_;
	my (@trsptTextLine, $trsptTextLine, $trsptIds, @tempArr);

	$trsptIds = "";
	@trsptTextLine = split(/\n/, $multiTrsptText);
	foreach $trsptTextLine(@trsptTextLine){
		@tempArr = ();
		@tempArr = split(/:/, $trsptTextLine);
		$trsptIds.=$tempArr[0] . "," if(index($tempArr[1], $altConcatExons)>=0);
	}
	return substr($trsptIds, 0, length($trsptIds) - 1);
}
