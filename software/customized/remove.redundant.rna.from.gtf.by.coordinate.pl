#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV <0){
	print "\n$0 \\\n" . 
		"\t\t --inputGtfFile refseq.gtf  \\\n" . 
		"\t\t --outputSrtGtfFile new.refseq.gtf \\\n" . 
		"\t\t --outputTmpDir srtTmp \n";
	print "\nThis script is used to remove redundant rna from gtf file according to coordinate.\n\n";
	exit;
}
my ($inputGtfFile, $outputSrtGtfFile, $outputTmpDir);

GetOptions(
	'inputGtfFile=s'=>\$inputGtfFile,
	'outputSrtGtfFile=s'=>\$outputSrtGtfFile,
	'outputTmpDir=s'=>\$outputTmpDir,
);

my (%exon, %CDS);
my ($featureLine, @cols, @attr);
my (@geneId );

my ($geneId, $transcriptId);

my ($i, $j, $exonId, $exonNum, $cdsNum);

my (%gene, $geneId, $geneBiotype, $geneName, $geneAttrString, $geneSource);

#mkdir temporary directory for generate temporate data file
system("mkdir -p " . $outputTmpDir);

my (%trsptIdToExonCoordMltipleLineText);
my (%concatExonCoordToTrsptId);
#read exon coordinates of each transcript and generate sorted concated exon coordinates
open FF, "<$inputGtfFile";
while($featureLine = <FF>){
	@cols = ();
	@cols = split(/\t/, $featureLine);
	
	#discard non-gene feature
	next if($cols[2] ne "exon");
	
	#get transcript Id
	$transcriptId = &getTranscriptIdInAttrs($cols[8]);
	
	#put coordinates into two dim array: chr, start, end, chain
	$trsptIdToExonCoordMltipleLineText{$transcriptId} .= join("\t", $cols[0], $cols[3], $cols[4], $cols[6]) . "\n";
}
close FF;

#generate trsptIdToConcatexonCoord hash
my (@exonCoordTwoDimArray, @exonCoordLine, $exonCoordLine);
my ($concatExonCoordString, @concatExonCoordString);
my @transcriptId = ();
@transcriptId = keys(%trsptIdToExonCoordMltipleLineText);
foreach $transcriptId(@transcriptId){
	
	#get exon coordinate line in multiple text
	@exonCoordLine = ();
	@exonCoordLine = split(/\n/, $trsptIdToExonCoordMltipleLineText{$transcriptId});

	#put exon coordinates into @exonCoordTwoDimArray
	@exonCoordTwoDimArray = ();
	foreach $exonCoordLine(@exonCoordLine){
		push @exonCoordTwoDimArray, [split ' ', $exonCoordLine];		
	}

	#sort exon coord in two dim array
	@exonCoordTwoDimArray = sort{$a->[1]<=>$b->[1]}@exonCoordTwoDimArray;

	#generate concated exon coord string
	$concatExonCoordString = "";
	for(my $i=0; $i<=$#exonCoordTwoDimArray; $i++){
		$concatExonCoordString.=join(" ", $exonCoordTwoDimArray[$i][0], $exonCoordTwoDimArray[$i][1], $exonCoordTwoDimArray[$i][2], $exonCoordTwoDimArray[$i][3]) . "#";
	}
	
	#register concated exon coordinates into hash. If multiple transcripts with same exons, only last transcript will be
	#registered into hash.
	$concatExonCoordToTrsptId{$concatExonCoordString}=$transcriptId;		
}

#register non-redundant transcriptId into hash
my (%transcriptIdExistHash);
@concatExonCoordString = ();
@concatExonCoordString = keys(%concatExonCoordToTrsptId);
foreach $concatExonCoordString(@concatExonCoordString){
	$transcriptId = $concatExonCoordToTrsptId{$concatExonCoordString};
	$transcriptIdExistHash{$transcriptId}=1;
}


#open gene coordiante file
open GWW, ">" . $outputTmpDir . "/gene.coordinate.tsv";
#read gene feature into gene hash and output gene coordiante into temporary for sort later.
open FF, "<$inputGtfFile";
while($featureLine = <FF>){
	chomp($featureLine);
	next if($featureLine=~/#/);
	@cols = ();
	@cols = split(/\t/, $featureLine);

	#discard non-gene feature
	next if($cols[2] ne "gene");
	
	#get geneId
	$geneId = &getGeneIdInAttrs($cols[8]);

	#register gene feature into gene hash
	${$gene{$geneId}}{"geneFeature"} = $featureLine;
	${$gene{$geneId}}{"chain"} = $cols[6];

	#output gene coordinate
	print GWW join(" ",$geneId, $cols[0], $cols[3], $cols[6]) . "\n";
}
close FF;
close GWW;

#gather transcript feature into hash and build relationship between transcript and gene.
my (%transcript, $transcriptId, $transcriptBiotype, $transcriptName, $transcriptAttrString);
open FF, "<$inputGtfFile";
while($featureLine = <FF>){
	@cols = split(/\t/, $featureLine);
	#discard non-transcript feature
	next if($cols[2] ne "transcript");

	#get gene and transcript id
	$geneId = &getGeneIdInAttrs($cols[8]);
	$transcriptId = &getTranscriptIdInAttrs($cols[8]);

	#register transcript id into gene
	if(-s $outputTmpDir . "/gene.coordinate.tsv"){

		#if gene feature exists, transcript are register into cooresponding gene
		${$gene{$geneId}}{"transcriptList"} .= join(" ", $transcriptId, $cols[0], $cols[3]) . "\n";

	}else{

		#if doesn't exist gene feature, all transcripts are registered into only one gene named "ComprehensiveSet"
		${$gene{"ComprehensiveSe"}}{"transcriptList"} .= join(" ", $transcriptId, $cols[0], $cols[3]) . "\n";

	}

	#register transcript chain
	${$transcript{$transcriptId}}{"chain"} = $cols[6];

	#append transcript feature into transcript hash
	${$transcript{$transcriptId}}{"transcriptFeature"} = $featureLine;
}
close FF;

#gather exon and cds feature and append into transcript exoncds feature
my $exonId;
my ($exonAttrString);
$exonId = 0;
open FF, "<$inputGtfFile";
while($featureLine = <FF>){

	#split feature into cols
	@cols = split(/\t/, $featureLine);

	#discard non-exon and non-cds feature
	next if($cols[2] ne "exon" and $cols[2] ne "CDS" );

	#get transcript id
	$transcriptId = &getTranscriptIdInAttrs($cols[8]);

	#append feature into transcript exonCdsMultipleTextFeature
	${$transcript{$transcriptId}}{"exonCdsFeatureText"}.=$featureLine
}
close FF;


#sort geneId by coordinate
my (@transcriptCoordLine, @transcriptCoordTwoDimArr);
my (@exonCdsFeatureLine, @exonCdsFeatureTwoDimArr);

open GWWW, ">$outputSrtGtfFile";

@geneId = ();
my $geneIdList = "";
if(-s $outputTmpDir . "/gene.coordinate.tsv"){
	$geneIdList = "sort -t\' \' -k2,2 -k3,3n " . $outputTmpDir . "/gene.coordinate.tsv | awk -F \' \' \'{print \$1}\'";
	$geneIdList = `$geneIdList`;
	@geneId = split("\n", $geneIdList);
}else{
	$geneId[0] = "ComprehensiveSe";
}

foreach $geneId(@geneId){
	chomp($geneId);
	#output gene feature
	print GWWW ${$gene{$geneId}}{"geneFeature"} . "\n";
	
	#extract transcripts
	@transcriptCoordTwoDimArr = ();

	@transcriptCoordLine = (); 
	@transcriptCoordLine = split(/\n/, ${$gene{$geneId}}{"transcriptList"});

	foreach my $transcriptCoordLine(@transcriptCoordLine){
		push @transcriptCoordTwoDimArr, [split ' ', $transcriptCoordLine];
	}


	#sort transcript by coordintate
	if($#geneId > 0){
		#if there are gene features in annotation, transcripts in a specified gene are sorted by ascend when chain is +
		#sorted by descend when chain is -
		if(${$gene{$geneId}}{"chain"} eq "+"){
			@transcriptCoordTwoDimArr = sort{$a->[1] cmp $b->[1] or $a->[2]<=>$b->[2]}@transcriptCoordTwoDimArr;
		}else{
			@transcriptCoordTwoDimArr = sort{$a->[1] cmp $b->[1] or $b->[2]<=>$a->[2]}@transcriptCoordTwoDimArr;
		}
	}else{
		#if no gene features in annotation, all transcript together are sorted by ascend
		@transcriptCoordTwoDimArr = sort{$a->[1] cmp $b->[1] or $a->[2]<=>$b->[2]}@transcriptCoordTwoDimArr;
	}

	#output each transcript by sorted transcritId
	for(my $i=0; $i<=$#transcriptCoordTwoDimArr; $i++){
		
		#get transcriptId
		$transcriptId = $transcriptCoordTwoDimArr[$i][0];

		#if transcript is an redundant one, then discard it
		next if(not $transcriptIdExistHash{$transcriptId});

		#output transcript feature		
		print GWWW ${$transcript{$transcriptId}}{"transcriptFeature"};

		#get exonCds feature from transcript hash
		@exonCdsFeatureLine = ();
		@exonCdsFeatureLine = split(/\n/, ${$transcript{$transcriptId}}{exonCdsFeatureText});
		
		#split 9 cols of exon and cds features
		@exonCdsFeatureTwoDimArr = ();
		foreach my $exonCdsFeatureLine(@exonCdsFeatureLine){
			push @exonCdsFeatureTwoDimArr, [split '\t', $exonCdsFeatureLine];
		}

		#sort features by coordinate
		if(${$transcript{$transcriptId}}{"chain"} eq "+"){
			@exonCdsFeatureTwoDimArr = sort{$a->[3]<=>$b->[3]}@exonCdsFeatureTwoDimArr;
		}else{
			@exonCdsFeatureTwoDimArr = sort{$b->[3]<=>$a->[3]}@exonCdsFeatureTwoDimArr;
		}

		#output cds and exon features of this transcript
		for(my $m=0; $m<=$#exonCdsFeatureTwoDimArr; $m++){
			print GWWW join("\t", $exonCdsFeatureTwoDimArr[$m][0], $exonCdsFeatureTwoDimArr[$m][1], $exonCdsFeatureTwoDimArr[$m][2], $exonCdsFeatureTwoDimArr[$m][3], $exonCdsFeatureTwoDimArr[$m][4], $exonCdsFeatureTwoDimArr[$m][5], $exonCdsFeatureTwoDimArr[$m][6], $exonCdsFeatureTwoDimArr[$m][7], $exonCdsFeatureTwoDimArr[$m][8]) . "\n";
		}
	}
}
close GWWW;

system("rm -rf ". $outputTmpDir);


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


sub getGeneIdInMultipleFeatures{
        my ($multipleFeatureText) = $_[0];
        my (@featureLines, $featureLine, @attrs);

        @featureLines = split(/\n/, $multipleFeatureText);

        foreach $featureLine(@featureLines){

                @attrs = ();
                @attrs = split(/;/, $featureLine);
                foreach my $attr(@attrs){
                        if($attr=~/gene_id "(.*)"/){
                                return $1;
                        }
                }

        }
        return "";
}
