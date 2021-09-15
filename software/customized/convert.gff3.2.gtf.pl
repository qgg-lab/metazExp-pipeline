#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV <0){
	print "\n$0 \\\n" . 
		"\t\t --inputGff3File refseq.gff3  \\\n" . 
		"\t\t --detectSourceFlag no \\\n" .
		"\t\t --sourceName refseq \\\n" . 
		"\t\t --outputGtfFile refseq.gtf   \n";
	print "\nThis script is used to convert gff3 to gtf.\n\n";
	exit;
}
my ($inputGff3File, $outputGtfFile, $detectSourceFlag, $sourceName);

$sourceName = "";
GetOptions(
	'inputGff3File=s'=>\$inputGff3File,
	'detectSourceFlag=s'=>\$detectSourceFlag,
	'sourceName=s'=>\$sourceName,
	'outputGtfFile=s'=>\$outputGtfFile,
);

if(uc($detectSourceFlag) eq "NO" and $sourceName eq ""){
	print STDERR "If you don't use gene source name in gff3, you must specify sourceName option\n";
	exit;
}

my (%exon, %CDS);
my ($featureLine, @cols, @attr);
my (@geneId );

my ($geneId);

my ($i, $j, $exonId, $exonNum, $cdsNum);

my (%gene, $geneId, $geneBiotype, $geneName, $geneAttrString, $geneSource);

$geneSource = "";
open FF, "<$inputGff3File";
while($featureLine = <FF>){
	chomp($featureLine);
	next if($featureLine=~/#/);
	@cols = split(/\t/, $featureLine);

	if($cols[2] eq "gene"){
		# gather gene_id, gene_name, gene_biotype and gene_source 
		# and then combine these attrs to build the 9th col in gene feature.
		$geneId = &getId($cols[8]);
		$geneBiotype = &getBiotype($cols[8]);
		$geneName = &getName($cols[8]);

		$geneAttrString  = "gene_id \"" . $geneId . "\";";
		$geneAttrString .= " gene_name \"" . $geneName . "\";" if($geneName ne "");
		$geneAttrString .= " gene_biotype \"" . $geneBiotype . "\";" if($geneBiotype ne "");

		if(uc($detectSourceFlag) eq "NO"){
			$geneSource = $sourceName;
		}else{
			$geneSource = &getGeneSource($cols[8]);
		}
		$geneAttrString .= " gene_source \"" . $geneSource . "\";";

		${$gene{$geneId}}{"attrs"} = $geneAttrString;
		${$gene{$geneId}}{"feature"} = join("\t", $cols[0], $geneSource, "gene", $cols[3], $cols[4], $cols[5], $cols[6], $cols[7], $geneAttrString); 
		${$gene{$geneId}}{"chain"} = $cols[6];

		${$gene{$geneId}}{"transcriptNum"} = 0;
		${$gene{$geneId}}{"transcriptIdList"} = "";
	}
}
close FF;



#gather transcript level feature, such as lnc_RNA, C_gene_seqment, V_gene_segment, transcript, mRNA and so on.
my (%transcript, $transcriptId, $transcriptBiotype, $transcriptName, $transcriptAttrString);
open FF, "<$inputGff3File";
while($featureLine = <FF>){
	chomp($featureLine);
	next if($featureLine=~/#/);
	@cols = split(/\t/, $featureLine);
	if($cols[2] eq "antisense_RNA" or $cols[2] eq "C_gene_segment" or $cols[2] eq "guide_RNA" 
		or $cols[2] eq "lnc_RNA" or $cols[2] eq "primary_transcript" or $cols[2] eq "RNase_MRP_RNA"
		or $cols[2] eq "rRNA" or $cols[2] eq "snoRNA" or $cols[2] eq "snRNA" or $cols[2] eq "SRP_RNA" 
		or $cols[2] eq "telomerase_RNA" or $cols[2] eq "tRNA" or $cols[2] eq "V_gene_segment"
		or $cols[2] eq "transcript" or $cols[2] eq "mRNA"){
		
		$transcriptId = &getId($cols[8]);
		$transcriptBiotype = &getBiotype($cols[8]);
		$transcriptName = &getName($cols[8]);
		$geneId = &getParent($cols[8]);

		if(exists(${$gene{$geneId}}{"feature"})){

			${$gene{$geneId}}{"transcriptNum"}++; 
			${$gene{$geneId}}{"transcriptIdList"} .= $transcriptId . "#";			

			$transcriptAttrString = "transcript_id \"" . $transcriptId . "\";";
			$transcriptAttrString .= " transcript_name \"" . $transcriptName . "\";" if($transcriptName ne "");
			$transcriptAttrString .= " transcript_biotype \"" . $transcriptBiotype . "\";" if($transcriptBiotype ne "");
			$transcriptAttrString .= " transcript_source \"" . $geneSource . "\";";

			#gene attrs should be put in front of transcript attrs
			${$transcript{$transcriptId}}{"attrs"} = ${$gene{$geneId}}{"attrs"} . " " . $transcriptAttrString;	

			#generate full transcript feature including 9 cols
			${$transcript{$transcriptId}}{"feature"} = join("\t", $cols[0], $geneSource, "transcript", $cols[3], $cols[4], $cols[5], $cols[6], $cols[7], ${$transcript{$transcriptId}}{"attrs"});

			#specify the chain of transcript
			${$transcript{$transcriptId}}{"chain"} = $cols[6];

			#initate exonNum, exon features and cds feature text
			${$transcript{$transcriptId}}{"exonNum"} = 0;
			${$transcript{$transcriptId}}{"exonFeatures"} = "";
			${$transcript{$transcriptId}}{"cdsFeatures"} = "";
		}
	}
}
close FF;



#gather exon feature and append into transcript exon feature
my $exonId;
my ($exonAttrString);
$exonId = 0;
open FF, "<$inputGff3File";
while($featureLine = <FF>){
	chomp($featureLine);
	next if($featureLine=~/#/);
	@cols = split(/\t/, $featureLine);

	if($cols[2] eq "exon"){

		$exonId++;
		$transcriptId = &getParent($cols[8]);

		#discard orphan exon
		if(exists(${$transcript{$transcriptId}}{"feature"})){

			${$transcript{$transcriptId}}{"exonNum"}++;

			#exon_id is numbered in whole genome scale
			$exonAttrString = "exon_id \"exon" . $exonId . "\";";

			#generate exon the 9th col. The transcript attrs should be put in front of exon attrs.
			$exonAttrString = ${$transcript{$transcriptId}}{"attrs"} . " " . $exonAttrString;

			#generate exon feature
			${$transcript{$transcriptId}}{"exonFeatures"} .= join("\t", $cols[0], $geneSource, $cols[2], $cols[3], $cols[4], $cols[5], $cols[6], $cols[7], $exonAttrString) . "\n";
		}
	}
}
close FF;


#gather cds feature 
my $cdsId;
$cdsId = 0;
my ($cdsAttrString);
open FF, "<$inputGff3File";
while($featureLine = <FF>){
	chomp($featureLine);
	next if($featureLine=~/#/);
	@cols = split(/\t/, $featureLine);

	if($cols[2] eq "CDS"){
		$cdsId++;
		$transcriptId = &getParent($cols[8]);

		#discard orphan cds
		if(exists(${$transcript{$transcriptId}}{"feature"})){

			#cds_id is numbered in whole genome scale
			$cdsAttrString = "CDS_id \"CDS" . $cdsId . "\";";

			#transcript attrs should be put in front of cds attrs
			$cdsAttrString = ${$transcript{$transcriptId}}{"attrs"} . " " . $cdsAttrString;

			#generate cds feature including 9 cols
			${$transcript{$transcriptId}}{"cdsFeatures"} .= join("\t", $cols[0], $geneSource, $cols[2], $cols[3], $cols[4], $cols[5], $cols[6], $cols[7], $cdsAttrString) . "\n";
		}		
	}
}
close FF;

#extract geneId by coordinate
my $cmd = "grep -P \"\\tgene\\t\" $inputGff3File |sort -t'	' -k1,1 -k4,4n |awk -F \'\\t\' \'{print \$9}\' |awk -F \';\' \'{print \$1}\' |awk -F \'ID=\' \'{print \$2}\'";
my $geneIdList = `$cmd`;
@geneId = split("\n", $geneIdList);

#output sorted gene
open WW, ">$outputGtfFile";
my @transcriptId = ();
my (@exonCdsFields, @exonFeatures, @cdsFeatures);

#list each geneId
foreach $geneId (@geneId){

	#discard gene without transcripts
	next if(${$gene{$geneId}}{"transcriptIdList"} eq "");

	#first output gene feature
	print WW ${$gene{$geneId}}{"feature"} . "\n";

	#extract each transcript in this gene
	@transcriptId = ();
	@transcriptId = split(/#/, ${$gene{$geneId}}{"transcriptIdList"});
	foreach $transcriptId(@transcriptId){

		#discard transcript without exon
		next if(${$transcript{$transcriptId}}{"exonNum"}==0);
		
		#output transcript feature
		print WW ${$transcript{$transcriptId}}{"feature"} . "\n";

		#combine exon and cds together and sort them by coordinate

		@exonCdsFields = ();
		
		#read exon coordinate into two dim array
		@exonFeatures = ();
		@exonFeatures = split(/\n/, ${$transcript{$transcriptId}}{"exonFeatures"});

		foreach my $featureLine(@exonFeatures){
			push @exonCdsFields, [split '\t', $featureLine];
		}

		#read cds coordinate into two dim array
		@cdsFeatures = ();
		@cdsFeatures = split(/\n/, ${$transcript{$transcriptId}}{"cdsFeatures"});
		foreach my $featureLine(@cdsFeatures){
			push @exonCdsFields, [split '\t', $featureLine];
		}

		#sort cds and exon according to coordinate
		if(${$transcript{$transcriptId}}{"chain"} eq "+"){
			@exonCdsFields = sort{$a->[3]<=>$b->[3]}@exonCdsFields;
		}else{
			@exonCdsFields = sort{$b->[3]<=>$a->[3]}@exonCdsFields;
		}

		#give exon_number within transcript and then output full feature for exon and cds
		my $exonNum;
		$exonNum = 0;
		for(my $ll=0; $ll<=$#exonCdsFields; $ll++){
			if($exonCdsFields[$ll][2] eq "exon"){
				$exonNum++;
				$exonCdsFields[$ll][8] .= " exon_number \"" . $exonNum . "\";";
			}
			print WW join("\t", $exonCdsFields[$ll][0], $exonCdsFields[$ll][1], $exonCdsFields[$ll][2], $exonCdsFields[$ll][3], $exonCdsFields[$ll][4], $exonCdsFields[$ll][5], $exonCdsFields[$ll][6], $exonCdsFields[$ll][7], $exonCdsFields[$ll][8]) . "\n";
		}
	}
}
close WW;


sub getId{
	my $col9 = $_[0];
	my (@attr, $i, $id);
	$id = "";
	@attr = ();
	@attr = split(/;/, $col9);
	for($i=0; $i<=$#attr; $i++){
		if($attr[$i]=~/ID=(.*)/){
			$id = $1;
			last;
		}
	}
	return $id;
}

sub getName{
	my $col9 = $_[0];
	my (@attr, $i, $name);
	$name = "";
	@attr = ();
	@attr = split(/;/, $col9);
	for($i=0; $i<=$#attr; $i++){
		if($attr[$i]=~/Name=(.*)/){
			$name = $1;
			last;
		}
	}
	return $name;
}

sub getParent{
	my $col9 = $_[0];
	my (@attr, $i, $parent);
	$parent = "";
	@attr = ();
	@attr = split(/;/, $col9);
	for($i=0; $i<=$#attr; $i++){
		if($attr[$i]=~/Parent=(.*)/){
			$parent = $1;
		}
	}
	return $parent;
}
sub getBiotype{
	my $col9 = $_[0];
	my (@attr, $i, $biotype);
	$biotype = "";
	@attr = ();
	@attr = split(/;/, $col9);
	for($i=0; $i<=$#attr; $i++){
		if($attr[$i]=~/_biotype=(.*)/){
			$biotype = $1;
			last;
		}
	}
	return $biotype;
}
sub getGeneSource{
	my $col9 = $_[0];
	my (@attr, $i, $source);
	$source = "";
	@attr = ();
	@attr = split(/;/, $col9);
	for($i=0; $i<=$#attr; $i++){
		if($attr[$i]=~/source=(.*)/){
			$source = $1;
			last;
		}
	}
	return $source;

}
