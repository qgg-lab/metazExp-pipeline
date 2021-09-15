#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
	" \t\t --referenceGtf /mnt/home/liujind1/workAS/01-cattle/000-prepare-local-datasource/ensembl.annotation.gtf \\\n" . 
	" \t\t --assemblyGtf /mnt/home/liujind1/workAS/01-cattle/000-prepare-local-datasource/refSeq.annotation.with_ensemblSeqId.gtf \\\n" . 
	" \t\t --cuffcompare  /mnt/home/liujind1/software/cufflinks-2.2.1.Linux_x86_64/cuffcompare \\\n" . 
	" \t\t --transcriptMinCov      2 \\\n" .
	" \t\t --transcriptMinLen      2 \\\n" . 
	" \t\t --exonMinCov            0.5 \\\n" .
	" \t\t --minExonNum  2 \\\n" .
	" \t\t --outputTmpDir  /mnt/home/liujind1/workAS/01-cattle/000-prepare-local-datasource/tmpDir \\\n" . 
	" \t\t --outputCombinedGtf /mnt/home/liujind1/workAS/01-cattle/000-prepare-local-datasource/combined.ensembl.refseq.gtf \\\n" . 
	" \t\t --outputCombineStat /mnt/home/liujind1/workAS/01-cattle/000-prepare-local-datasource/combined.ensembl.refseq.sta\n\n\n";

	print "Only append the U and J assembled transcripts into combined annoation.\n" . 
		"In addition, the appended transcripts can be filered with specified minCov, minLen and minExonNum\n\n\n";

	exit;
}

my ($referenceGtf, $assemblyGtf, $cuffcompare, $transcriptMinCov, $transcriptMinLen, $minExonNum, $outputCombinedGtf, $outputCombineStat, $outputTmpDir, $exonMinCov);

GetOptions(
        'referenceGtf=s'=>\$referenceGtf,
        'assemblyGtf=s'=>\$assemblyGtf,
        'cuffcompare=s'=>\$cuffcompare,
	'transcriptMinCov=s'=>\$transcriptMinCov,
	'transcriptMinLen=s'=>\$transcriptMinLen,
	'exonMinCov=s'=>\$exonMinCov,
	'minExonNum=s'=>\$minExonNum,
	'outputTmpDir=s'=>\$outputTmpDir,
	'outputCombinedGtf=s'=>\$outputCombinedGtf,
	'outputCombineStat=s'=>\$outputCombineStat,
);

my (%referenceTrsptFullFeature, %assemblyTrsptFullFeature);
my (%referenceGeneIdToTrsptIdList, %assemblyGeneIdToTrsptIdList);
my (%referenceTrsptIdToExonNum, %assemblyTrsptIdToExonNum);

my (@fields, @attrs, $featureLine);
my ($transcriptId, $geneId);
my (@transcriptId, @exonCoordinateFullLine, @srtExon);

system("mkdir -p " . $outputTmpDir);

#read reference tanscript and exon into hash %referenceTranscript
open STA, ">$outputCombineStat";
open FF, "<$referenceGtf";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	if($fields[2] eq "exon"){
		$transcriptId = &getTranscriptId($fields[8]);
		#register feature into hash %referenceTrsptFullFeature
		$referenceTrsptFullFeature{$transcriptId} .= $featureLine . "\n";
		$referenceTrsptIdToExonNum{$transcriptId} += 1;
	}elsif($fields[2] eq "transcript"){
		#register feature into hash %referenceTrsptFullFeature
		$transcriptId = &getTranscriptId($fields[8]);
		$geneId = &getGeneId($fields[8]);
		$referenceTrsptFullFeature{$transcriptId} .= $featureLine . "\n";
		$referenceGeneIdToTrsptIdList{$geneId}.=$transcriptId . "#";
	}
}
close FF;

#register gene, transcript and exon num of reference annotation into stat file
my @tt = ();
@tt = keys(%referenceGeneIdToTrsptIdList);
my $totalGeneNum = $#tt + 1;
@tt = ();
@tt = keys(%referenceTrsptFullFeature);
my $totalTrsptNum = $#tt + 1;
my $cmd = "grep -P \"\\texon\\t\" $referenceGtf |awk -F \'\\t\' \'{print \$1\"#\"\$4\"#\"\$5\"#\"\$7}\' |sort -u |wc -l";
my $totalExonNum = `$cmd`;
chomp($totalExonNum);
print STA join("\t", $referenceGtf, $totalGeneNum, $totalTrsptNum, $totalExonNum) . "\n";





#read assembly tanscript and exon into hash %assemblyTranscript
my (%transcriptIdToExonMinCov, $exonCov);

open FF, "<$assemblyGtf";
while($featureLine =<FF>){
	chomp($featureLine);
	@fields = ();
	@fields = split(/\t/, $featureLine);
	if($fields[2] eq "exon"){
		$transcriptId = &getTranscriptId($fields[8]);
		#register feature into hash %assemblyTrsptFullFeature
		$assemblyTrsptFullFeature{$transcriptId} .= $featureLine . "\n";
		$assemblyTrsptIdToExonNum{$transcriptId} += 1;

		#get exon cov
		$exonCov = &getExonCov($fields[8]);
		#print "fields[8]:" . $fields[8] . "\n";
		#print "exonCov:" . $fields[8] . "\n";

		#First exon covmust be registered into ExonMinCov hash.
		if(not exists($transcriptIdToExonMinCov{$transcriptId})){
			$transcriptIdToExonMinCov{$transcriptId} = $exonCov;
		}else{
			if($exonCov < $transcriptIdToExonMinCov{$transcriptId}){
				$transcriptIdToExonMinCov{$transcriptId} = $exonCov;
			}
		}

	}elsif($fields[2] eq "transcript"){
		#register feature into hash %assemblyTrsptFullFeature
		$transcriptId = &getTranscriptId($fields[8]);
		$geneId = &getGeneId($fields[8]);
		$assemblyTrsptFullFeature{$transcriptId} .= $featureLine . "\n";
		$assemblyGeneIdToTrsptIdList{$geneId}.=$transcriptId . "#";
	}
}
close FF;

#register gene, transcript and exon num of assembly into sta file
my @tt = ();
@tt = keys(%assemblyGeneIdToTrsptIdList);
my $totalGeneNum = $#tt + 1;
@tt = ();
@tt = keys(%assemblyTrsptFullFeature);
my $totalTrsptNum = $#tt + 1;
my $cmd = "grep -P \"\\texon\\t\" $assemblyGtf |awk -F \'\\t\' \'{print \$1\"#\"\$4\"#\"\$5\"#\"\$7}\' |sort -u |wc -l";
my $totalExonNum = `$cmd`;
chomp($totalExonNum);
print STA join("\t", $assemblyGtf, $totalGeneNum, $totalTrsptNum, $totalExonNum) . "\n";







#compare -----
#register gene, transcript and exon num of assembly annotation into stat file
my $cmd = $cuffcompare . 
	" -r " . $referenceGtf .
	" -o " . $outputTmpDir . "/cmp" .
	" -C " . $assemblyGtf . 
	" &> " . $outputTmpDir . "/log.cuffcompare";
#print $cmd . "\n";
system($cmd);

#<STDIN>;
#open cmp.xxx.tmap file to load j and u transcripts into reference hash
my (@tmpArr, $cmpRltPrefix, $tmapFile, $line);
my ($referenceGeneId, $referenceTrsptId, $assemblyGeneId, $assemblyTrsptId, $cmpCode);

@tmpArr = split(/\//, $assemblyGtf);
pop(@tmpArr);
$cmpRltPrefix = join("/", @tmpArr);

@tmpArr = ();
@tmpArr = split(/\//, $assemblyGtf);
my $fileNameWithoutPath = pop(@tmpArr);
$tmapFile = "";
if($#tmpArr>=0){
	$tmapFile = $cmpRltPrefix . "/cmp." . $fileNameWithoutPath . ".tmap";
}else{
	$tmapFile = "./cmp." . $fileNameWithoutPath . ".tmap";
}
#print $tmapFile . "\n";
#<STDIN>;

my ($cov, $len);
my ($referenceGeneName, $assemblyGeneName);
#read tmap
open FF, "<$tmapFile";
#ref_gene_id ref_id class_code cuff_gene_id cuff_id FMI FPKM FPKM_conf_lo FPKM_conf_hi cov len    major_iso_id ref_match_len
# -          -      u          LOC112447072 rna0    0   0.0  0.00         0.00         0.0 1723    rna0        -
#RCAN1   ENSBTAT00000074355      =       STRG.2  STRG.2.1        100     4.975749        0.000000        0.000000        34.629589       1971
#    STRG.2.2        1971
<FF>;
while($line=<FF>){
	@tmpArr = ();
	@tmpArr = split(/\t/, $line);

	$referenceTrsptId = $tmpArr[1];
	$referenceGeneId = &getGeneId($referenceTrsptFullFeature{$referenceTrsptId});
	$cmpCode = $tmpArr[2];
#	$assemblyGeneId = $tmpArr[3];
	$assemblyTrsptId = $tmpArr[4];
	$assemblyGeneId = &getGeneId($assemblyTrsptFullFeature{$assemblyTrsptId});
	$cov = $tmpArr[9];
	$len = $tmpArr[10];
	
	#filter
	next if(uc($cmpCode) ne "U" and uc($cmpCode) ne "J");
	
	#print $assemblyTrsptId . "\n";
	#print "TranscriptCov:" . $cov . "\n";
	#print "TranscriptExonNum:" . $assemblyTrsptIdToExonNum{$assemblyTrsptId} . "\n";
	#print "exonMinCov:" . $transcriptIdToExonMinCov{$assemblyTrsptId} . "\n";
	#print "Type:" . uc($cmpCode) . "\n";
	#<STDIN>;
	next if($cov < $transcriptMinCov or $len < $transcriptMinLen);
	next if($assemblyTrsptIdToExonNum{$assemblyTrsptId} < $minExonNum);
	next if($transcriptIdToExonMinCov{$assemblyTrsptId} < $exonMinCov);

	#Replace assembly gene Id with reference gene Id
	if(uc($cmpCode) eq "J"){
		$assemblyGeneId =~s/\./\\\./g;
		$assemblyTrsptFullFeature{$assemblyTrsptId}=~s/gene_id "$assemblyGeneId";/gene_id "$referenceGeneId";/g;
		#Append assembly transcript full features into reference hash
		$referenceTrsptFullFeature{$assemblyTrsptId} = $assemblyTrsptFullFeature{$assemblyTrsptId};
		$referenceGeneIdToTrsptIdList{$referenceGeneId}.=$assemblyTrsptId . "#";
	}else{
		$referenceTrsptFullFeature{$assemblyTrsptId} = $assemblyTrsptFullFeature{$assemblyTrsptId};
		$referenceGeneIdToTrsptIdList{$assemblyGeneId}.=$assemblyTrsptId . "#";
	}
}
close FF;

#output all reference transcripts including appended assembly transcripts
open WW, ">$outputCombinedGtf";
@transcriptId = ();
@transcriptId = keys(%referenceTrsptFullFeature);
foreach $transcriptId(@transcriptId){
	print WW $referenceTrsptFullFeature{$transcriptId};
}
close WW;

#register gene, transcript and exon num of combined annotation into stat file
my @tt = ();
@tt = keys(%referenceGeneIdToTrsptIdList);
my $totalGeneNum = $#tt + 1;
@tt = ();
@tt = keys(%referenceTrsptFullFeature);
my $totalTrsptNum = $#tt + 1;
my $cmd = "grep -P \"\\texon\\t\" $outputCombinedGtf |awk -F \'\\t\' \'{print \$1\"#\"\$4\"#\"\$5\"#\"\$7}\' |sort -u |wc -l";
my $totalExonNum = `$cmd`;
chomp($totalExonNum);
print STA join("\t",$outputCombinedGtf, $totalGeneNum, $totalTrsptNum, $totalExonNum) . "\n";

close STA;
#clean temporary output data
#system("rm -rf " . $outputTmpDir);
#@tmpArr = split(/\//, $assemblyGtf);
#if($#tmpArr>0){
#	system("rm -rf " . $cmpRltPrefix . "/cmp.*");
#}else{
#	system("rm -rf cmp.*");
#}

sub getTranscriptId{
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
sub getGeneId{
	my ($singleTrsptAnnoText) = $_[0];
	my (@featureLines, $featureLine, @attrs);
	@featureLines = split(/\n/, $singleTrsptAnnoText);
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

# cov "0.527027";
sub getExonCov{
	my $attrString=$_[0];
	my (@atts, $cov, $att);
	$cov = 0;
	@atts = split(/;/, $attrString);
	foreach $att(@atts){
		if($att=~/cov "(.*)"/){
			$cov = $1;
		}
	}
	return $cov;
}
