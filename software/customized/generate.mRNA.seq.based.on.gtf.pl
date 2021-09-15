#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputGenomeFile ../001-prepare-local-datasource/genome.fa \\\n" .
		"--inputGtfFile ./final.complete.trspt.anno.gtf \\\n" . 
		"--outputCdnaSeqFile ./final.complete.trspt.cDNA.fa \n";

	exit;
}

my (%genomeSeq);
my ($inputGenomeFile, $inputGtfFile);
my ($outputCdnaSeqFile);

GetOptions(
        'inputGenomeFile=s'=>\$inputGenomeFile,
        'inputGtfFile=s'=>\$inputGtfFile,
        'outputCdnaSeqFile=s'=>\$outputCdnaSeqFile,
);

# load genome sequence into hash
my ($line, @tt, $id);
open FF, "<$inputGenomeFile";
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
print "Finish load genome into hash.\n";

# load gtf into hash
my (%isoform, $geneId, $geneSymbol, $isoformId, $isoformSymbol);
my ($featureLine, @cols);

open FF, "<$inputGtfFile";

while($featureLine = <FF>){
	chomp($featureLine);
	next if($featureLine=~/#/);
	@cols = split(/\t/, $featureLine);

	$geneId = &getGeneIdInAttrs($cols[8]);
	$geneSymbol = &getGeneSymbolInAttrs($cols[8]);
	$isoformId = &getIsoformIdInAttrs($cols[8]);
	$isoformSymbol = &getIsoformSymbolInAttrs($cols[8]);

	#register isoform chr, strand, start, end and origin
	if($cols[2] eq "transcript"){
		${$isoform{$isoformId}}{"chr"} = $cols[0];
		${$isoform{$isoformId}}{"strand"} = $cols[6];
		${$isoform{$isoformId}}{"start"} = $cols[3];
		${$isoform{$isoformId}}{"end"} = $cols[4];
		${$isoform{$isoformId}}{"geneId"} = $geneId;
		${$isoform{$isoformId}}{"geneSymbol"} = $geneSymbol;
		${$isoform{$isoformId}}{"isoformSymbol"} = $isoformSymbol;

	}elsif($cols[2] eq "exon"){

		${$isoform{$isoformId}}{"exonSeries"} .= $cols[3] . ".." . $cols[4] . ",";
	}
}
close FF;

#my @trsptId = keys(%isoform);
#print $#trspId;
#exit;
# generate DNA seq for isoform and distinguish intron and exon
open RNA, ">$outputCdnaSeqFile";

my (@isoformId, @exonIntron, @tmpExon, @exon, $exon, $cDNAseq, $exonSeq);
@isoformId = keys(%isoform);

print "finish load gtf\n";
foreach $isoformId(@isoformId){
#	print "isoformId:" . $isoformId;
#	<STDIN>;
#	print "exonSeries:" . ${$isoform{$isoformId}}{"exonSeries"} . "\n";
#	<STDIN>;
	# process exonSeries to cut final comma
	${$isoform{$isoformId}}{"exonSeries"} = substr(${$isoform{$isoformId}}{"exonSeries"}, 0, length(${$isoform{$isoformId}}{"exonSeries"})-1);

	# to generate exon
	@exon = ();
	@tmpExon = ();
	@tmpExon = split(/,/, ${$isoform{$isoformId}}{"exonSeries"});
	my $exonNum=0;
	foreach $exon(@tmpExon){
		if($exon=~/(\d+)\.\.(\d+)/){
			$exon[$exonNum][0]= $1;
			$exon[$exonNum][1]= $2;
#			print $exon[$exonNum][0] . "\t" . $exon[$exonNum][1] . "\n";
			$exonNum++;
		}
	}
#	<STDIN>;
	# output RNA sequence(cDNA)
	$cDNAseq = "";
	$exonSeq = "";
	for(my $i=0; $i<=$#exon; $i++){
		$exonSeq = uc(substr($genomeSeq{${$isoform{$isoformId}}{"chr"}}, $exon[$i][0]-1, $exon[$i][1] - $exon[$i][0] + 1));
		if(${$isoform{$isoformId}}{"strand"} eq "-"){
			$exonSeq = reverse($exonSeq);
			$exonSeq=~tr/ACGT/TGCA/;
		}
		$cDNAseq .= $exonSeq;
	}
	print RNA ">" . $isoformId . " transcript_name:". ${$isoform{$isoformId}}{"isoformSymbol"} . " gene_id:" . ${$isoform{$isoformId}}{"geneId"} . " gene_name:" . ${$isoform{$isoformId}}{"geneSymbol"} . "\n";
	print RNA $cDNAseq . "\n";


}

close RNA;

sub getGeneIdInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/gene_id "(.*)"/){
                        return $1;
                }
        }
        return "NA";
}

sub getIsoformIdInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/transcript_id "(.*)"/){
                        return $1;
                }
        }
        return "NA";
}


sub getGeneSymbolInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/gene_name "(.*)"/){
                        return $1;
                }
        }
        return "NA";
}

sub getIsoformSymbolInAttrs{
        my ($attrsString) = $_[0];
        my (@attrs, $attr);
        @attrs = split(/;/, $attrsString);
        foreach $attr(@attrs){
                if($attr=~/transcript_name "(.*)"/){
                        return $1;
                }
        }
        return "NA";
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
