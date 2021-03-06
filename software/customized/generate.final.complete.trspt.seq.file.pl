#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputGenomeFile ../001-prepare-local-datasource/genome.fa \\\n" .
		"--inputFinalAnnotGtfFile ./final.complete.trspt.anno.gtf \\\n" . 
		"--inputOriginGtfFileList ../001-prepare-local-datasource/ensembl.gtf \\\n" .
		"--inputOriginGffFileList ../001-prepare-local-datasource/refseq.gff3 \\\n" .
		"--inputProteinFileList ../001-prepare-local-datasource/ensembl.pep.fa,../001-prepare-local-datasource/refseq.pep.fa \\\n" .
		"--outputDNAseqFile ./final.complete.trspt.DNA.fa \\\n" .
		"--outputProteinSeqFile ./final.complete.trspt.pep.fa \\\n" .
		"--outputCdnaSeqFile ./final.complete.trspt.cDNA.fa \n";


	exit;
}

my (%genomeSeq, %proteinSeq, %trsptIdToPepId);
my ($inputGenomeFile, $inputFinalAnnotGtfFile, $inputOriginGtfFileList, $inputOriginGffFileList, $inputProteinFileList);
my ($outputProteinSeqFile, $outputCdnaSeqFile, $outputDNAseqFile);

GetOptions(
        'inputGenomeFile=s'=>\$inputGenomeFile,
        'inputFinalAnnotGtfFile=s'=>\$inputFinalAnnotGtfFile,
	'inputOriginGtfFileList=s'=>\$inputOriginGtfFileList,
	'inputOriginGffFileList=s'=>\$inputOriginGffFileList,
	'inputProteinFileList=s'=>\$inputProteinFileList,
	'outputDNASeqFile=s'=>\$outputDNAseqFile,
	'outputProteinSeqFile=s'=>\$outputProteinSeqFile,
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

# load protein sequence into hash
my @proteinFiles = split(/,/, $inputProteinFileList);
foreach my $proteinFile(@proteinFiles){
	open FF, "<$proteinFile";
	while($line = <FF>){
        	chomp($line);
	        if($line=~/>/){
                	@tt = ();
        	        @tt = split(/ /, $line);
	                if($tt[0]=~/>ref\|(.*?)\|/ or $tt[0]=~/>gi.*\|ref\|(.*?)\|/ or  $tt[0]=~/>(.*)/){
                        	$id = $1;
                	}
        	}else{
	               $proteinSeq{$id} .=uc($line);
#			print ">" . $id . "\n";
#			print $proteinSeq{$id} . "\n";
#			<STDIN>;
        	}
	}
	close FF;
}

print "Finish load protein into hash.\n";

# load the relationship between transcriptId to proteinId in gff file
my (@attr, $transcriptId, $tmpProteinId, $proteinId, $proteinVersion);
my @gffFile=split(/,/, $inputOriginGffFileList);
foreach my $gffFile(@gffFile){
	open FF, "<$gffFile";
	while($line=<FF>){
		chomp($line);
		@tt = split(/\t/, $line);
		next if($tt[2] ne "CDS");

		$transcriptId = "";
		$tmpProteinId = "";
		$proteinId = "";

		@attr = ();
		@attr = split(/;/, $tt[8]);
		foreach my $attr(@attr){
			if($attr=~/Parent=(.*)/){
				$transcriptId = $1;
			}
			if($attr=~/protein_id=(.*)/){
				$proteinId = $1;
			}
		}

		$trsptIdToPepId{$transcriptId} = $proteinId;

#		print $transcriptId . "\t" . $trsptIdToPepId{$transcriptId} ;
#		<STDIN>;

	}
	close FF;
}

print "Finish gff3 into hash.\n";


# load the relationship between transcriptId to proteinId in gff file
my (@attr, $transcriptId, $tmpProteinId, $proteinId, $proteinVersion);
my @gtfFile=split(/,/, $inputOriginGtfFileList);
foreach my $gtfFile(@gtfFile){
	open FF, "<$gtfFile";
	while($line=<FF>){
		chomp($line);
		@tt = split(/\t/, $line);
		next if($tt[2] ne "CDS");

		$transcriptId = "";
		$tmpProteinId = "";
		$proteinId = "";

		@attr = ();
		@attr = split(/;/, $tt[8]);
		foreach my $attr(@attr){
			if($attr=~/transcript_id \"(.*)\"/){
				$transcriptId = $1;
			}
			if($attr=~/protein_id \"(.*)\"/){
				$proteinId = $1;
			}
			if($attr=~/protein_version \"(.*)\"/){
				$proteinVersion = $1;
			}
		}

		$trsptIdToPepId{$transcriptId} = $proteinId . "." . $proteinVersion;
		#print $transcriptId . "\t" . $trsptIdToPepId{$transcriptId} ;
		#<STDIN>;
	}
	close FF;
}

print "Finish gtf into hash.\n";

##############################
#			     #
# load isoform into hash     #
#			     #
##############################
my (%isoform, $geneId, $geneSymbol, $isoformId, $isoformSymbol);
my ($featureLine, @cols);

open FF, "<$inputFinalAnnotGtfFile";

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
 
# generate DNA seq for isoform and distinguish intron and exon
open RNA, ">$outputCdnaSeqFile";
open PEP, ">$outputProteinSeqFile";
open DNA, ">$outputDNAseqFile";

my (@isoformId, @exonIntron, @tmpExon, @exon, $exon, $cDNAseq, $exonSeq);
@isoformId = keys(%isoform);

foreach $isoformId(@isoformId){

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
			$exonNum++;
		}
	}

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


	# output DNA sequence(including intron and exon)
	my $DNAseq = "";
	$DNAseq = uc(substr($genomeSeq{${$isoform{$isoformId}}{"chr"}}, ${$isoform{$isoformId}}{"start"} -1 , ${$isoform{$isoformId}}{"end"} - ${$isoform{$isoformId}}{"start"} + 1));
	if(${$isoform{$isoformId}}{"strand"} eq "-"){
		$DNAseq = reverse($DNAseq);
		$DNAseq=~tr/ACGT/TGCA/;
	}
	print DNA ">" . $isoformId . " transcript_name:". ${$isoform{$isoformId}}{"isoformSymbol"} . " gene_id:" . ${$isoform{$isoformId}}{"geneId"} . " gene_name:" . ${$isoform{$isoformId}}{"geneSymbol"} . "\n";
	print DNA $DNAseq . "\n";

	# output protein sequence
	if(exists($proteinSeq{$trsptIdToPepId{$isoformId}})){
		print PEP ">" . $isoformId . " transcript_name:". ${$isoform{$isoformId}}{"isoformSymbol"} . " gene_id:" . ${$isoform{$isoformId}}{"geneId"} . " gene_name:" . ${$isoform{$isoformId}}{"geneSymbol"} . " protein_id:" . $trsptIdToPepId{$isoformId} . "\n";
		print PEP $proteinSeq{$trsptIdToPepId{$isoformId}} . "\n";
	}
}

close RNA;
close PEP;
close DNA;

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
