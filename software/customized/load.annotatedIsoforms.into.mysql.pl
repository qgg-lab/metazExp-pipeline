#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--dbName asdb \\\n" . 
		"--dbUser lsas \\\n" . 
		"--dbPWD njaulsas2019 \\\n" . 
		"--species \"Bos taurus\" \\\n" .
		"--inputAnnoGtfFile ./final.complete.trspt.anno.gtf \\\n" .
		"--inputCdnaSeqFile ./final.complete.trspt.cDNA.fa \\\n".
		"--inputDNAseqFile ./final.complete.trspt.DNA.fa \\\n" .
		"--inputPepSeqFile ./final.complete.trspt.pep.fa \\\n" . 
		"--inputPfamGoAnnoFile ./final.complete.trspt.pfam.Go.anno.tsv \n";

	print "Note: the annotated gtf file must be sorted by exon coordinate\n\n";
	exit;
}

my ($dbName, $dbUser, $dbPWD);
my ($species, $inputAnnoGtfFile);
my ($inputAnnoGtfFile, $inputDNAseqFile, $inputPepSeqFile, $inputCdnaSeqFile);
my ($inputPfamGoAnnoFile);

GetOptions(
        'dbName=s'=>\$dbName,
        'dbUser=s'=>\$dbUser,
        'dbPWD=s'=>\$dbPWD,
        'species=s'=>\$species,
        'inputAnnoGtfFile=s'=>\$inputAnnoGtfFile,
	'inputCdnaSeqFile=s'=>\$inputCdnaSeqFile,
	'inputDNAseqFile=s'=>\$inputDNAseqFile, 
	'inputPepSeqFile=s'=>\$inputPepSeqFile,
	'inputPfamGoAnnoFile=s'=>\$inputPfamGoAnnoFile,
);

#my $dbh = DBI->connect("DBI:mysql:database=$dbName", $dbUser, $dbPWD);
my $dbh = DBI->connect("DBI:mysql:database=liujind1_db1;host=db-01;mysql_skip_secure_auth=1", "liujind1", "DHQn4TdaUE=9h#9");
$dbh->{mysql_auto_reconnect} = 1;
my $query;
##############################
#			     #
# load gtf into hash         #
#			     #
##############################
my (%isoform, $geneId, $geneSymbol, $isoformId, $isoformSymbol);
my ($featureLine, @cols);

# print "Begin load gtf into hash.\n";
open FF, "<$inputAnnoGtfFile";

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
		if($cols[1] eq "StringTie"){
			${$isoform{$isoformId}}{"annoOrigin"} = "RNAseq";
		}else{
			${$isoform{$isoformId}}{"annoOrigin"} = $cols[1];
		}
	}elsif($cols[2] eq "exon"){

		${$isoform{$isoformId}}{"exonSeries"} .= $cols[3] . ".." . $cols[4] . ",";

	}
}
close FF;

# print "Finish load gtf into hash.\n";

# load cDNA sequence into hash
# print "Begin load cDNA sequence into hash.\n";
my ($isoformId, $isoformSymbol, $geneId, $geneSymbol, $proteinId);
open FF, "<$inputCdnaSeqFile";
# >rna65341 transcript_name:XM_024984141.1 gene_id:ENSBTAG00000049619 gene_name:GFOD1
while(my $line =<FF>){
	chomp($line);
	if($line=~/>(.*) transcript_name:(.*) gene_id:(.*) gene_name:(.*)/){
		$isoformId = $1;
	}else{
		${$isoform{$isoformId}}{"cDNAseq"}.=$line;
	}
}
close FF;

# print "Finish load cDNA sequence into hash.\n";


# print "Begin load DNA sequence into hash.\n";
open FF, "<$inputDNAseqFile";
# >rna65341 transcript_name:XM_024984141.1 gene_id:ENSBTAG00000049619 gene_name:GFOD1
while(my $line =<FF>){
	chomp($line);
	if($line=~/>(.*) transcript_name:(.*) gene_id:(.*) gene_name:(.*)/){
		$isoformId = $1;
	}else{
		${$isoform{$isoformId}}{"DNAseq"}.=$line;
	}
}
close FF;

# print "Finish load DNA sequence into hash.\n";

# print "Begin load pep sequence into hash.\n";

open FF, "<$inputPepSeqFile";
# >rna65341 transcript_name:XM_024984141.1 gene_id:ENSBTAG00000049619 gene_name:GFOD1 protein_id:XP_024839909.1
while(my $line =<FF>){
	chomp($line);
	if($line=~/>(.*) transcript_name:(.*) gene_id:(.*) gene_name:(.*) protein_id:(.*)/){
		$isoformId = $1;
		${$isoform{$isoformId}}{"proteinId"} = $5;
	}else{
		${$isoform{$isoformId}}{"pepSeq"}.=$line;
	}
}
close FF;

# print "Finish load pep sequence into hash.\n";


# print "Begin load Pfam and Go annotation into hash.\n";
open FF, "<$inputPfamGoAnnoFile";
# ENSBTAT00000073238      PF12130,PF00412,PF00307,PF01494 GO:0005515,GO:0071949
while(my $line=<FF>){
	chomp($line);
	my @tt = ();
	@tt = split(/\t/, $line);
	${$isoform{$tt[0]}}{"Pfam"} = $tt[1];
	${$isoform{$tt[0]}}{"Go"} = $tt[2];
}
close FF;

# print "Finish load Pfam and Go annotation into hash.\n";

##################################################################
# generate DNA seq for isoform and distinguish intron and exon

# print "Begin isoform into mysql.\n";
my (@isoformId, @exonIntron, @tmpExon, @exon, $exon, $query);
@isoformId = keys(%isoform);
foreach $isoformId(@isoformId){

	# process exonSeries to cut final comma
	${$isoform{$isoformId}}{"exonSeries"} = substr(${$isoform{$isoformId}}{"exonSeries"}, 0, length(${$isoform{$isoformId}}{"exonSeries"})-1);

	# load isoform into mysql database
	my $sql ="insert into isoformTable (isoformId, species, isoformSymbol, geneId, geneSymbol, chr, strand, start, end, exonSeries, annoOrigin, cDNAseq, DNAseq, pepSeq, isoformPfam, isoformGo, proteinId) values (\"$isoformId\", \"$species\", \"" . ${$isoform{$isoformId}}{"isoformSymbol"} . "\", \"" . ${$isoform{$isoformId}}{"geneId"} . "\", \"" . ${$isoform{$isoformId}}{"geneSymbol"} . "\", \"" . ${$isoform{$isoformId}}{"chr"} . "\", \"" . ${$isoform{$isoformId}}{"strand"} . "\", " . ${$isoform{$isoformId}}{"start"} . ", " . ${$isoform{$isoformId}}{"end"} . ", \"" . ${$isoform{$isoformId}}{"exonSeries"} . "\", \"" . ${$isoform{$isoformId}}{"annoOrigin"} . "\", \"". ${$isoform{$isoformId}}{"cDNAseq"} . "\", \"" . ${$isoform{$isoformId}}{"DNAseq"} . "\", \"" . ${$isoform{$isoformId}}{"pepSeq"} . "\", \"" . ${$isoform{$isoformId}}{"Pfam"} . "\", \"" . ${$isoform{$isoformId}}{"Go"} . "\", \"" . ${$isoform{$isoformId}}{"proteinId"}  . "\")";
	#print $sql;
	#<STDIN>;
	$query = $dbh->prepare($sql);
	$query->execute();

}

# print "Finish isoform into mysql.\n";

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
