#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--origGtf \\\n" .
                "--newGtf \n";
	exit;
}

my ($origGtf, $newGtf);

GetOptions(
        'origGtf=s'=>\$origGtf,
        'newGtf=s'=>\$newGtf,
);

# 从origGtf读入gene->trspt->text
# NC_029679.1     Gnomon  gene    8753    10205   .       -       .       gene_id "LOC107415205"; db_xref "GeneID:107415205"; gbkey "Gene"; gene "LOC107415205"; gene_biotype "protein_coding";
# NC_029679.1     Gnomon  exon    8753    10205   .       -       .       gene_id "LOC107415205"; transcript_id "XM_016023852.2"; db_xref "GeneID:107415205"; gbkey "mRNA"; gene "LOC107415205"; model_evidence "Supporting evidence includes similarity to: 1 Protein, and 100% coverage of the annotated genomic feature by RNAseq alignments"; product "uncharacterized LOC107415205"; exon_number "1"; 
my (%gtfHash, $gtfHref, $line, @field, @geneId, $geneId, $trsptId, @trsptId, %geneFeature);
$gtfHref = \%gtfHash;
open FF, "<$origGtf";
while($line=<FF>){
	@field = ();
	@field = split(/\t/, $line);
	next if($#field!=8);

	($geneId, $trsptId) = ("", "");
	&getGeneIdAndTrsptId($field[8], \$geneId, \$trsptId);
	if($field[2] eq "gene"){
		$geneFeature{$geneId} = $line;
	}else{
		if(not exists($gtfHref->{$geneId}->{$trsptId})){
			$gtfHref->{$geneId}->{$trsptId}->{"chr"} = $field[0];
			$gtfHref->{$geneId}->{$trsptId}->{"strand"} = $field[6];
			$gtfHref->{$geneId}->{$trsptId}->{"start"} = $field[3];
			$gtfHref->{$geneId}->{$trsptId}->{"stop"} = $field[4];
			$gtfHref->{$geneId}->{$trsptId}->{"text"}=$line;
		}else{
			$gtfHref->{$geneId}->{$trsptId}->{"start"} = $field[3] if($gtfHref->{$geneId}->{$trsptId}->{"start"} > $field[3]);
			$gtfHref->{$geneId}->{$trsptId}->{"stop"} = $field[4] if($gtfHref->{$geneId}->{$trsptId}->{"stop"} < $field[4]);
			$gtfHref->{$geneId}->{$trsptId}->{"text"} .= $line;
		}

	}
}
close FF;


# 输出所有gene
open WW, ">$newGtf";
@geneId = keys(%gtfHash);
foreach $geneId(@geneId){
	print WW $geneFeature{$geneId};
	@trsptId = ();
	@trsptId = keys(%{$gtfHref->{$geneId}});
	foreach $trsptId(@trsptId){
		print WW join("\t", $gtfHref->{$geneId}->{$trsptId}->{"chr"}, "Gnomon", "transcript", $gtfHref->{$geneId}->{$trsptId}->{"start"}, $gtfHref->{$geneId}->{$trsptId}->{"stop"}, ".", $gtfHref->{$geneId}->{$trsptId}->{"strand"}, ".", "gene_id \"$geneId\"; transcript_id \"$trsptId\"\n");
		print WW $gtfHref->{$geneId}->{$trsptId}->{"text"};
	}
}
close WW;

sub getGeneIdAndTrsptId{
	my ($attrString, $geneId, $trsptId) = @_;
	my (@attr, $attr);
	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$geneId = $1;
		}
		if($attr=~/transcript_id "(.*)"/){
			$$trsptId = $1;
		}
	}
}
