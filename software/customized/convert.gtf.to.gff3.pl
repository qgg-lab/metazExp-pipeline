#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputGtf \\\n" .
                "--outputGff \n";
	exit;
}

my ($inputGtf, $outputGff);

GetOptions(
        'inputGtf=s'=>\$inputGtf,
        'outputGff=s'=>\$outputGff,
);

# 以gene单位将gtf读入
# 1       araport11       transcript      3631    5899    .       +       .       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; gene_name "NAC001"; gene_source "araport11"; gene_biotype "protein_coding"; transcript_source "araport11"; transcript_biotype "protein_coding";
# 1       araport11       exon    3631    3913    .       +       .       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; exon_number "1"; gene_name "NAC001"; gene_source "araport11"; gene_biotype "protein_coding"; transcript_source "araport11"; transcript_biotype "protein_coding"; exon_id "AT1G01010.1.exon1";
my ($gene_id, $gene_name, $gene_source, $gene_biotype, $transcript_id, $transcript_name, $transcript_source, $transcript_biotype, $exon_id, $exon_number);
my (@field, $attrString);
my (%geneGtf, $geneGtfHref);
$geneGtfHref = \%geneGtf;
my ($line);
open FF, "<$inputGtf";
while($line=<FF>){
	chomp($line);
	($gene_id, $gene_name, $gene_source, $gene_biotype, $transcript_id, $transcript_name, $transcript_source, $transcript_biotype, $exon_id, $exon_number) = ("", "", "", "", "", "", "", "", "", "");
	@field = split(/\t/, $line);
	&getAttrValue($field[8], \$gene_id, \$gene_name, \$gene_source, \$gene_biotype, \$transcript_id, \$transcript_name, \$transcript_source, \$transcript_biotype, \$exon_id, \$exon_number);
	if($field[1] eq "StringTie"){
		$field[1] = "ASplant";
	}

	# 登记geneFeature
	$geneGtfHref->{$gene_id}->{"chr"} = $field[0];
	$geneGtfHref->{$gene_id}->{"source"} = $field[1];
	$geneGtfHref->{$gene_id}->{"strand"} = $field[6];
	if($gene_name ne ""){
		$geneGtfHref->{$gene_id}->{"gene_name"} = $gene_name;
	}
	if($gene_source ne ""){
		$geneGtfHref->{$gene_id}->{"gene_source"} = $gene_source;
	}
	if($gene_biotype ne ""){
		$geneGtfHref->{$gene_id}->{"gene_biotype"} = $gene_biotype;
	}
	if(not(exists($geneGtfHref->{$gene_id}->{"start"}))){
		$geneGtfHref->{$gene_id}->{"start"} = $field[3];
	}elsif($field[3] < $geneGtfHref->{$gene_id}->{"start"}){
		$geneGtfHref->{$gene_id}->{"start"} = $field[3];
	}

	if(not(exists($geneGtfHref->{$gene_id}->{"stop"}))){
		$geneGtfHref->{$gene_id}->{"stop"} = $field[4];
	}elsif($field[4] > $geneGtfHref->{$gene_id}->{"stop"}){
		$geneGtfHref->{$gene_id}->{"stop"} = $field[4];
	}

	# 登记
	if($field[2] eq "transcript"){
		$field[2] = "mRNA";
		pop(@field);
		$attrString = "ID=" . $transcript_id . ";Parent=" . $gene_id;
		if($transcript_name ne ""){
			$attrString .= ";Name=" . $transcript_name;
		}
		if($transcript_source ne ""){
			$attrString .= ";Source=" . $transcript_source;
		}
		if($transcript_biotype ne ""){
			$attrString .= ";Biotype=" . $transcript_biotype;
		}
		$geneGtfHref->{$gene_id}->{"transcript"}->{$transcript_id}->{"mRNA"} = join("\t", @field, $attrString);
	}
	if($field[2] eq "exon"){
		pop(@field);
		if($exon_id eq ""){
			$attrString = "ID=exon." . $transcript_id . "." . $exon_number . ";Parent=" . $transcript_id;
		}else{
			$attrString = "ID=" . $exon_id . ";Parent=" . $transcript_id;
		}
		$geneGtfHref->{$gene_id}->{"transcript"}->{$transcript_id}->{"exonText"} .= join("\t", @field, $attrString) . "\n";
	}
}
close FF;

# 按照染色体位置对gene进行排序
my (@transcript_id);
my @gene_id =keys(%geneGtf);
open WW, ">$outputGff";
foreach $gene_id(@gene_id){
	print WW join("\t", $gene_id, $geneGtfHref->{$gene_id}->{"chr"}, $geneGtfHref->{$gene_id}->{"start"}) . "\n";
}
close WW;
my $srtGeneIdList =`sort -k2,2 -k3,3n $outputGff | awk -F \'\\t\' \'{print \$1}\'`;
@gene_id = split(/\n/, $srtGeneIdList);
open WW, ">$outputGff";
foreach $gene_id(@gene_id){
	$attrString = "";
	$attrString = "ID=" . $gene_id;
	if(exists($geneGtfHref->{$gene_id}->{"gene_name"})){
		$attrString .= ";Name=" . $geneGtfHref->{$gene_id}->{"gene_name"};
	}
	if(exists($geneGtfHref->{$gene_id}->{"gene_source"})){
		$attrString .= ";Source=" . $geneGtfHref->{$gene_id}->{"gene_source"};
	}
	if(exists($geneGtfHref->{$gene_id}->{"gene_biotype"})){
		$attrString .= ";Biotype=" . $geneGtfHref->{$gene_id}->{"gene_biotype"};
	}
	# 输出gene feature
	print WW join("\t", $geneGtfHref->{$gene_id}->{"chr"}, $geneGtfHref->{$gene_id}->{"source"}, "gene", $geneGtfHref->{$gene_id}->{"start"}, $geneGtfHref->{$gene_id}->{"stop"}, ".", $geneGtfHref->{$gene_id}->{"strand"}, ".", $attrString) . "\n";

	# 输出transcript
	@transcript_id = ();
	@transcript_id = keys(%{$geneGtfHref->{$gene_id}->{"transcript"}});
	# gene_name gene_source gene_biotype start stop chr source strand
	foreach $transcript_id(@transcript_id){
		# 输出transcript feature
		print WW $geneGtfHref->{$gene_id}->{"transcript"}->{$transcript_id}->{"mRNA"} . "\n";
		# 输出exon feature
		print WW $geneGtfHref->{$gene_id}->{"transcript"}->{$transcript_id}->{"exonText"};
	}
}
close WW;


sub getAttrValue{
	my ($attrString, $gene_id, $gene_name, $gene_source, $gene_biotype, $transcript_id, $transcript_name, $transcript_source, $transcript_biotype, $exon_id, $exon_number) = @_;
	my (@attr, $attr);
	@attr = split(/; /, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$gene_id = $1;
		}elsif($attr=~/gene_name "(.*)"/){
			$$gene_name = $1;
		}elsif($attr=~/gene_source "(.*)"/){
			$$gene_source = $1;
		}elsif($attr=~/gene_biotype "(.*)"/){
			$$gene_biotype = $1;
		}elsif($attr=~/transcript_id "(.*)"/){
			$$transcript_id = $1;
		}elsif($attr=~/transcript_name "(.*)"/){
			$$transcript_name = $1;
		}elsif($attr=~/transcript_source "(.*)"/){
			$$transcript_source = $1;
		}elsif($attr=~/transcript_biotype "(.*)"/){
			$$transcript_biotype = $1;
		}elsif($attr=~/exon_id "(.*)"/){
			$$exon_id = $1;
		}elsif($attr=~/exon_number "(.*)"/){
			$$exon_number = $1;
		}
	}
}

