#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--origTranscriptome \\\n" .
		"--newTranscriptome \n";
	exit;
}

my ($origTranscriptome, $newTranscriptome);
my ($line, @field, @attr, $attrString, $attr);
GetOptions(
        'origTranscriptome=s'=>\$origTranscriptome,
        'newTranscriptome=s'=>\$newTranscriptome,
);

my ($refGeneId, $refTrsptId, $cov, $tpm, $fpkm, $exonNumber);
open FF, "<$origTranscriptome";
open WW, ">$newTranscriptome";
while($line=<FF>){
	chomp($line);
	next if($line=~/^#/);
	@field = split(/\t/, $line);

	$attrString = pop(@field);
	$attrString = substr($attrString, 0, length($attrString) - 1) if(substr($attrString, length($attrString) - 1, 1) eq ";");

	($refGeneId, $refTrsptId, $cov, $tpm, $fpkm, $exonNumber) = ("", "", "", "", "", "");
	&getRefId($attrString, \$refGeneId, \$refTrsptId, \$cov, \$tpm, \$fpkm, \$exonNumber);

	# 检测发现当前转录本（及其exon）来自参考注释，则输出，否则放弃
	if($refGeneId ne ""){
		if($field[2] eq "transcript"){
			$attrString = "gene_id \"$refGeneId\"; transcript_id \"$refTrsptId\"; cov \"$cov\"; FPKM \"$fpkm\"; TPM \"$tpm\";";
		}elsif($field[2] eq "exon"){
			$attrString = "gene_id \"$refGeneId\"; transcript_id \"$refTrsptId\"; exon_number \"$exonNumber\"; cov \"$cov\";";
		}
		print WW join("\t", @field, $attrString) . "\n";
	}

}
close FF;
close WW;

sub getRefId{
	my ($attring, $refGeneId, $refTrsptId, $cov, $tpm, $fpkm, $exonNumber) = @_;
	#transcript: gene_id "ERX2120647.24900"; transcript_id "ERX2120647.24900.1"; reference_id "SRX826605.29097.1"; ref_gene_id "Os11g0116200"; cov "230.583450"; FPKM "58.123074"; TPM "87.043785";
	#exon: gene_id "ERX2120647.24900"; transcript_id "ERX2120647.24900.4"; exon_number "1"; reference_id "Os11t0116200-02"; ref_gene_id "Os11g0116200"; ref_gene_name "OsLTP1.13"; cov "2.901192";
	($$refGeneId, $$refTrsptId, $$cov, $$tpm, $$fpkm, $$exonNumber) = ("", "", "", "", "", "");
	my (@attr, $attr);
	@attr = split(/; /, $attring);
	foreach $attr(@attr){
		if($attr =~/ref_gene_id "(.*)"/){
			$$refGeneId = $1;
		}elsif($attr =~/reference_id "(.*)"/){
			$$refTrsptId = $1;
		}elsif($attr =~/FPKM "(.*)"/){
			$$fpkm = $1;
		}elsif($attr =~/TPM "(.*)"/){
			$$tpm = $1;
		}elsif($attr =~/cov "(.*)"/){
			$$cov = $1;
		}elsif($attr =~/exon_number "(.*)"/){
			$$exonNumber = $1;
		}
	}
}
