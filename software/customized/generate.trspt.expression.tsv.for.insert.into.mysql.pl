#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--finalTrsptIdListTsv final.trsptId.list.tsv \\\n" .
		"--inputAssemblyDir \\\n" . 
		"--filteredExptInfoTsv \\\n" .
		"--outputTrsptExpFile trspt.exp.mysql.tsv \n\n";
	exit;
}

my ($inputAssemblyDir, $filteredExptInfoTsv, $outputTrsptExpFile, $finalTrsptIdListTsv);

GetOptions(
	'finalTrsptIdListTsv=s'=>\$finalTrsptIdListTsv,
        'inputAssemblyDir=s'=>\$inputAssemblyDir,
        'filteredExptInfoTsv=s'=>\$filteredExptInfoTsv,
        'outputTrsptExpFile=s'=>\$outputTrsptExpFile,
);

my (@fieldName, @field, $i, $exptId, $exptIdPos, $geneExpFile);
my ($experimentLine, %experiment, $trsptLine, $geneId, $cov, $fpkm, $tpm, $trsptId, $tissue, $treatment);

# 转录本序列编号读入hash
my (%trsptId);
open FF, "<$finalTrsptIdListTsv";
# >SRX1660584.17147.2 transcript_name:NA gene_id:AT5G22830 gene_name:NA
# >AT5G45110.2 transcript_name:NA gene_id:AT5G45110 gene_name:NPR3 
while(my $line=<FF>){
	if($line=~/>(.*?) .*/){
		$trsptId{$1}=1;
	}
}
close FF;

open FF, "<$filteredExptInfoTsv";
open WW, ">$outputTrsptExpFile";
$experimentLine=<FF>;
chomp($experimentLine);

# 获得首行字段
@fieldName = split(/\t/, $experimentLine);
my $exptIdPos = 0;
for($exptIdPos=0; $exptIdPos<=$#fieldName; $exptIdPos++){
        last if($fieldName[$exptIdPos] eq "Experiment");
}

my ($cmd, $grepRltText, @trsptLine, $trsptLine, %tmpExptHash, $tmpExptHashHref);
$tmpExptHashHref=\%tmpExptHash;
# 依次读取每个实验的编号和对应的表达文件
while($experimentLine=<FF>){
        chomp($experimentLine);
        @field = ();
        @field = split(/\t/, $experimentLine);

	for(my $j=0; $j<=$#field; $j++){
		$tmpExptHashHref->{$fieldName[$j]} = $field[$j];
	}

	$tissue = $tmpExptHashHref->{"Tissue"};
	$tissue=~tr/ /_/;
	$treatment=$tmpExptHashHref->{"Treatment"};
	$treatment=~tr/ /_/;
        $exptId = $field[$exptIdPos];
	$geneExpFile = $inputAssemblyDir . "/" . $exptId . "/transcriptomeByStringtie.gtf";
        if(not(-e $geneExpFile)){
		print $exptId . " have no assembled transcriptome gtf file.\n";   
                next;
        }
	#print $exptId;
	#<STDIN>;
	# 读取表达水平文件，只保留geneId, exptId, coverage, fpkm, tpm
	# 注意： 表达transcriptome中提供的gene_id和transcript_id不是final中的编号的，而是重新又组装过的编号
	#        我们应该用reference_id(第3个)，ref_gene_id（第4个）

	# 1       StringTie       transcript      3631    5899    1000    +       .       gene_id "SRX1426062.1"; transcript_id "SRX1426062.1.1"; reference_id "AT1G01010.1"; ref_gene_id "AT1G01010"; ref_gene_name "NAC001"; cov "51.289845"; FPKM "5.365530"; TPM "8.648281";
	$cmd = "grep -P \"\\ttranscript\\t.*.*reference_id .* ref_gene_id\" " . $geneExpFile . " | awk -F \'\\t\' \'{print \$9}\'";

	# gene_id "SRX1426062.1"; transcript_id "SRX1426062.1.1"; reference_id "AT1G01010.1"; ref_gene_id "AT1G01010"; ref_gene_name "NAC001"; cov "51.289845"; FPKM "5.365530"; TPM "8.648281";
	$grepRltText = `$cmd`;
	@trsptLine = ();
	@trsptLine = split(/\n/, $grepRltText);
	foreach $trsptLine(@trsptLine){
		&getGeneTrsptIdAndExpr($trsptLine, \$geneId, \$trsptId, \$cov, \$fpkm, \$tpm);
		# 在finalGtf中不存在这个转录本
		next if(not exists($trsptId{$trsptId}));
		# 直接将treatment设置为NA
		$treatment = "NA";
		print WW "GeneId, TrsptId, ExptId, Tissue, Treatment, Cov, FPKM, TPM___" . $geneId . ", " . $trsptId . ", " .  $exptId . ", " . $tissue . ", " . $treatment . ", " . sprintf("%.3f", $cov) . ", " . sprintf("%.3f", $fpkm) . ", " . sprintf("%.3f", $tpm) . "\n";
	}
}
close WW;

sub getGeneTrsptIdAndExpr{
	my ($attrString, $geneId, $trsptId, $cov, $fpkm, $tpm)=@_;
	my (@field, $field);
	@field = split(/; /, $attrString);
	foreach $field(@field){
		if($field=~/reference_id "(.*)"/){
			$$trsptId = $1;
		}
		if($field=~/ref_gene_id "(.*)"/){
			$$geneId = $1;
		}
		if($field=~/cov "(.*)"/){
			$$cov = $1;
		}
		if($field=~/FPKM "(.*)"/){
			$$fpkm = $1;
		}
		if($field=~/TPM "(.*)"/){
			$$tpm = $1;
		}
	}
}
