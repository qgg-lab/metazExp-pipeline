#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw/sum/;
use List::Util qw/max min/;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputTrsptGtf \\\n" .
                "--inputAsFile \\\n" .
                "--inputAsPsiFile \\\n" .
                "--inputExperimentFile \\\n" .
		"--outputAsGff3 \n";
	exit;
}

my ($inputAsFile, $inputAsPsiFile, $inputExperimentFile, $outputAsGff3, $inputTrsptGtf);

GetOptions(
        'inputTrsptGtf=s'=>\$inputTrsptGtf,
        'inputAsFile=s'=>\$inputAsFile,
        'inputAsPsiFile=s'=>\$inputAsPsiFile,
        'inputExperimentFile=s'=>\$inputExperimentFile,
        'outputAsGff3=s'=>\$outputAsGff3,
);

# 将基因的span读入
my (@field, $geneId, $trsptId, %geneSpan, $geneSpanHref);
$geneSpanHref = \%geneSpan;
open FF, "<$inputTrsptGtf";
while(my $line = <FF>){
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "transcript");
	($geneId, $trsptId) = ("", "");
	&getGeneIdAndTrsptId($field[8], \$geneId, \$trsptId);

#	print "geneId:" . $geneId . "," . "trsptId:". $trsptId;
#	<STDIN>;
	if(not exists($geneSpanHref->{$geneId})){
		$geneSpanHref->{$geneId}->{"chr"} = $field[0];
		$geneSpanHref->{$geneId}->{"strand"} = $field[6];
		$geneSpanHref->{$geneId}->{"start"} = $field[3];
		$geneSpanHref->{$geneId}->{"stop"} = $field[4];
	}else{
		$geneSpanHref->{$geneId}->{"start"} = $field[3] if($field[3]<$geneSpanHref->{$geneId}->{"start"});
		$geneSpanHref->{$geneId}->{"stop"} = $field[4] if($field[4] > $geneSpanHref->{$geneId}->{"stop"});
	}
}
close FF;
print "Finish load gene span\n";
my @geneId = keys(%geneSpan);

# 将experiment信息读入，获得experiment->tissue的hash
my (%exptToTissue, $exptToTissueHref, @tt, $exptId, $tissue);
my ($fieldNameString, @fieldName, $valueString, @value, $i, %tmpHash, $tmpHash, $line);
$exptToTissueHref = \%exptToTissue;
$tmpHash = \%tmpHash;
open FF, "<$inputExperimentFile";
while($line=<FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);
	@fieldName = ();
	@fieldName = split(/, /, $fieldNameString);
	@value = ();
	@value = split(/, /, $valueString);
	%tmpHash = ();
	for($i=0; $i<=$#fieldName; $i++){
		$tmpHash->{$fieldName[$i]} = $value[$i];
	}
	@tt = ();
	@tt = split(/"/, $tmpHash->{"experimentId"});
	$exptId = $tt[1];
	@tt = ();
	@tt = split(/"/, $tmpHash->{"tissue"});
	$tissue = $tt[1];
	$tissue=~s/ /_/;
	$exptToTissueHref->{$exptId} = $tissue;
}
close FF;

print "Finish load experiment information.\n";

# 将psi值按照as->tissue = psiList
my (%asToTissueToPsiList, $asToTissueToPsiList, $asId,  %geneToAsGff3, $geneToAsGff3Href);
$geneToAsGff3Href=\%geneToAsGff3;
$asToTissueToPsiList = \%asToTissueToPsiList;
open FF, "<$inputAsPsiFile";
while($line=<FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);
        @fieldName = ();
        @fieldName = split(/, /, $fieldNameString);
        @value = ();
        @value = split(/, /, $valueString);
        %tmpHash = ();
        for($i=0; $i<=$#fieldName; $i++){
                $tmpHash->{$fieldName[$i]} = $value[$i];
        }
	@tt = ();
#	print "asId: " .$tmpHash->{"asId"};
#	#<STDIN>;
	@tt= split(/"/, $tmpHash->{"asId"});
	$asId = $tt[1];
#	print "experiment: " . $tmpHash->{"experiment"};
#	#<STDIN>;
	@tt = ();
	@tt = split(/"/, $tmpHash->{"experiment"});
	$exptId = $tt[1];
#	print join("\t", "asId:", $asId, "experiment:", $exptId, "tissue:", $exptToTissueHref->{$exptId}, "psi:", $tmpHash->{"JCECpsi"});
#	#<STDIN>;
	if(not exists($asToTissueToPsiList->{$asId}->{$exptToTissueHref->{$exptId}})){
		$asToTissueToPsiList->{$asId}->{$exptToTissueHref->{$exptId}} =  $tmpHash->{"JCECpsi"};
	}else{
		$asToTissueToPsiList->{$asId}->{$exptToTissueHref->{$exptId}} .= "," . $tmpHash->{"JCECpsi"};
	}
}

print "Finish load psi information.\n";

# 读取as，并且提取每个 tissue，计算每个tissue的psi最大值，最小值和平均值
# 生成外显子结构
my (%tissueToPsi, $tissueToPsi, $asType, $strand, $chr, @tissue, $tissue, @psi, $gtfTissuePsiAttrs, $gffTissuePsiAttrs, $asStart, $asStop);
open FF, "<$inputAsFile";	
while($line=<FF>){
#	print $line;
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);
        @fieldName = ();
        @fieldName = split(/, /, $fieldNameString);
        @value = ();
        @value = split(/, /, $valueString);
        %tmpHash = ();
        for($i=0; $i<=$#fieldName; $i++){
		if($fieldName[$i]=~/"(.*)"/){
			$fieldName[$i] = $1;
		}
		if($value[$i]=~/"(.*)"/){
			$value[$i] = $1;
		}
                $tmpHash->{$fieldName[$i]} = $value[$i];
        }
	%tissueToPsi = ();
	$tissueToPsi = \%tissueToPsi;
	$asId = $tmpHash->{"asId"};
	$asType = $tmpHash->{"asType"};
	$chr = $tmpHash->{"chr"};
	$strand = $tmpHash->{"strand"};
#	print "asId:" . $asId . "\n";
#	print "asType:" . $asType . "\n";
#	print "chr:" . $chr . "\n";
#	print "strand:" . $strand . "\n";
	#<STDIN>;
	# 计算AS被检测到psi值在各个tissue中最大值，最小值，平均值情况
	$gffTissuePsiAttrs = "";
	$gtfTissuePsiAttrs = "";
	@tissue = ();
	if(exists($asToTissueToPsiList->{$asId})){
		@tissue = keys($asToTissueToPsiList->{$asId});
		my $tissueNum =0;
		foreach $tissue(@tissue){
			$tissueNum++;
			@psi = ();
			@psi = split(/,/, $asToTissueToPsiList->{$asId}->{$tissue});
#			print "@psi";
#			print "\n";
			$tissueToPsi->{$tissue}->{"maxPsi"} = max @psi;
			$tissueToPsi->{$tissue}->{"minPsi"} = min @psi;
			$tissueToPsi->{$tissue}->{"meanPsi"} = sum(@psi)/($#psi + 1);
			$tissueToPsi->{$tissue}->{"exptNum"} = $#psi + 1;
#			print "$tissue max:" . $tissueToPsi->{$tissue}->{"maxPsi"};
			#<STDIN>;
#			print "$tissue min:" . $tissueToPsi->{$tissue}->{"minPsi"};
			#<STDIN>;
#			print $tissue . "mean:" .$tissueToPsi->{$tissue}->{"meanPsi"};
			#<STDIN>;
			if($gffTissuePsiAttrs eq ""){
				$gffTissuePsiAttrs = "PSI_" . $tissue . "(" . $tissueToPsi->{$tissue}->{"exptNum"} . "exp.)=" . "max:" . sprintf("%.2f", $tissueToPsi->{$tissue}->{"maxPsi"}) . ",min:" . sprintf("%.2f", $tissueToPsi->{$tissue}->{"minPsi"}) . ",mean:" . sprintf("%.2f", $tissueToPsi->{$tissue}->{"meanPsi"});
			}else{
				$gffTissuePsiAttrs .= ";" . "PSI_" . $tissue . "(" . $tissueToPsi->{$tissue}->{"exptNum"} . "exp.)=" . "max:" . sprintf("%.2f", $tissueToPsi->{$tissue}->{"maxPsi"}) . ",min:" . sprintf("%.2f", $tissueToPsi->{$tissue}->{"minPsi"}) . ",mean:" . sprintf("%.2f", $tissueToPsi->{$tissue}->{"meanPsi"});
			}
		}
	}else{
		print $asId . " has no psi.\n";
#		<STDIN>;
	}
#	print "gffTissuePsiAttrs:$gffTissuePsiAttrs\n";
#	print "gtfTissuePsiAttrs:$gtfTissuePsiAttrs\n";
	#<STDIN>;
	# 生成AS的gtf格式
	$asStart = 0;
	$asStop = 0;


#################### A3SS #############################
	if($asType eq "A3SS"){
		$asStart = min(($tmpHash->{"longExonStart_0base"}+1, $tmpHash->{"longExonEnd"}, $tmpHash->{"shortES"}+1, $tmpHash->{"shortEE"}, $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}));
		$asStop = max(($tmpHash->{"longExonStart_0base"}+1, $tmpHash->{"longExonEnd"}, $tmpHash->{"shortES"}+1, $tmpHash->{"shortEE"}, $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}));


################ GFF3格式  ###############
	######### 1: 输出AS事件  #########
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "gene", $asStart, $asStop, ".", $strand, ".", "ID=$asId;AsType=$asType;GeneId=" . $tmpHash->{"geneID"} . ";$gffTissuePsiAttrs\n");

	######### 2.1: 输出longAlt #######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.long" . ";Parent=" . $asId . ";Alt=long;\n");
		# 2.2: 输出long exon
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}, ".", $strand, ".", "ID=$asId.long.flankingExon" . ";Parent=$asId.long;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}, ".", $strand, ".", "ID=CDS:$asId.long.flankingExon" . ";Parent=$asId.long;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"longExonStart_0base"}+1, $tmpHash->{"longExonEnd"}, ".", $strand, ".", "ID=$asId.longExon" . ";Parent=$asId.long;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"longExonStart_0base"}+1, $tmpHash->{"longExonEnd"}, ".", $strand, ".", "ID=CDS:$asId.longExon" . ";Parent=$asId.long;\n");

	######### 3.1: 输出shortAlt ######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.short" . ";Parent=" . $asId . ";Alt=short;\n");
		# 3.2: 输出shortexon
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}, ".", $strand, ".", "ID=$asId.short.flankingExon" . ";Parent=$asId.short;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}, ".", $strand, ".", "ID=CDS:$asId.short.flankingExon" . ";Parent=$asId.short;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"shortES"}+1, $tmpHash->{"shortEE"}, ".", $strand, ".", "ID=$asId.shortExon" . ";Parent=$asId.short;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"shortES"}+1, $tmpHash->{"shortEE"}, ".", $strand, ".", "ID=CDS:$asId.shortExon" . ";Parent=$asId.short;\n");
			
	}


##################  A5SS #############################################
	if($asType eq "A5SS"){
		$asStart = min(($tmpHash->{"longExonStart_0base"}+1, $tmpHash->{"longExonEnd"}, $tmpHash->{"shortES"}+1, $tmpHash->{"shortEE"}, $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}));
		$asStop = max(($tmpHash->{"longExonStart_0base"}+1, $tmpHash->{"longExonEnd"}, $tmpHash->{"shortES"}+1, $tmpHash->{"shortEE"}, $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}));


################ GFF3格式  ###############
	######### 1: 输出AS事件  #########
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "gene", $asStart, $asStop, ".", $strand, ".", "ID=$asId;AsType=$asType;geneId=" . $tmpHash->{"geneID"} . ";$gffTissuePsiAttrs\n");

	######### 2.1: 输出longAlt #######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.long" . ";Parent=" . $asId . ";Alt=long;\n");
		# 2.2: 输出long exon
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"longExonStart_0base"}+1, $tmpHash->{"longExonEnd"}, ".", $strand, ".", "ID=$asId.longExon" . ";Parent=$asId.long;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .=join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"longExonStart_0base"}+1, $tmpHash->{"longExonEnd"}, ".", $strand, ".", "ID=CDS:$asId.longExon" . ";Parent=$asId.long;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}, ".", $strand, ".", "ID=$asId.long.flankingExon" . ";Parent=$asId.long;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}, ".", $strand, ".", "ID=CDS:$asId.long.flankingExon" . ";Parent=$asId.long;\n");

	######### 3.1: 输出shortAlt ######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.short" . ";Parent=" . $asId . ";Alt=short;\n");
		# 3.2: 输出shortexon
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"shortES"}+1, $tmpHash->{"shortEE"}, ".", $strand, ".", "ID=$asId.shortExon" . ";Parent=$asId.short;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"shortES"}+1, $tmpHash->{"shortEE"}, ".", $strand, ".", "ID=CDS:$asId.shortExon" . ";Parent=$asId.short;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}, ".", $strand, ".", "ID=$asId.short.flankingExon" . ";Parent=$asId.short;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"flankingES"}+1, $tmpHash->{"flankingEE"}, ".", $strand, ".", "ID=CDS:$asId.short.flankingExon" . ";Parent=$asId.short;\n");


	}



###############  RI #################################
	if($asType eq "RI"){
		$asStart = $tmpHash->{"riExonStart_0base"}+1;
		$asStop = $tmpHash->{"riExonEnd"};

################ GFF3格式  ###############
	######### 1: 输出AS事件  #########
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "gene", $asStart, $asStop, ".", $strand, ".", "ID=$asId;AsType=$asType;geneId=" . $tmpHash->{"geneID"} . ";$gffTissuePsiAttrs\n");

	######### 2.1: 输出includedIntron #######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.includedIntron" . ";Parent=" . $asId . ";Alt=includedIntron;\n");
		# 2.2: 输出exon
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"riExonStart_0base"}+1, $tmpHash->{"riExonEnd"}, ".", $strand, ".", "ID=$asId.includedIntron.riExon" . ";Parent=$asId.includedIntron;\n");
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"riExonStart_0base"}+1, $tmpHash->{"riExonEnd"}, ".", $strand, ".", "ID=CDS:$asId.includedIntron.riExon" . ";Parent=$asId.includedIntron;\n");

	######### 3.1: 输出excludedIntron ######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.excludedIntron" . ";Parent=" . $asId . ";Alt=excludedIntron;\n");
		# 3.2: 输出exon
		if($strand eq "+"){
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.excludedIntron.upstreamExon" . ";Parent=$asId.excludedIntron;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.excludedIntron.upstreamExon" . ";Parent=$asId.excludedIntron;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.excludedIntron.downstreamExon" . ";Parent=$asId.excludedIntron;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.excludedIntron.downstreamExon" . ";Parent=$asId.excludedIntron;\n");
		}else{
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.excludedIntron.downstreamExon" . ";Parent=$asId.excludedIntron;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.excludedIntron.downstreamExon" . ";Parent=$asId.excludedIntron;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.excludedIntron.upstreamExon" . ";Parent=$asId.excludedIntron;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.excludedIntron.upstreamExon" . ";Parent=$asId.excludedIntron;\n");
		}	

	}


############### SE  ###################
	if($asType eq "SE"){
		$asStart = $tmpHash->{"upstreamES"} + 1;
		$asStop = $tmpHash->{"downstreamEE"} + 1;

################ GFF3格式  ###############
	######### 1: 输出AS事件  #########
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "gene", $asStart, $asStop, ".", $strand, ".", "ID=$asId;AsType=$asType;geneId=" . $tmpHash->{"geneID"} . ";$gffTissuePsiAttrs\n");

	######### 2.1: 输出inclusion Alt #######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.inclusion" . ";Parent=" . $asId . ";Alt=inclusion;\n");
		# 2.2: 输出exon
		if($strand eq "+"){
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.inclusion.upstreamExon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.inclusion.upstreamExon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"exonStart_0base"}+1, $tmpHash->{"exonEnd"}, ".", $strand, ".", "ID=$asId.inclusion.exon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"exonStart_0base"}+1, $tmpHash->{"exonEnd"}, ".", $strand, ".", "ID=CDS:$asId.inclusion.exon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.inclusion.downstreamExon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.inclusion.downstreamExon" . ";Parent=$asId.inclusion;\n");
		}else{
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.inclusion.downstreamExon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.inclusion.downstreamExon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"exonStart_0base"}+1, $tmpHash->{"exonEnd"}, ".", $strand, ".", "ID=$asId.inclusion.exon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"exonStart_0base"}+1, $tmpHash->{"exonEnd"}, ".", $strand, ".", "ID=CDS:$asId.inclusion.exon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.inclusion.upstreamExon" . ";Parent=$asId.inclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.inclusion.upstreamExon" . ";Parent=$asId.inclusion;\n");
		}

	######### 3.1: 输出exclusionAlt ######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.exclusion" . ";Parent=" . $asId . ";Alt=exclusion;\n");
		# 3.2: 输出exon
		if($strand eq "+"){
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.exclusion.upstreamExon" . ";Parent=$asId.exclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.exclusion.upstreamExon" . ";Parent=$asId.exclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.exclusion.downstreamExon" . ";Parent=$asId.exclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.exclusion.downstreamExon" . ";Parent=$asId.exclusion;\n");
		}else{
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.exclusion.downstreamExon" . ";Parent=$asId.exclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.exclusion.downstreamExon" . ";Parent=$asId.exclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.exclusion.upstreamExon" . ";Parent=$asId.exclusion;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.exclusion.upstreamExon" . ";Parent=$asId.exclusion;\n");

		}
			

	}




############### MXE  ###################
	if($asType eq "MXE"){
		$asStart = $tmpHash->{"upstreamES"} + 1;
		$asStop = $tmpHash->{"downstreamEE"} + 1;

################ GFF3格式  ###############
	######### 1: 输出AS事件  #########
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "gene", $asStart, $asStop, ".", $strand, ".", "ID=$asId;AsType=$asType;geneId=" . $tmpHash->{"geneID"} . ";$gffTissuePsiAttrs\n");

	######### 2.1: 输出1st Alt #######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.1st" . ";Parent=" . $asId . ";Alt=1st;\n");
		# 2.2: 输出exon
		if($strand eq "+"){
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.1st.upstreamExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.1st.upstreamExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"1stExonStart_0base"}+1, $tmpHash->{"1stExonEnd"}, ".", $strand, ".", "ID=$asId.1st.1stExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"1stExonStart_0base"}+1, $tmpHash->{"1stExonEnd"}, ".", $strand, ".", "ID=CDS:$asId.1st.1stExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.1st.downstreamExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.1st.downstreamExon" . ";Parent=$asId.1st;\n");
		}else{
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.1st.downstreamExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.1st.downstreamExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"1stExonStart_0base"}+1, $tmpHash->{"1stExonEnd"}, ".", $strand, ".", "ID=$asId.1st.1stExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"1stExonStart_0base"}+1, $tmpHash->{"1stExonEnd"}, ".", $strand, ".", "ID=CDS:$asId.1st.1stExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.1st.upstreamExon" . ";Parent=$asId.1st;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.1st.upstreamExon" . ";Parent=$asId.1st;\n");
		}

	######### 3.1: 输出exclusionAlt ######
		$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "mRNA", $asStart, $asStop, ".", $strand, ".", "ID=$asId.2nd" . ";Parent=" . $asId . ";Alt=2nd;\n");
		# 3.2: 输出exon
		if($strand eq "+"){
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.2nd.upstreamExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.2nd.upstreamExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"2ndExonStart_0base"}+1, $tmpHash->{"2ndExonEnd"}, ".", $strand, ".", "ID=$asId.2nd.2ndExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"2ndExonStart_0base"}+1, $tmpHash->{"2ndExonEnd"}, ".", $strand, ".", "ID=CDS:$asId.2nd.2ndExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.2nd.downstreamExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.2nd.downstreamExon" . ";Parent=$asId.2nd;\n");
		}else{
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=$asId.2nd.downstreamExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"downstreamES"}+1, $tmpHash->{"downstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.2nd.downstreamExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"2ndExonStart_0base"}+1, $tmpHash->{"2ndExonEnd"}, ".", $strand, ".", "ID=$asId.2nd.2ndExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"2ndExonStart_0base"}+1, $tmpHash->{"2ndExonEnd"}, ".", $strand, ".", "ID=CDS:$asId.2nd.2ndExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "exon", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=$asId.2nd.upstreamExon" . ";Parent=$asId.2nd;\n");
			$geneToAsGff3Href->{$tmpHash->{"geneID"}} .= join("\t", $tmpHash->{"chr"}, "aslive", "CDS", $tmpHash->{"upstreamES"}+1, $tmpHash->{"upstreamEE"}, ".", $strand, ".", "ID=CDS:$asId.2nd.upstreamExon" . ";Parent=$asId.2nd;\n");

		}
			

	}
}
close FF;

# 输出
my @geneId = keys(%geneToAsGff3);
open GFF, ">$outputAsGff3";
foreach $geneId(@geneId){
#	print GFF join("\t", $geneSpanHref->{$geneId}->{"chr"}, "aslive", "gene", $geneSpanHref->{$geneId}->{"start"}, $geneSpanHref->{$geneId}->{"stop"}, ".", $geneSpanHref->{$geneId}->{"strand"}, ".", "ID=" . $geneId) . "\n";
	print GFF $geneToAsGff3Href->{$geneId};
}
close GFF;

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
