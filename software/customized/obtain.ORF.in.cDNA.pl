#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--getorf getorf \\\n" .
                "--cDNAFasta new.assembled.trspt.cDNA.fa \\\n" .
                "--pepFasta  new.assembled.trspt.pep.fa \\\n" .
		"--tmpDir currDir/" .
		"--orfInCdnaTsv orf.in.new.assembled.trspt.cDNA.Tsv \n";
	exit;
}

my ($getorf, $cDNAFasta, $pepFasta, $orfInCdnaTsv, $tmpDir);

GetOptions(
        'getorf=s'=>\$getorf,
        'cDNAFasta=s'=>\$cDNAFasta,
	'tmpDir=s'=>\$tmpDir,
        'pepFasta=s'=>\$pepFasta,
        'orfInCdnaTsv=s'=>\$orfInCdnaTsv,
);

# 将转录本cDNA序列读入
my (%trsptCdnaSeq, @seqId, $seqId, $line);
open FF, "<$cDNAFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$seqId = $1;
		@seqId = split(/ /, $seqId);
		$seqId = $seqId[0];
	}else{
		$trsptCdnaSeq{$seqId}.=$line;
	}
}
close FF;

# getorf -sequence new.assembled.traspt.cDNA.fa -outseq new.assembled.traspt.pep.fa -reverse N -methionine Y -find 1
# -reverse N 不在反义链上找orf
# -methionine Y 强制起始密码子为ATG
# -find 1 报告的序列为为起始密码子和终止密码之间的序列
my $cmd;
$cmd = $getorf . 
	" -sequence " . $cDNAFasta . 
	" -outseq " . $tmpDir . "/tmp.rlt.of.getorf.in.new.assembled.trspt.fa " .
	" -reverse N " . 
	" -methionine Y " . 
	" -find 1";
system($cmd);

# 每个trspt上识别的orf有多个，只挑选最长orf作为最终的pep
#>SRX2909953.37098.2_1 [122 - 466]
#MEPSAWAVAFAALCIVALAPPASGFYLPGVAPNDFEKVRVRQPPFVSPPFLLLGSWHEIS
#SLALCVFGDRVRFGSAVGRPRSANARVGRIEIWLRVFATILALRAISGAQNVSVL
#>SRX2909953.37098.2_2 [466 - 516]
#MSRAKSFICASPVTFSS

my (%pep, $pepHref, $trsptId, $pepId, $startCodonBegin, $startCodonEnd, $stopCodonBegin, $stopCodonEnd, $codonPostOrf);
$pepHref=\%pep;
open FF, "<$tmpDir" . "/tmp.rlt.of.getorf.in.new.assembled.trspt.fa";
while(my $line=<FF>){
	chomp($line);
	if($line=~/>((.*?)_\d+) \[(\d+) \- (\d+)\]/){
		$trsptId = $2;
		$pepId = $1;
	
		# 判断ORF的第1个密码子是否为起始密码子
		if(uc(substr($trsptCdnaSeq{$trsptId}, $3-1, 3)) eq "ATG"){
			$startCodonBegin = $3;
			$startCodonEnd = $startCodonBegin + 2;
		}else{
			$startCodonBegin = -1;
			$startCodonEnd = -2;
		}

		# 判断ORF框之后的密码子是否为终止密码子，进一步获得终止密码子位置	
		if(length($trsptCdnaSeq{$trsptId}) < $4 + 3){
			# cDNAseq长度不足以支撑ORF之后存在终止密码子
			$stopCodonBegin = -1;
			$stopCodonEnd = -1;
		}else{
			# cDNAseq长度可以支撑ORF之后存在终止密码子，检查该密码子是否为终止密码子
			$codonPostOrf = substr($trsptCdnaSeq{$trsptId}, $4, 3);
			if(uc($codonPostOrf) eq "TAA" or uc($codonPostOrf) eq "TGA" or uc($codonPostOrf) eq "TAG"){
				$stopCodonBegin = $4 + 1;
				$stopCodonEnd = $stopCodonBegin + 2;		
			}else{
				$stopCodonBegin = -1;
				$stopCodonEnd = -1;
			}
		}

		$pepHref->{$trsptId}->{$pepId}->{"pepSeq"} = "";
		$pepHref->{$trsptId}->{$pepId}->{"orfBegin"} = $3;
		$pepHref->{$trsptId}->{$pepId}->{"orfEnd"} = $4;
		$pepHref->{$trsptId}->{$pepId}->{"startCodonBegin"} = $startCodonBegin;
		$pepHref->{$trsptId}->{$pepId}->{"startCodonEnd"} = $startCodonEnd;
		$pepHref->{$trsptId}->{$pepId}->{"stopCodonBegin"} = $stopCodonBegin;
		$pepHref->{$trsptId}->{$pepId}->{"stopCodonEnd"} = $stopCodonEnd;
	}else{
		$pepHref->{$trsptId}->{$pepId}->{"pepSeq"} .= $line;
	}
}
close FF;

# 在每个trspt内部找最长pep，然后将其放在maxPep的hash中
my (%maxPep, $maxPepHref, @trsptId, $trsptId, @pepId, $pepId);
$maxPepHref=\%maxPep;
@trsptId=keys(%$pepHref);
foreach $trsptId(@trsptId){
	@pepId = ();
	@pepId = keys(%{$pepHref->{$trsptId}});
#	print $trsptId . ":\n";
#	print join("\n",@pepId);
#	<STDIN>;
	foreach $pepId(@pepId){
		# 如果该转录本不存在最大pep，那么直接将第一次碰到的pep加入
#		print "trsptId:$trsptId". "\ttrsptCdnaSeq length:" . length($trsptCdnaSeq{$trsptId});
#		<STDIN>;
		if(not exists($maxPepHref->{$trsptId}->{"pepSeq"})){
			$maxPepHref->{$trsptId}->{"pepSeq"} = $pepHref->{$trsptId}->{$pepId}->{"pepSeq"};
			$maxPepHref->{$trsptId}->{"orfBegin"} = $pepHref->{$trsptId}->{$pepId}->{"orfBegin"};
			$maxPepHref->{$trsptId}->{"orfEnd"} = $pepHref->{$trsptId}->{$pepId}->{"orfEnd"};
			$maxPepHref->{$trsptId}->{"startCodonBegin"} = $pepHref->{$trsptId}->{$pepId}->{"startCodonBegin"};
			$maxPepHref->{$trsptId}->{"startCodonEnd"} = $pepHref->{$trsptId}->{$pepId}->{"startCodonEnd"};
			$maxPepHref->{$trsptId}->{"stopCodonBegin"} = $pepHref->{$trsptId}->{$pepId}->{"stopCodonBegin"};
			$maxPepHref->{$trsptId}->{"stopCodonEnd"} = $pepHref->{$trsptId}->{$pepId}->{"stopCodonEnd"};
			$maxPepHref->{$trsptId}->{"cDNAseqLen"} = length($trsptCdnaSeq{$trsptId});
			
		}elsif(length($pepHref->{$trsptId}->{$pepId}->{"pepSeq"}) > length($maxPepHref->{$trsptId}->{"pepSeq"})){
			$maxPepHref->{$trsptId}->{"pepSeq"} = $pepHref->{$trsptId}->{$pepId}->{"pepSeq"};
			$maxPepHref->{$trsptId}->{"orfBegin"} = $pepHref->{$trsptId}->{$pepId}->{"orfBegin"};
			$maxPepHref->{$trsptId}->{"orfEnd"} = $pepHref->{$trsptId}->{$pepId}->{"orfEnd"};
			$maxPepHref->{$trsptId}->{"startCodonBegin"} = $pepHref->{$trsptId}->{$pepId}->{"startCodonBegin"};
			$maxPepHref->{$trsptId}->{"startCodonEnd"} = $pepHref->{$trsptId}->{$pepId}->{"startCodonEnd"};
			$maxPepHref->{$trsptId}->{"stopCodonBegin"} = $pepHref->{$trsptId}->{$pepId}->{"stopCodonBegin"};
			$maxPepHref->{$trsptId}->{"stopCodonEnd"} = $pepHref->{$trsptId}->{$pepId}->{"stopCodonEnd"};
			$maxPepHref->{$trsptId}->{"cDNAseqLen"} = length($trsptCdnaSeq{$trsptId});
		}
	}
}

# 输出每个trspt中最长pep对应的orf
open WW, ">$pepFasta";
open POSITION, ">$orfInCdnaTsv";
print POSITION join("\t", "trsptId", "orfBegin", "orfEnd", "startCodonBegin", "startCodonEnd", "stopCodonBegin", "stopCodonEnd", "cDNAseqLen") . "\n";
@trsptId = keys(%maxPep);
foreach $trsptId (@trsptId){
	print WW ">$trsptId\n";
	print WW $maxPepHref->{$trsptId}->{"pepSeq"} . "\n";
	print POSITION join("\t", $trsptId, $maxPepHref->{$trsptId}->{"orfBegin"}, $maxPepHref->{$trsptId}->{"orfEnd"}, $maxPepHref->{$trsptId}->{"startCodonBegin"}, $maxPepHref->{$trsptId}->{"startCodonEnd"}, $maxPepHref->{$trsptId}->{"stopCodonBegin"}, $maxPepHref->{$trsptId}->{"stopCodonEnd"}, $maxPepHref->{$trsptId}->{"cDNAseqLen"}) . "\n";
}
close WW;
close POSITION;
