use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
		"--prepDE \\\n" .
                "--psiDir \\\n" .
		"--exptInfoTsv \\\n" .
                "--exptIdListFile \\\n" .
		"--outputDir \n";
        exit;
}

my ($psiDir, $exptIdListFile, $exptInfoTsv, $outputDir, $prepDE);

GetOptions(
	'prepDE=s'=>\$prepDE,
        'psiDir=s'=>\$psiDir,
	'exptInfoTsv=s'=>\$exptInfoTsv,
        'exptIdListFile=s'=>\$exptIdListFile,
	'outputDir=s'=>\$outputDir,
);


system("mkdir -p $outputDir/gene");
system("mkdir -p $outputDir/trspt");

# 读取实验信息获得实验测序的read长度
my (%exptIdToReadLen);
my (@exptId, $exptId, $abundanceFile, $assemblyFile);
my (@nameField, @valueField, $i, $line, %tmpHash);
open FF, "<$exptInfoTsv";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}
	$exptIdToReadLen{$tmpHash{"Experiment"}} = $tmpHash{"ReadLen"};
}
close FF;
# 创建gene目录，用于保存基因上的表达信息
# geneId  		readCount       Coverage        FPKM    	TPM
# EPlORYSAT000373851      1       	0.598199        0.888400        1.354598
# ENSRNA049446479 	0       	0.069444        0.103133        0.157253
# Os03g0586800    	1618    	36.534538       54.258366       82.731010
#
#transcriptId    readCount      Coverage        FPKM    	TPM
#Os12t0130300-00 81      	1.622807        0.565278        0.859528
#Os09t0122500-01 158     	12.023621       4.188230        6.368369

my (%geneExpr, $geneExprHref, %trsptExpr, $trsptExprHref);
my ($line, $geneId, $geneName, $ref, $strand, $start, $end, $cov, $fpkm, $tpm, $readLen, @field, $trsptId, @trsptId, $cmd, @geneId);
my ($geneSymbol, $readCount, $symbol);
open FF, "<$exptIdListFile";
@exptId = <FF>;
close FF;
foreach $exptId(@exptId){
	# 获得实验编号
	chomp($exptId);
	
	# 测序的read长度
	$readLen= $exptIdToReadLen{$exptId};

	# 清空表达数据hash
	%geneExpr = ();
	$geneExprHref = \%geneExpr;
	%trsptExpr = ();
	$trsptExprHref = \%trsptExpr;

	# 打开abundance文件，读取因的coverage, FPKM, TPM
	$abundanceFile = $psiDir . "$exptId/geneAbundanceByStringtie.tab";
	open FF, "<$abundanceFile";
	# Gene ID 	Gene Name       Reference       Strand  Start   End     Coverage        FPKM    	TPM
	# Os01g0100100  -       	1       	+       2983    10815   63.793186       6.685385        10.758520
	# Os01g0100200  -       	1       	+       11218   12435   4.926498        0.515535        0.829630
	<FF>;
	while($line=<FF>){
		chomp($line);
		($geneId, $geneName, $ref, $strand, $start, $end, $cov, $fpkm, $tpm) = split(/\t/, $line);
		$geneExprHref->{$geneId}->{"cov"} = $cov;
		$geneExprHref->{$geneId}->{"fpkm"} = $fpkm;
		$geneExprHref->{$geneId}->{"tpm"} = $tpm;
		$geneExprHref->{$geneId}->{"readCount"} = 0;
	}
	close FF;

	$assemblyFile = $psiDir . "$exptId/transcriptomeByStringtie.gtf";
	# 打开transcriptome文件，将transcript的coverage, FPKM, TPM读入
	open FF, "<$assemblyFile";
	while($line=<FF>){
		chomp($line);
		next if($line=~/^#/);
	
		@field = split(/\t/, $line);
		next if($field[2] ne "transcript");
		($geneId, $trsptId, $cov, $fpkm, $tpm) = ("", "", 0, 0, 0);
		&getAttr($field[8], \$geneId, \$trsptId, \$cov, \$fpkm, \$tpm);
		if($trsptId ne ""){
			$trsptExprHref->{$trsptId}->{"cov"} = $cov;
			$trsptExprHref->{$trsptId}->{"fpkm"} = $fpkm;
			$trsptExprHref->{$trsptId}->{"tpm"} = $tpm;
			$trsptExprHref->{$trsptId}->{"readCount"} = 0;
		}
	}
	close FF;

	# 创建指向transcriptom的文件列表
	open WW, ">file_position.txt.$exptIdListFile";
	print WW $exptId . " " . $assemblyFile;
	close WW;

	# 执行prepDE程序，获得gene和trspt上read数量
	$cmd = "python $prepDE -i file_position.txt.$exptIdListFile -g gene.count.tsv.$exptIdListFile -t trspt.count.tsv.$exptIdListFile -l $readLen";
	system($cmd);

	# 打开gene.count.tsv，将gene上的read数量累加进hash
	system("dos2unix gene.count.tsv.$exptIdListFile");
	open FF, "<gene.count.tsv.$exptIdListFile";
	<FF>;
	while($line = <FF>){
		chomp($line);
		# Os09g0127800|OsWD40-165,1810
		# Os03g0760800|GASR1,190
		($geneSymbol, $readCount) = split(/,/, $line);
		if($geneSymbol=~/(.*)\|(.*)/){
			$geneId = $1;
		}else{
			$geneId = $geneSymbol;
		}
		$geneExprHref->{$geneId}->{"readCount"} += $readCount;
	}
	close FF;

	# 打开trspt.count.tsv，将trspt上的read数读入hash
	system("dos2unix trspt.count.tsv.$exptIdListFile");
	open FF, "<trspt.count.tsv.$exptIdListFile";
	<FF>;
	while($line = <FF>){
		chomp($line);
		# transcript_id,SRX5808211
		# SRX5173465.28905.3,45
		($trsptId, $readCount) = split(/,/, $line);
		$trsptExprHref->{$trsptId}->{"readCount"} =$readCount;
	}
	close FF;

	# 将gene上表达数据写入到文件中
	@geneId = ();
	@geneId = keys(%geneExpr);
	open WW, ">$outputDir/gene/$exptId";
	print WW join("\t", "geneId", "readCount", "Coverage", "FPKM", "TPM") . "\n";
	foreach $geneId(@geneId){
		print WW join("\t", $geneId, $geneExprHref->{$geneId}->{"readCount"}, $geneExprHref->{$geneId}->{"cov"}, $geneExprHref->{$geneId}->{"fpkm"}, $geneExprHref->{$geneId}->{"tpm"}) . "\n";
	}
	close WW;

	# 将trspt上表达数据写入到文件中
	@trsptId = ();
	@trsptId = keys(%trsptExpr);
	open WW, ">$outputDir/trspt/$exptId";	
	print WW join("\t", "trsptId", "readCount", "Coverage", "FPKM", "TPM") . "\n";
	foreach $trsptId(@trsptId){
		print WW join("\t", $trsptId, $trsptExprHref->{$trsptId}->{"readCount"}, $trsptExprHref->{$trsptId}->{"cov"}, $trsptExprHref->{$trsptId}->{"fpkm"}, $trsptExprHref->{$trsptId}->{"tpm"}) . "\n";
	}
	close WW;
}

sub getAttr{
	my ($attrString, $geneId, $trsptId, $cov, $fpkm, $tpm) = @_;
	my (@attr, $attr);
	@attr = split(/;/, $attrString);
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$$geneId = $1;
		}elsif($attr=~/transcript_id "(.*)"/){
			$$trsptId = $1;
		}elsif($attr=~/cov "(.*)"/){
			$$cov = $1;
		}elsif($attr=~/FPKM "(.*)"/){
			$$fpkm = $1;
		}elsif($attr=~/TPM "(.*)"/){
			$$tpm = $1;
		}
	}
}
