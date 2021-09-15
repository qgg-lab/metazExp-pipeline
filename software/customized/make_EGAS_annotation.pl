#!/usr/bin/perl
use strict;
if($#ARGV<0){
	my $promotion=<<content;
	
	程序功能：
	
	该程序用于标准化GFF3格式文件。
	
	具体功能：
	
		（1）按照位置排序gene、mRNA、exon、CDS、5UTR、3UTR
		
		（2）计算CDS的phase（不管原来有没有提供phase）

		（3）统一编号gene、mRNA、exon、CDS、5UTR、3UTR，并且建立相互关联关系

			gene在全基因组范围内统一编号：ID=gene1;Parent=gene1

			mRNA在gene范围内统一编号：ID=mRNA1.gene1;Parent=gene1;Name=mRNA1.gene1

			exon在mRNA范围内统一编号：ID=exon1.mRNA1.gene1;Parent=mRNA1.gene1;Name=exon1.mRNA1.gene1

			CDS在mRNA范围内统一编号：ID=cds1.mRNA1.gene1;Parent=mRNA1.gene1;Name=cds1.mRNA1.gene1

			5utr在mRNA范围内统一编号：ID=5UTR1.mRNA1.gene1;Parent=mRNA1.gene1;Name=5UTR.mRNA1.gene1

			3utr在mRNA范围内统一编号：ID=3UTR1.mRNA1.gene1;Parent=mRNA1.gene1;Name=3UTR1.mRNA1.gene1

	输入gff3格式文件要求：

		（1）必须要有gene特征行，第9列属性中必须有ID=xxx

		（2）必须要求“mRNA”特征行，但不要求必须用“mRNA”标签（可以通过第4个输入参数指定）；第9列属性中必须有ID=xxx;Parent=yyy。mRNA可以指定用别的词汇（如transcript）替换。
		
		（3）计算方法标签（第2列）可以任意（不起任何标识作用）

		（4）exon部分可有可无，标签名可以为exon，也可以不用“exon”标签（可以通过第5个输入参数指定）。一旦有exon行，那么必须有ID=xxx;Parent=yyy。

		（5）CDS部分可有可无，标签名可以为CDS，也可不用“CDS”标签（可以通过第6个输入参数指定），必须有ID=xxx;Parent=yyy。

		（6）5utr部分可有可无，标签名可以为5utr，也可以不用“5utr”标签（可以通过第7个输入参数指定）。一旦有5utr行，那么必须有ID=xxx;Parent=yyy。

		（7）3utr部分可有可无，标签名可以为3utr，也可以不用“3utr”标签（可以通过第8个输入参数指定）。一旦有3utr行，那么必须有ID=xxx;Parent=yyy。

	输入参数格式：

		（1）输入的GFF3格式文件；
		
		（2）输出的GFF3格式文件名；
		
		（3）输入的GFF3格式文件中计算方法名称（对应第2列），此参数将不起任何作用；

		（4）输入的GFF3格式文件中的mRNA标签名（对应第3列）；
	
		（5）输入的GFF3格式文件中的exon标签名（对应第3列）；

		（6）输入的GFF3格式文件中的CDS标签名（对应第3列）；

		（7）输入的GFF3格式文件中的5UTR标签名（对应第3列）；

		（8）输入的GFF3格式文件中的3UTR标签名（对应第3列）；

		（9）输出的GFF3格式文件中计算方法名称（对应第2列），此参数将会作用；

		（10）输出的GFF3格式文件中的mRNA标签名（对应第3列）；
	
		（11）输出的GFF3格式文件中的exon标签名（对应第3列）；

		（12）输出的GFF3格式文件中的CDS标签名（对应第3列）；

		（13）输出的GFF3格式文件中的5UTR标签名（对应第3列）；

		（14）输出的GFF3格式文件中的3UTR标签名（对应第3列）；


	处理流程：

	1、提取gene行，按照位置进行排序，产生一个由有序gene行构成的文件。因此：输入文件必须有gene特征

	2、在有序gene行文件中提取基因编号，并赋予新的统一编号。因此：输入gene中必须有ID=特征


	3、提取mRNA行，并按照位置进行排序，产生一个由有序mRNA行构成的文件。因此：必须有mRNA特征。

	4、在排序后mRNA行文件中提取mRNA编号，并在所属gene下对当前mRNA进行统一编号。因此：mRNA必须有ID=xx;Parent=yy;


	5、提取exon行，并且按照位置进行排序（正链从小到大，负链从大到小）。因此：exon特征可有可无；

	6、在排序后的exon行文件中提取exon编号，并且把exon行关联到所属的mRNA中。因此，如果有exon特征，那么属性必须有ID=xx;Parent=yy


	7、提取CDS行，并且按照位置进行排序（正链从小到大，负链从大到小）。因此：CDS特征可有可无；

	8、在排序后的CDS行文件中提取CDS编号，并且把CDS行关联到所属的mRNA中。因此，如果有CDS特征，那么属性必须有ID=xx;Parent=yy


	9、提取3utr行，并且按照位置进行排序（正链从小到大，负链从大到小）。因此：3utr特征可有可无；

	10、在排序后的3utr行文件中提取3utr编号，并且把3utr行关联到所属的mRNA中。因此，如果有3utr特征，那么属性必须有ID=xx;Parent=yy


	11、提取5utr行，并且按照位置进行排序（正链从小到大，负链从大到小）。因此：5utr特征可有可无；

	12、在排序后的5utr行文件中提取5utr编号，并且把5utr行关联到所属的mRNA中。因此，如果有5utr特征，那么属性必须有ID=xx;Parent=yy


	13、提取所有mRNA的CDS，计算每个CDS的phase。

	14、对所有mRNA进行：
	
		（1）提取当前mRNA的exon文本，然后对其exon进行内部统一编号。
			编号格式：ID=exon[mRNA内部编号].mRNA[gene内部编号].gene[整个基因组内部gene编号]
		（2）提取当前mRNA的CDS文本，然后对其CDS进行内部统一编号。
			编号格式：ID=cds[mRNA内部编号].mRNA[gene内部编号].gene[整个基因组内部gene编号]
		（3）提取当前mRNA的5UTR文本，然后对其5UTR进行内部统一编号。
			编号格式：ID=5UTR[mRNA内部编号].mRNA[gene内部编号].gene[整个基因组内部gene编号]
		（4）提取当前mRNA的3UTR文本，然后对其3UTR进行内部统一编号。
			编号格式：ID=3UTR[mRNA内部编号].mRNA[gene内部编号].gene[整个基因组内部gene编号]
		（5）合并exon、CDS、5utr、3utr文本到一起，然后按照位置进行排序。排序后再次关联到mRNA中。

	15、最后输出所有的gene信息。

content
;
print $promotion . "\n\n";
print "./make_EGAS_annotation.pl 输入gff3文件名 输出gff3文件名 不起作用的原方法名称 原mRNA对应的标签名 原exon对应的标签 原CDS对应的标签 原5utr对应的标签 原3utr对应的标签 新方法名称（有效） 新mRNA标签 新exon标签 新CDS标签 新5UTR标签 新3UTR标签\n\n\n";
print "如果输入GFF3格式文件中缺少CDS、exon、5UTR、3UTR等元件时，那么用任意不冲突的字符串占用相关位置\n\n\n";
exit;

}

my $input_gff3 = $ARGV[0];
my $output_gff3= $ARGV[1];

my $old_label_source = $ARGV[2];
my $old_label_mRNA = $ARGV[3];
my $old_label_exon = $ARGV[4];
my $old_label_cds = $ARGV[5];
my $old_label_5utr = $ARGV[6];
my $old_label_3utr = $ARGV[7];

my $label_source = $ARGV[8];
my $label_mRNA = $ARGV[9];
my $label_exon = $ARGV[10];
my $label_cds = $ARGV[11];
my $label_5utr = $ARGV[12];
my $label_3utr = $ARGV[13];
my $j;
#扫描输入文件，读取所有类型为gene的特征
my $cmd = "grep -P \"\\tgene\\t\" $input_gff3 > gene.gff3";
system($cmd);
#首先按照参考序列排列特征，相同参考序列时再按开始start对特征进行排列
#排序的目的是输出时按照位置先后输出,且同一个scaffold_id上的gene聚集在一起
$cmd = "sort -t\'	\' -k1,1 -k4,4n gene.gff3 > sorted_gene.gff3";
system($cmd);
#依次读取sorted_gene.gff3，然后给每个基因进行编号
my ($line, %gene, @sorted_geneid, $i, $chain, $start, $end);
my ($gene_id, $gene_name, $scaffold_id, $attribute);
my (@attribute, $gene_num, $new_gene_id);
open FSGF, "<sorted_gene.gff3";
while(<FSGF>){
	$line = $_;
	if($line=~/^(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/){
		$scaffold_id = $1;
		$start = $4;
		$end = $5;
		$chain = $7;
		$attribute = $9;
		#用";"分解attribute属性部分，然后提取其中的基因编号和基因名称
		#这样按照次序保存基因编号，后输出的基因就是按照位置排序了
		@attribute = ();
		@attribute = split(/;/, $attribute);
		for($i=0; $i<=$#attribute; $i++){
			if($attribute[$i]=~/ID=(.*)/){
				$gene_id = $1;
				$gene_num++;
##				#新基因的名称设置位置
				$new_gene_id = "gene" . $gene_num;
			}
			if($attribute[$i]=~/Name=(.*)/){
				$gene_name=$1;
			}
		}
		${$gene{$gene_id}}{"new_gene_id"}= $new_gene_id;
		#基因名采用新的基因编号
		${$gene{$gene_id}}{"gene_name"}= $new_gene_id;
		${$gene{$gene_id}}{"scaffold_id"}= $scaffold_id;
		${$gene{$gene_id}}{"start"}= $start;
		${$gene{$gene_id}}{"end"}= $end;
		${$gene{$gene_id}}{"chain"}= $chain;
		#按照原来次序登记基因的编号，以便后续的输出
		$sorted_geneid[$#sorted_geneid+1]=$gene_id;
	}
}
close FSGF;
#
#print "gene\n";
#<STDIN>;
#输出所有gene feature
#for($i=0; $i<=$#sorted_geneid; $i++){
#	print ${$gene{$sorted_geneid[$i]}}{"scaffold_id"} . "\t" . $label_source . "\t" . "gene" . "\t" . ${$gene{$sorted_geneid[$i]}}{"start"} . "\t" . ${$gene{$sorted_geneid[$i]}}{"end"} . "\t" . "." . "\t" . ${$gene{$sorted_geneid[$i]}}{"chain"} . "\t" . "." . "\t" . "ID=" . ${$gene{$sorted_geneid[$i]}}{"new_gene_id"} . ";Name=" . ${$gene{$sorted_geneid[$i]}}{"gene_name"} . "\n";
#}
#<STDIN>;
#
#

#扫描输入文件，读取所有类型为mRNA的特征
$cmd = "grep -P \"\\t$old_label_mRNA\\t\" $input_gff3 > mRNA.gff3";
system($cmd);
#排序的目的是同一个gene的mRNA可以按照位置先后进行编号
$cmd = "sort -t\'	\' -k1,1 -k4,4n mRNA.gff3 > sorted_mRNA.gff3";
system($cmd);

my (%mRNA, $mRNA_id, $mRNA_name, @sorted_mRNAid);
open FF, "<sorted_mRNA.gff3";
while(<FF>){
	$line=$_;
	if($line=~/^(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/){
		$scaffold_id = $1;
		$start = $4;
		$end = $5;
		$chain = $7;
		$attribute = $9;
		#分解attribute部分，提取ID、Name和Parent
		#这样可以利用Parent把当前mRNA关联到gene中去
		@attribute = ();
		@attribute = split(/;/, $attribute);
		for($i=0; $i<=$#attribute; $i++){
			if($attribute[$i]=~/ID=(.*)/){
				$mRNA_id = $1;
			}
			if($attribute[$i]=~/Name=(.*)/){
				$mRNA_name=$1;
			}
			if($attribute[$i]=~/Parent=(.*)/){
				$gene_id=$1;
			}
		}
		#按照位置登记mRNAid，以便后续对mRNA的处理
		$sorted_mRNAid[$#sorted_mRNAid+1]=$mRNA_id;
		
		#登记此mRNA所在的基因编号
		${$mRNA{$mRNA_id}}{"gene_id"}= ${$gene{$gene_id}}{"new_gene_id"};
		#对gene内部mRNA进行统一编号
		${$gene{$gene_id}}{"mRNA_num"}=	${$gene{$gene_id}}{"mRNA_num"} + 1;
		${$mRNA{$mRNA_id}}{"new_mRNA_id"}= "mRNA" . ${$gene{$gene_id}}{"mRNA_num"} . "." . ${$gene{$gene_id}}{"new_gene_id"};
		${$mRNA{$mRNA_id}}{"mRNA_name"}=${$mRNA{$mRNA_id}}{"new_mRNA_id"};

		#把当前mRNA的原始编号添加到对应基因的mRNAid_list中，以便于从基因能够找到它所包含的mRNA
		${$gene{$gene_id}}{"mRNAid_list"}=${$gene{$gene_id}}{"mRNAid_list"} . $mRNA_id . "#";
		
		#登记其他信息
		${$mRNA{$mRNA_id}}{"scaffold_id"}= $scaffold_id;
		${$mRNA{$mRNA_id}}{"start"}= $start;
		${$mRNA{$mRNA_id}}{"end"}= $end;
		${$mRNA{$mRNA_id}}{"chain"}= $chain;
		
	}
}
close FF;
#
#print "mRNA\n";
#<STDIN>;
#for($i=0; $i<=$#sorted_mRNAid; $i++){
#	 print ${$mRNA{$sorted_mRNAid[$i]}}{"scaffold_id"} . "\t" . $label_source . "\t" . "$label_mRNA" . "\t" . ${$mRNA{$sorted_mRNAid[$i]}}{"start"} . "\t" . ${$mRNA{$sorted_mRNAid[$i]}}{"end"} . "\t" . "." . "\t" . ${$mRNA{$sorted_mRNAid[$i]}}{"chain"} . "\t" . "." . "\t" . "ID=" . ${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"} . ";Parent=" . ${$mRNA{$sorted_mRNAid[$i]}}{"gene_id"} . ";Name=" . ${$mRNA{$sorted_mRNAid[$i]}}{"mRNA_name"} . "\n";
#}
#<STDIN>;
#

#扫描输入文件，读取所有类型列为exon的特征
$cmd = "grep -P \"\\t$old_label_exon\\t\" $input_gff3 > exon.gff3";
system($cmd);
#按照坐标排序外显子，这样可以使得同一个mRNA的外显子能够有序：正链从小到大、负链从大到小（特殊处理）
$cmd = "sort -t\'	\' -k1,1 -k4,4n exon.gff3 > sorted_exon.gff3";
system($cmd);

my ($exon_id, $exon_name, $exon_line);
open FF, "<sorted_exon.gff3";
while(<FF>){
	$line=$_;
	if($line=~/^(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/){
		$scaffold_id = $1;
		$start = $4;
		$end = $5;
		$chain = $7;
		$attribute = $9;
		#用";"分解attribute属性，提取ID、Name和Parent
		#根据Parent可以把外显子关联到mRNA中去
		@attribute = ();
		@attribute = split(/;/, $attribute);
		for($i=0; $i<=$#attribute; $i++){
			if($attribute[$i]=~/ID=(.*)/){
				$exon_id = $1;
			}
			if($attribute[$i]=~/Name=(.*)/){
				$exon_name=$1;
			}
			if($attribute[$i]=~/Parent=(.*)/){
				$mRNA_id=$1;
			}
		}
		#外显子特征串,没有包含最后属性列，主要原因是外显子编号无法确定。确定外显子编号工作单独处理。
		$exon_line = $scaffold_id . "\t" . $label_source . "\t" . $label_exon . "\t" . $start . "\t" . $end . "\t" . "." . "\t" . $chain . "\t" . "." . "\n";
		#把外显子添加到相应mRNA的外显子文本exon_text中,exon_text是多行文本。添加工作就是关联工作
		#正链正向排列，负链反向排列
		if($chain eq "+"){
			${$mRNA{$mRNA_id}}{"exon_text"}=${$mRNA{$mRNA_id}}{"exon_text"} . $exon_line;
		}else{
			${$mRNA{$mRNA_id}}{"exon_text"}=$exon_line .  ${$mRNA{$mRNA_id}}{"exon_text"};
		}
	}
}
close FF;

#
#输出每个mRNA的外显子
#
#print "输出mRNA的外显子\n";
#<STDIN>;
#for($i=0; $i<=$#sorted_mRNAid; $i++){
#	print ${$mRNA{$sorted_mRNAid[$i]}}{"exon_text"};
#}
#<STDIN>;

#扫描输入文件，读取所有类型列为CDS的特征
$cmd = "grep -P \"\\t$old_label_cds\\t\" $input_gff3 > cds.gff3";
system($cmd);
#按照坐标排序CDS，这样可以使得同一个mRNA的CDS能够有序：正链从小到大、负链从大到小（另外处理）
$cmd = "sort -t\'	\' -k1,1 -k4,4n cds.gff3 > sorted_cds.gff3";
system($cmd);
my ($cds_id, $cds_name, $cds_line);
open FF, "<sorted_cds.gff3";
while(<FF>){
	$line=$_;
	if($line=~/^(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/){
		$scaffold_id = $1;
		$start = $4;
		$end = $5;
		$chain = $7;
		#用";"分解attribute，获得ID、Name和Parent
		#这样可以利用Parent把CDS特征关联到相应的mRNA中
		$attribute = $9;
		@attribute = ();
		@attribute = split(/;/, $attribute);
		for($i=0; $i<=$#attribute; $i++){
			if($attribute[$i]=~/ID=(.*)/){
				$cds_id = $1;
			}
			if($attribute[$i]=~/Name=(.*)/){
				$cds_name=$1;
			}
			if($attribute[$i]=~/Parent=(.*)/){
				$mRNA_id=$1;
			}
		}
		#cds特征串,没有包含最后属性列，主要原因是cds编号无法确定。确定cds编号工作单独处理。
		$cds_line = $scaffold_id . "\t" . $label_source . "\t" . $label_cds . "\t" . $start . "\t" . $end . "\t" . "." . "\t" . $chain . "\t" . "." . "\n";
		#把CDS添加到相应mRNA的CDS文本cds_text中,cds_text是多行文本。添加工作就是关联工作
		#正链正向排列，负链反向排列
		if($chain eq "+"){
			${$mRNA{$mRNA_id}}{"cds_text"}=${$mRNA{$mRNA_id}}{"cds_text"} . $cds_line;
		}else{
			${$mRNA{$mRNA_id}}{"cds_text"}=$cds_line . ${$mRNA{$mRNA_id}}{"cds_text"};
		}
	}
}
close FF;


#扫描输入文件，读取所有类型列为3utr的特征

###***************根据gff3中的3UTR表示符号做相应改变**************************###
$cmd = "grep -P \"\\t$old_label_3utr\\t\" $input_gff3 > utr3.gff3";
system($cmd);
#按照坐标排序3UTR，这样可以使得同一个mRNA的3UTR能够有序：正链从小到大、负链从大到小（特殊处理）
$cmd = "sort -t\'	\' -k1,1 -k4,4n utr3.gff3 > sorted_utr3.gff3";
system($cmd);
my ($utr3_id, $utr3_name, $utr3_line);
open FF, "<sorted_utr3.gff3";
while(<FF>){
	$line=$_;
	if($line=~/^(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/){
		$scaffold_id = $1;
		$start = $4;
		$end = $5;
		$chain = $7;
		$attribute = $9;
		#用";"分解attribute属性，获得ID、Name和Parent
		#利用Parent把3UTR和mRNA关联起来
		@attribute = ();
		@attribute = split(/;/, $attribute);
		for($i=0; $i<=$#attribute; $i++){
			if($attribute[$i]=~/ID=(.*)/){
				$utr3_id = $1;
			}
			if($attribute[$i]=~/Name=(.*)/){
				$utr3_name=$1;
			}
			if($attribute[$i]=~/Parent=(.*)/){
				$mRNA_id=$1;
			}
		}
		#utr3特征串,没有包含最后属性列，主要原因是utr3编号无法确定。确定utr3编号工作单独处理。
		$utr3_line = $scaffold_id . "\t" . $label_source . "\t" . $label_3utr . "\t" . $start . "\t" . $end . "\t" . "." . "\t" . $chain . "\t" . "." . "\n";
		#把3UTR添加到相应mRNA的3UTR文本utr3_text中,utr3_text是多行文本。添加工作就是关联工作
                #正链正向排列，负链反向排列
		if($chain eq "+"){
			${$mRNA{$mRNA_id}}{"utr3_text"}=${$mRNA{$mRNA_id}}{"utr3_text"} . $utr3_line;
		}else{
			${$mRNA{$mRNA_id}}{"utr3_text"}=$utr3_line . ${$mRNA{$mRNA_id}}{"utr3_text"};
		}
	}
}
close FF;

#print "输出mRNA的3UTR\n";
#<STDIN>;
#for($i=0; $i<=$#sorted_mRNAid; $i++){
#        print ${$mRNA{$sorted_mRNAid[$i]}}{"utr3_text"};
#}
#<STDIN>;

#扫描输入文件，读取所有类型列为5utr的特征
$cmd = "grep -P \"\\t$old_label_5utr\\t\" $input_gff3 > utr5.gff3";
system($cmd);
$cmd = "sort -t\'	\' -k1,1 -k4,4n utr5.gff3 > sorted_utr5.gff3";
system($cmd);
my ($utr5_id, $utr5_name, $utr5_line);
open FF, "<sorted_utr5.gff3";
while(<FF>){
	$line=$_;
	if($line=~/^(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/){
		$scaffold_id = $1;
		$start = $4;
		$end = $5;
		$chain = $7;
		$attribute = $9;
		@attribute = ();
		@attribute = split(/;/, $attribute);
		for($i=0; $i<=$#attribute; $i++){
			if($attribute[$i]=~/ID=(.*)/){
				$utr5_id = $1;
			}
			if($attribute[$i]=~/Name=(.*)/){
				$utr5_name=$1;
			}
			if($attribute[$i]=~/Parent=(.*)/){
				$mRNA_id=$1;
			}
		}
		#utr5特征串,没有包含最后属性列，主要原因是utr5编号无法确定。确定utr5编号工作单独处理。
		$utr5_line = $scaffold_id . "\t" . $label_source . "\t" . $label_5utr . "\t" . $start . "\t" . $end . "\t" . "." . "\t" . $chain . "\t" . "." . "\n";
		#把5UTR添加到相应mRNA的5UTR文本utr5_text中,utr5_text是多行文本。添加工作就是关联工作
		#正链正向排列，负链反向排列
		if($chain eq "+"){
			${$mRNA{$mRNA_id}}{"utr5_text"}=${$mRNA{$mRNA_id}}{"utr5_text"} . $utr5_line;
		}else{
			${$mRNA{$mRNA_id}}{"utr5_text"}=$utr5_line . ${$mRNA{$mRNA_id}}{"utr5_text"};
		}
	}
}
close FF;


#print "开始处理CDS的phase\n";
my ($cds_text, $phase);
#处理CDS上的phase(经过之前处理，正链CDS按照位置从小到大排列，负链CDS按照位置从大到小排列了。而且CDS自身编号也已经完成)
my ($cds_curr_len,$cds_previous_len, @column, @feature, $attribute);
for($i=0; $i<=$#sorted_mRNAid; $i++){
	@feature=();
        @feature=split(/\n/,${$mRNA{$sorted_mRNAid[$i]}}{"cds_text"});
	#列举所有CDS的特征
	$cds_text="";
	$cds_curr_len=0;
	$cds_previous_len=0;
	for($j=0; $j<=$#feature; $j++){
		#分解特征上的所有列，提取每个CDS特征长度，计算每个CDS特征的phase
		@column=();
		@column=split(/\t/, $feature[$j]);
		$scaffold_id = $column[0];
		$start = $column[3];
		$end = $column[4];
		$chain = $column[6];
		$cds_curr_len = $cds_previous_len + $end - $start + 1;
		#计算phase
		$phase = 3 - ($cds_previous_len - int($cds_previous_len / 3)*3);
		if($phase == 3){
			$phase = 0;
		}
		$cds_text = $cds_text . $scaffold_id . "\t" . $label_source . "\t" . $label_cds . "\t" . $start . "\t" . $end . "\t" . "." . "\t" . $chain . "\t" . $phase . "\n";
		$cds_previous_len = $cds_curr_len;
	}
	${$mRNA{$sorted_mRNAid[$i]}}{"cds_text"} = $cds_text;
}

#对exon, cds, 5utr, 3utr进行mRNA内部编号,并且添加上ID、Parent和Name
#最后把exon、cds、5utr、3utr合并一起，排序
my ($i, $exon_text, $cds_text, $utr5_text, $utr3_text, @temp, $j, $id, $parent, $name, $number, $total_number);
for($i=0; $i<=$#sorted_mRNAid; $i++){
	#处理exon
	@temp=();
	@temp=split(/\n/,${$mRNA{$sorted_mRNAid[$i]}}{"exon_text"});
	$number = 1;
	$total_number = $#temp +1;
	$exon_text = "";
	for($j=0; $j<=$#temp; $j++){
		$id= "exon" . $number . "." . ${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$number ++;
		$parent=${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$name = $id;
		#添加到字符串$exon_text
		$exon_text = $exon_text . $temp[$j] . "\t" . "ID=" . $id . ";Parent=" . $parent . ";Name=" . $name . "\n";
	}
	${$mRNA{$sorted_mRNAid[$i]}}{"exon_text"} = $exon_text;
	#print "exon\n";
	#<STDIN>;
	#print $exon_text;
	#<STDIN>;

	#处理CDS
	$cds_text="";
	@temp=();
	@temp=split(/\n/,${$mRNA{$sorted_mRNAid[$i]}}{"cds_text"});
	$number = 1;
	$total_number = $#temp +1;
	for($j=0; $j<=$#temp; $j++){
		$id= "cds" . $number . "." . ${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$number++;
		$parent=${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$name = $id;
		#添加到字符串$cds_text
		$cds_text = $cds_text . $temp[$j] . "\t" . "ID=" . $id . ";Parent=" . $parent . ";Name=" . $name . "\n";
	}
	${$mRNA{$sorted_mRNAid[$i]}}{"cds_text"} = $cds_text;
	
	#print "cds\n";
	#<STDIN>;
	#print $cds_text;
	#<STDIN>;

	#处理5utr
	$utr5_text="";
	@temp=();
	@temp=split(/\n/,${$mRNA{$sorted_mRNAid[$i]}}{"utr5_text"});
	$number = 1;
	$total_number = $#temp +1;
	for($j=0; $j<=$#temp; $j++){
		$id= "5UTR" . $number . "." . ${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$number++;
		$parent=${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$name = $id;
		#添加到字符串$utr5_text
		$utr5_text = $utr5_text . $temp[$j] . "\t" . "ID=" . $id . ";Parent=" . $parent . ";Name=" . $name . "\n";
	}
	${$mRNA{$sorted_mRNAid[$i]}}{"utr5_text"} = $utr5_text;
	#print "utr5\n";
	#<STDIN>;
	#print $utr5_text;
	#<STDIN>;

	#处理3utr
	$utr3_text="";
	@temp=();
	@temp=split(/\n/,${$mRNA{$sorted_mRNAid[$i]}}{"utr3_text"});
	$number = 1;
	$total_number = $#temp +1;
	for($j=0; $j<=$#temp; $j++){
		#$id= $number . "." . ${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$id= "3UTR" . $number . "." . ${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$number++;	
		$parent=${$mRNA{$sorted_mRNAid[$i]}}{"new_mRNA_id"};
		$name = $id;
		#添加到字符串$utr3_text
		$utr3_text = $utr3_text . $temp[$j] . "\t" . "ID=" . $id . ";Parent=" . $parent . ";Name=" . $name . "\n";
	}
	${$mRNA{$sorted_mRNAid[$i]}}{"utr3_text"} = $utr3_text;
	
	#print "utr3\n";
	#<STDIN>;
	#print $utr3_text;

	#把mRNA的exon、5utr、3utr以及cds特征合并到一起，重新按照坐标排序
	open WEUC, ">exon_utr5_3_cds.gff3";
	print WEUC $exon_text;
	print WEUC $cds_text;
	print WEUC $utr3_text;
	print WEUC $utr5_text;
	close WEUC;
	
	#按照坐标排序,正链正序，负链反序
	if(${$mRNA{$sorted_mRNAid[$i]}}{"chain"} eq "+"){
		system("sort -t\'	\' -k4,4n exon_utr5_3_cds.gff3 > sorted_exon_utr5_3_cds.gff3");
	}else{
		system("sort -t\'	\' -k4,4nr exon_utr5_3_cds.gff3 > sorted_exon_utr5_3_cds.gff3");
	}

	#把排序后的内容读入到mRNA中
	#print "排序完成\n";
	#<STDIN>;
	open FF, "<sorted_exon_utr5_3_cds.gff3";
	while(<FF>){
		${$mRNA{$sorted_mRNAid[$i]}}{"exon_cds_utr5_utr3"} = ${$mRNA{$sorted_mRNAid[$i]}}{"exon_cds_utr5_utr3"} . $_;
	}
	close FF;

	#print "exon_cds_utr5_utr3\n";
	#<STDIN>;
	#print ${$mRNA{$sorted_mRNAid[$i]}}{"exon_cds_utr5_utr3"};
#	<STDIN>;
}

#依次读取基因编号，然后根据基因编号输出基因的特征；
#根据基因编号从基因中提取mRNA编号，根据mRNA编号，输出mRNA的特征和mRNA下的所有exon、cds、5utr、3utr特征
#<STDIN>;
#print "输出最终结构\n";
#
open WSD, ">$output_gff3";
my ($mRNAid_list);
for($i=0; $i<=$#sorted_geneid; $i++){
	#输出基因结构
	print WSD ${$gene{$sorted_geneid[$i]}}{"scaffold_id"} . "\t" . $label_source . "\t" . "gene" . "\t" . ${$gene{$sorted_geneid[$i]}}{"start"} . "\t" . ${$gene{$sorted_geneid[$i]}}{"end"} . "\t" . "." . "\t" . ${$gene{$sorted_geneid[$i]}}{"chain"} . "\t" . "." . "\t" . "ID=" . ${$gene{$sorted_geneid[$i]}}{"new_gene_id"} . ";Name=" . ${$gene{$sorted_geneid[$i]}}{"gene_name"} . "\n";

	#提取mRNAid_list
	$mRNAid_list =  ${$gene{$sorted_geneid[$i]}}{"mRNAid_list"};
	@temp=();
	@temp=split(/#/,$mRNAid_list);

	#依次输出当前gene下的所有mRNA的特征
	for($j=0; $j<=$#temp; $j++){
		$mRNA_id = $temp[$j];
		#输出mRNA的结构
		print WSD ${$mRNA{$mRNA_id}}{"scaffold_id"} . "\t" . $label_source . "\t" . $label_mRNA . "\t" . ${$mRNA{$mRNA_id}}{"start"} . "\t" . ${$mRNA{$mRNA_id}}{"end"} . "\t" . "." . "\t" . ${$mRNA{$mRNA_id}}{"chain"} . "\t" . "." . "\t" . "ID=" . ${$mRNA{$mRNA_id}}{"new_mRNA_id"} . ";Parent=" . ${$mRNA{$mRNA_id}}{"gene_id"} . ";Name=" . ${$mRNA{$mRNA_id}}{"mRNA_name"} . "\n";
		#输出mRNA下的exon/CDS/5UTR/3UTR
		print WSD ${$mRNA{$mRNA_id}}{"exon_cds_utr5_utr3"};
	}
}
close WSD;

system("rm -rf cds.gff3 exon.gff3 exon_utr5_3_cds.gff3 gene.gff3 mRNA.gff3 sorted_cds.gff3 sorted_exon.gff3 sorted_exon_utr5_3_cds.gff3 sorted_gene.gff3 sorted_mRNA.gff3 sorted_utr3.gff3 sorted_utr5.gff3 utr3.gff3 utr5.gff3");
