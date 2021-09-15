#!/bin/bash

awk -F '\t' '{print $9}' ../001-prepare-local-datasource/curated.experiment.tsv > orig.id.list

sed -i '/Experiment/d' orig.id.list

cat ../008-pickup-psi-of-ASs-in-all-expers/SRA.error.* orig.id.list |sort |uniq -u > sample.id.list

rm -rf orig.id.list


perl /mnt/home/liujind1/software/customized/extract.readCount.cov.tpm.for.gene.and.trspt.in.each.experiment.pl \
	--prepDE /mnt/home/liujind1/software/stringtie-2.1.4/prepDE.py \
	--psiDir ../008-pickup-psi-of-ASs-in-all-expers/psiOutputDir/ \
	--exptInfoTsv ../001-prepare-local-datasource/curated.experiment.tsv \
	--exptIdListFile sample.id.list\
	--outputDir ./expts



export totalJCECFileList=../012-generate-AS-catalog/jcec.A3SS.txt,../012-generate-AS-catalog/jcec.A5SS.txt,../012-generate-AS-catalog/jcec.MXE.txt,../012-generate-AS-catalog/jcec.RI.txt,../012-generate-AS-catalog/jcec.SE.txt
export totalJCFileList=../012-generate-AS-catalog/jc.A3SS.txt,../012-generate-AS-catalog/jc.A5SS.txt,../012-generate-AS-catalog/jc.MXE.txt,../012-generate-AS-catalog/jc.RI.txt,../012-generate-AS-catalog/jc.SE.txt

perl /mnt/home/liujind1/software/customized/extract.readCount.of.AS.for.each.experiment.pl \
	--totalJCECFileList $totalJCECFileList\
	--totalJCFileList $totalJCFileList\
	--exptIdListFile sample.id.list\
	--outputDir ./expts
