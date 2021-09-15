#!/bin/bash

source ../../00-config/commonEnvironmentVars.cfg 

export taxonId=`pwd|awk -F '/' '{print $(NF-1)}'`
export getorf=$softwareDir/EMBOSS-6.6.0/bin/getorf
export currDir=`pwd`

export genomeFasta=$currDir/../001-prepare-local-datasource/ensembl.genome.fa
export finalGtf=$currDir/../006-form-final-trspt-annotation/final.complete.trspt.anno.gtf
export finalTrsptCdnaFasta=$currDir/../006-form-final-trspt-annotation/final.complete.trspt.cDNA.fa
export newAssembledTrsptCdnaFasta=$currDir/new.assembled.trspt.cDNA.fa
perl $customizedSoftwareDir/extract.new.assembled.trspt.from.final.cDNA.pl \
	--finalCdnaFasta $finalTrsptCdnaFasta \
	--finalGtf $finalGtf \
	--methodList StringTie \
	--newAssembledTrsptCdnaFasta  $newAssembledTrsptCdnaFasta

export newAssembledTrsptPepFasta=$currDir/new.assembled.trspt.pep.fa
export orfInNewAssembledTrsptCdnaTsv=$currDir/orf.in.new.assembled.trspt.cDNA.tsv
#:<<block2
perl $customizedSoftwareDir/obtain.ORF.in.cDNA.pl \
	--getorf $getorf \
	--cDNAFasta $newAssembledTrsptCdnaFasta \
	--pepFasta $newAssembledTrsptPepFasta \
	--tmpDir $currDir \
	--orfInCdnaTsv $orfInNewAssembledTrsptCdnaTsv
	#1> $currDir/log.o.obtain.ORF.startCodon.stopCodon.of.novelTrspt \
	#2> $currDir/log.e.obtain.ORF.startCodon.stopCodon.of.novelTrspt
rm -rf $currDir/log.*
#block2





export ensemblGtf=$currDir/../001-prepare-local-datasource/ensembl.gtf
export ensemblCdnaSeqFile=$currDir/../001-prepare-local-datasource/ensembl.cDNA.fa
export orfInEnsemblTrscriptCdnaTsv=$currDir/orf.in.ensembl.trspt.cDNA.tsv
#:<<block3
perl $customizedSoftwareDir/obtaion.ORF.in.cDNA.by.ensemblGtf.pl \
	--ensemblGtf $ensemblGtf \
	--cDNAseqFile $ensemblCdnaSeqFile \
	--orfInCdnaTsv $currDir/tmp.orf.in.ensembl.trspt.withRedundant.cDNA.tsv \
#	1> $currDir/log.o.obtaion.ORF.startCodon.stopCodon.of.ensemblTrspt \
#	2> $currDir/log.e.obtaion.ORF.startCodon.stopCodon.of.ensemblTrspt
grep -P "\tCDS\t" $currDir/../001-prepare-local-datasource/nonRedundant.ensembl.gtf \
	| awk -F '; ' '{print $2}' |awk -F '"' '{print $2}' \
	| sort -u > $currDir/nonRedunant.ensembl.trsptId.list
grep -Ff $currDir/nonRedunant.ensembl.trsptId.list $currDir/tmp.orf.in.ensembl.trspt.withRedundant.cDNA.tsv > $currDir/tmp.orf.in.ensembl.trspt.cDNA.tsv
#block3




#########  ==4 ===###############################
export pepFasta=$currDir/../006-form-final-trspt-annotation/final.complete.trspt.pep.fa
export cDNAfast=$currDir/../006-form-final-trspt-annotation/final.complete.trspt.cDNA.fa
module load BioPerl/1.7.2-Perl-5.26.1
#:<<block4
perl $customizedSoftwareDir/correct.orf.in.ensembl.trspt.pl \
	--origOrfInCdna $currDir/tmp.orf.in.ensembl.trspt.cDNA.tsv \
	--pepFasta $pepFasta \
	--cDNAFasta $cDNAfast \
	--correctOrfInCda $orfInEnsemblTrscriptCdnaTsv

module unload BioPerl/1.7.2-Perl-5.26.1
#block4






############  == 5 ==###############################
# 合并ensembl上的orf和stringTie新组装的orf
export orfInFinalTrsptCdnaTsv=$currDir/orf.in.final.ensembl.and.new.assembled.trspt.cDNA.tsv
#:<<block5
cat $orfInEnsemblTrscriptCdnaTsv $orfInNewAssembledTrsptCdnaTsv > $orfInFinalTrsptCdnaTsv
#block5



####  == 6 ==###############################
export ensemblTrsptPepFasta=$currDir/../006-form-final-trspt-annotation/final.complete.trspt.pep.fa
export finalEnsemblNewAssembledTrsptPepFasta=$currDir/final.ensembl.and.new.assembled.trspt.pep.fa
export interproscanAnnoDir=$currDir/interproscan
export pfamInFinalPepTsv=$currDir/pfam.in.final.ensembl.and.new.assembled.trspt.pep.tsv
#:<<block6
cat $newAssembledTrsptPepFasta $ensemblTrsptPepFasta > $finalEnsemblNewAssembledTrsptPepFasta
mkdir -p $interproscanAnnoDir
perl $customizedSoftwareDir/interproscan.annotate.pep.pl \
	--inputPepSeqFile $finalEnsemblNewAssembledTrsptPepFasta \
	--interproscan $interproscan \
	--outputDir $interproscanAnnoDir \
	--taxonId $taxonId
#	1> $currDir/log.o.annot.protein.by.interproscan \
#	2> $currDir/log.e.annot.protein.by.interproscan
#block6






#### ===== 7 将pfam在蛋白序列上的位置转换为在cDNA序列上的位置====################################
export pfamInFinalTrsptCdnaTsv=$currDir/pfam.in.final.ensembl.and.new.assembled.trspt.cDNA.tsv
:<<block7
cat $interproscanAnnoDir/*.tsv | grep -P "\tPfam\t" > $pfamInFinalPepTsv
perl $customizedSoftwareDir/convert.pfamPepPosition.to.pfamCdnaPositon.pl \
	--orfInCdnaTsv $orfInFinalTrsptCdnaTsv \
	--pfamInPepTsv $pfamInFinalPepTsv \
	--pfamInCdnaTsv $pfamInFinalTrsptCdnaTsv
#	1> $currDir/log.o.obtaion.pfamPosition.in.transcript \
#	2> $currDir/log.e.obtaion.pfamPosition.in.transcript
block7



export pfamSrtListBySiteInTrsptTsv=$currDir/pfam.sorted.list.by.site.in.trspt.tsv
:<<block8
perl $customizedSoftwareDir/obtain.pfam.sorted.list.by.site.in.trspt.pl \
	--pfamInCdnaTsv $pfamInFinalTrsptCdnaTsv \
	--pfamSrtListBySiteInTrsptTsv $pfamSrtListBySiteInTrsptTsv
#	1> $currDir/log.o.obtain.pfam.sorted.list.by.site.in.trspt
#	2> $currDir/log.e.obtain.pfam.sorted.list.by.site.in.trspt
block8



export goTermInTrsptTsv=$currDir/go.term.in.trspt.tsv
export goUniqTermInTrsptTsv=$currDir/go.uniq.term.trspt.tsv
:<<block9
cat $interproscanAnnoDir/*.tsv | grep GO: |awk -F '\t' '{print $1"\t"$(NF)}' > $goTermInTrsptTsv
perl $customizedSoftwareDir/extract.go.uniq.term.for.trspt.pl \
	--goTermInTrspt $goTermInTrsptTsv \
	--goUniqTermInTrspt $goUniqTermInTrsptTsv
#	1> $currDir/log.o.extract.goTerm.for.each.trspt \
#	2> $currDir/log.e.extract.goTerm.for.each.trspt
block9
