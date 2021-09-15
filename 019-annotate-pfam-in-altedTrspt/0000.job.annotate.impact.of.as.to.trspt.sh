#!/bin/bash

source ../../00-config/commonEnvironmentVars.cfg 

export taxonId=`pwd|awk -F '/' '{print $(NF-1)}'`
export getorf=$softwareDir/EMBOSS-6.6.0/bin/getorf
export currDir=`pwd`


export process=$currDir/../018-generate-process-of-AS-to-trspt/process.of.as.in.cDNA.tsv
export cDNAfasta=$currDir/../006-form-final-trspt-annotation/final.complete.trspt.cDNA.fa
export altedCdnaFasta=$currDir/alted.trspt.cDNA.fa
#:<<block1
perl $customizedSoftwareDir/generate.alted.trsptSeq.pl \
	--process $process \
	--cDNAfasta $cDNAfasta \
	--altedCdnaFasta $altedCdnaFasta
#	1>$currDir/log.o.generate.alted.trsptSeq \
#	2>$currDir/log.e.generate.alted.trsptSeq
#block1

export altedTrsptPepFasta=$currDir/alted.trspt.pep.fa
export orfInAltedTrsptCdnaTsv=$currDir/orf.in.alted.trspt.cDNA.tsv
#:<<block2
perl $customizedSoftwareDir/obtain.ORF.in.cDNA.pl \
        --getorf $getorf \
        --cDNAFasta $altedCdnaFasta \
        --pepFasta $altedTrsptPepFasta \
        --tmpDir $currDir \
        --orfInCdnaTsv $orfInAltedTrsptCdnaTsv
        #1> $currDir/log.o.obtain.ORF.startCodon.stopCodon.of.novelTrspt \
        #2> $currDir/log.e.obtain.ORF.startCodon.stopCodon.of.novelTrspt
#block2


export interproscanAnnoDir=$currDir/interproscan
:<<block3
mkdir -p $interproscanAnnoDir
perl $customizedSoftwareDir/interproscan.annotate.pep.pl \
        --inputPepSeqFile $altedTrsptPepFasta \
        --interproscan $interproscan \
        --outputDir $interproscanAnnoDir \
        --taxonId $taxonId \
	--jobName 019_ \
#       1> $currDir/log.o.annot.protein.by.interproscan \
#       2> $currDir/log.e.annot.protein.by.interproscan
block3



export pfamInFinalPepTsv=$currDir/pfam.in.alted.trspt.pep.tsv
:<<block4
cat $interproscanAnnoDir/*.tsv | grep -P "\tPfam\t" > $pfamInFinalPepTsv
block4


export orfInFinalTrsptCdnaTsv=$currDir/orf.in.alted.trspt.cDNA.tsv
export pfamInFinalTrsptCdnaTsv=$currDir/pfam.in.alted.trspt.cDNA.tsv
:<<block5
perl $customizedSoftwareDir/convert.pfamPepPosition.to.pfamCdnaPositon.pl \
        --orfInCdnaTsv $orfInAltedTrsptCdnaTsv \
        --pfamInPepTsv $pfamInFinalPepTsv \
        --pfamInCdnaTsv $pfamInFinalTrsptCdnaTsv
#       1> $currDir/log.o.obtaion.pfamPosition.in.transcript \
#       2> $currDir/log.e.obtaion.pfamPosition.in.transcript
block5



export goTermInTrspt=$currDir/go.term.in.alted.trspt.tsv
:<<block6
cat $interproscanAnnoDir/*.tsv | grep GO: |awk -F '\t' '{print $1"\t"$(NF)}' > $goTermInTrspt
block6


export goUniqTermInTrspt=$currDir/go.uniq.term.alted.trspt.tsv
:<<block7
perl $customizedSoftwareDir/extract.go.uniq.term.for.trspt.pl \
        --goTermInTrspt $goTermInTrspt \
        --goUniqTermInTrspt $goUniqTermInTrspt
#       1> $currDir/log.o.extract.goTerm.for.each.trspt \
#       2> $currDir/log.e.extract.goTerm.for.each.trspt
block7



export pfamSrtListBySiteInTrsptTsv=$currDir/pfam.sorted.list.by.site.in.trspt.tsv
:<<block8
perl $customizedSoftwareDir/obtain.pfam.sorted.list.by.site.in.trspt.pl \
        --pfamInCdnaTsv $pfamInFinalTrsptCdnaTsv \
        --pfamSrtListBySiteInTrsptTsv $pfamSrtListBySiteInTrsptTsv
#       1> $currDir/log.o.obtain.pfam.sorted.list.by.site.in.trspt
#       2> $currDir/log.e.obtain.pfam.sorted.list.by.site.in.trspt
block8
