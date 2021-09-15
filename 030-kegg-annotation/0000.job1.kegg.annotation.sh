#!/bin/bash
source ../../00-config/commonEnvironmentVars.cfg
export taxonId=`pwd|awk -F '/' '{print $(NF-1)}'`
export currDir=`pwd`

perl $customizedSoftwareDir/extract.longestPep.from.eachGene.pl.20200426 \
	--pepFasta $currDir/../006-form-final-trspt-annotation/final.complete.trspt.pep.fa \
	--gtf $currDir/../006-form-final-trspt-annotation/final.complete.trspt.anno.gtf \
	--longestPepFasta $currDir/$taxonId.longest.proteome.fa

