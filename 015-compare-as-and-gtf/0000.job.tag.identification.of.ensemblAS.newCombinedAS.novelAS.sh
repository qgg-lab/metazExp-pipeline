#!/bin/bash
#source $commonEnvir
#source $paras

export currDir=`pwd`
source $currDir/../../00-config/commonEnvironmentVars.cfg

export outputASTagFile=$currDir/AS.occur.in.ensembl.improvement.novel.tsv
export outputIsoformTagFile=$currDir/AS.two.exonSeries.to.isoformIdList.tsv


export ASTypeList=A3SS,A5SS,MXE,RI,SE
export ASCatalogList=$currDir/../012-generate-AS-catalog/A3SS.catalog.txt,$currDir/../012-generate-AS-catalog/A5SS.catalog.txt,$currDir/../012-generate-AS-catalog/MXE.catalog.txt,$currDir/../012-generate-AS-catalog/RI.catalog.txt,$currDir/../012-generate-AS-catalog/SE.catalog.txt
export ensemblGtf=$currDir/../001-prepare-local-datasource/nonRedundant.ensembl.gtf
export improvementGtf=$currDir/../006-form-final-trspt-annotation/final.complete.trspt.anno.gtf
export originOfAsTsv=$currDir/origin.of.as.tsv


perl $customizedSoftwareDir/assign.ensembl.improvement.novel.to.AS.pl \
		--ASTypeList $ASTypeList \
                --ASCatalogList $ASCatalogList \
                --ensemblGtf $ensemblGtf \
		--improvementGtf $improvementGtf \
		--outputTagFile $originOfAsTsv \
		1> $currDir/log.o.origin.of.as \
		2> $currDir/log.e.origin.of.as

export asWithLongAndShortAltTrsptListTsv=$currDir/as.with.long.short.alt.trsptIdList.tsv
perl $customizedSoftwareDir/assign.isoformIdList.to.AS.pl \
		--ASTypeList $ASTypeList \
		--ensemblGtf $ensemblGtf \
		--ASCatalogList $ASCatalogList \
		--improvementGtf $improvementGtf \
		--outputTagFile $asWithLongAndShortAltTrsptListTsv \
		1> log.o.assign.isoformIdList.to.AS \
		2> log.e.assign.isoformIdList.to.AS


export outputCoordinateTagFile=$currDir/AS.coordinate.in.isoform.coordinate.tsv
perl $customizedSoftwareDir/assign.symbol.whether.acceptor.and.donor.in.isoform.pl \
                --ASTypeList $ASTypeList \
                --ASCatalogList $ASCatalogList \
                --improvementGtf $improvementGtf \
                --outputTagFile $outputCoordinateTagFile

