#!/bin/bash

source ../../00-config/commonEnvironmentVars.cfg 

export taxonId=`pwd|awk -F '/' '{print $(NF-1)}'`
export getorf=$softwareDir/EMBOSS-6.6.0/bin/getorf
export currDir=`pwd`


export finalGtf=$currDir/../006-form-final-trspt-annotation/final.complete.trspt.anno.gtf
export asFileList=$currDir/../012-generate-AS-catalog/A3SS.catalog.txt,$currDir/../012-generate-AS-catalog/A5SS.catalog.txt,$currDir/../012-generate-AS-catalog/MXE.catalog.txt,$currDir/../012-generate-AS-catalog/RI.catalog.txt,$currDir/../012-generate-AS-catalog/SE.catalog.txt
export asTypeList=A3SS,A5SS,MXE,RI,SE
export asMappingTrsptTsv=$currDir/mapping.between.as.and.trspt.tsv

perl $customizedSoftwareDir/build.relation.between.AS.and.transcript.pl.20200725 \
	--gtf $finalGtf \
	--asFileList $asFileList \
	--asTypeList $asTypeList \
	--asMappingTrsptTsv $asMappingTrsptTsv

export processOfAsInCdnaTsv=$currDir/process.of.as.in.cDNA.tsv
export genomeFasta=$currDir/../001-prepare-local-datasource/ensembl.genome.fa

perl $customizedSoftwareDir/generate.as.edit.process.in.cDNA.pl.112509 \
	--asMappingTrsptTsv $asMappingTrsptTsv \
	--genomeFasta $genomeFasta \
	--asProcessInCdnaTsv $processOfAsInCdnaTsv

