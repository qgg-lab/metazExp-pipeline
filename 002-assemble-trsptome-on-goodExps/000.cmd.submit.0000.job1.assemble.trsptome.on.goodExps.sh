# ~/software/hisat2-2.1.0/extract_splice_sites.py ../001-prepare-local-datasource/ensembl.gtf |awk -F '\t' '{print $3-$2}' |sort -nr |head -n 1
# ~/software/hisat2-2.1.0/extract_splice_sites.py ../001-prepare-local-datasource/ensembl.gtf |awk -F '\t' '{print $3-$2}' |sort -n |head -n 1
export currDir=`pwd`
export taxonId=`pwd | awk -F '/' '{print $(NF-1)}'`

ln -sf $currDir/../../00-config/commonEnvironmentVars.cfg $currDir/comEnvirVars


sbatch --job-name="002_"$taxonId \
	--export=commonEnvir=$currDir/comEnvirVars,paras=$currDir/00000.parameter.of.0000.job1.cfg \
	--nodes=1 \
	--ntasks-per-node=1 \
	--cpus-per-task=1 \
	--mem=64G \
	--time=00:05:00 \
	0000.job1.assemble.transcriptome.on.goodExps.sb
