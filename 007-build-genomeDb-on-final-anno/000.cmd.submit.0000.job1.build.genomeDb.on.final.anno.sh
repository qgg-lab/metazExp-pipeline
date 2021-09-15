export currDir=`pwd`
export taxonId=`pwd | awk -F '/' '{print $(NF-1)}'`

ln -sf $currDir/../../00-config/commonEnvironmentVars.cfg $currDir/comEnvirVars

sbatch --job-name="007_"$taxonId \
	--export=commonEnvir=$currDir/comEnvirVars,threadNum=18 \
	--nodes=1 \
	--ntasks-per-node=1 \
	--cpus-per-task=18 \
	--mem=64G \
	--time=03:50:00 \
	0000.job1.build.genomeDb.on.final.anno.sb
