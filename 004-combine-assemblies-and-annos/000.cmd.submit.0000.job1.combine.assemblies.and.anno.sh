export currDir=`pwd`
export taxonId=`pwd | awk -F '/' '{print $(NF-1)}'`

ln -sf $currDir/../../00-config/commonEnvironmentVars.cfg $currDir/comEnvirVars

sbatch --job-name=004_$taxonId \
	--export=commonEnvir=$currDir/comEnvirVars,paras=$currDir/00000.parameter.of.0000.job1.cfg \
	--nodes=1 \
	--ntasks-per-node=1 \
	--cpus-per-task=1 \
	--mem=60G \
	--time=12:00:00 \
	0000.job1.combine.assemblies.and.anno.sb
