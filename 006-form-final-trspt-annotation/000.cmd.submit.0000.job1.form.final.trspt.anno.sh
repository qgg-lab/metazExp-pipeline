export currDir=`pwd`
export taxonId=`pwd | awk -F '/' '{print $(NF-1)}'`

ln -sf $currDir/../../00-config/commonEnvironmentVars.cfg $currDir/comEnvirVars

sbatch --job-name="006_"$taxonId \
	--export=commonEnvir=$currDir/comEnvirVars,paras=$currDir/00000.parameter.of.0000.job1.cfg,taxonId=$taxonId \
	--nodes=1 \
	--ntasks-per-node=1 \
	--cpus-per-task=1 \
	--mem=64G \
	--time=03:58:00 \
	0000.job1.form.final.trspt.anno.sb
