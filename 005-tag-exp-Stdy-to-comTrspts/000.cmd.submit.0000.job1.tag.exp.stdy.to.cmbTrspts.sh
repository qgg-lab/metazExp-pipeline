export currDir=`pwd`
export taxonId=`pwd | awk -F '/' '{print $(NF-1)}'`

ln -sf $currDir/../../00-config/commonEnvironmentVars.cfg $currDir/comEnvirVars

sbatch --job-name=005_$taxonId \
	--export=commonEnvir=$currDir/comEnvirVars \
	--nodes=1 \
	--ntasks-per-node=1 \
	--cpus-per-task=1 \
	--mem=110G \
	--time=24:00:00 \
	0000.job1.tag.exp.stdy.to.cmbTrspts.sb
