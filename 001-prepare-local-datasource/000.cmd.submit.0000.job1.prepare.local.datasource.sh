export currDir=`pwd`
export taxonId=`pwd | awk -F '/' '{print $(NF-1)}'`
ln -sf $currDir/../../00-config/commonEnvironmentVars.cfg $currDir/comEnvirVars
ln -sf $currDir/00000.parameter.of.0000.job1.cfg $currDir/speEnvirVars

sbatch --job-name="001_"$taxonId \
	--export=commonEnvir=$currDir/comEnvirVars,speciesEnvir=$currDir/speEnvirVars,threadNum=18 \
	--nodes=1 \
	--ntasks-per-node=1 \
	--cpus-per-task=18 \
	--mem=110G \
	--time=03:50:00 \
	0000.job1.prepare.local.datasource.sb
