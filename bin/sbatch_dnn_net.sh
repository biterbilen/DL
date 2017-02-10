#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_dnn_net.sh
# 
#         USAGE: ./sbatch_dnn_net.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Biter Bilen (), biterbilen@yahoo.com
#  ORGANIZATION: 
#       CREATED: 01/01/2017 03:28:50 PM
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
source ~/Projects/biter_biter_shared/bin/auxslurm.sh

# label stats
DNNib() {
	fle=$1
	echo FUNCTION $(date):  $FUNCNAME $*
	samtools index $fle
}

DNNpi() {
	name=$1
	echo FUNCTION $(date):  $FUNCNAME $*
	source ~/PI_HOME/PythonSandbox/theano/bin/activate; 
	export THEANO_FLAGS=mode=FAST_RUN,device=cpu
	python -c "import sys; sys.path.append('/home/biter/Projects/biter_biter_stemCellFateRegulation/bin'); from dnn_net.datasets.ENCODEDREAMchallengeTFBS import load_data; load_data('$name', force_read=False, verbose=2)"
}

DNN() {
	local analysisdirbase=$1

	ln -sf ${o[datadir]} .
	local indir=${o[datadir]##*/}
	local indir=${o[datadir]##*/}

	local qcmd
	local jobname=pi
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then 
		names=($(ls $indir/labels/*.train.labels.tsv.gz))
		names=(${names[@]##*/})
		names=(${names[@]/.train.labels.tsv.gz/})
		qcmd=$(cat <<- END
			names=(${names[@]})
			name=\${names[\$${o[stiv]}]}
			$(declare -f DNNpi)
			DNNpi \$name
		END
		)
		jobnames+=($jobname$analysisdirbase)
		qcmds[$jobname$analysisdirbase]="$qcmd"
		qarrs[$jobname$analysisdirbase]=$(( ${#names[@]} - 1 ))
		#qarrs[$jobname$analysisdirbase]=0
		qwpts[$jobname$analysisdirbase]=F
	fi

}

slurmMain DNN --jobnames=pi --datadir=/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484 --t=6 --tmin=20 --nc=4 --queue=normal


