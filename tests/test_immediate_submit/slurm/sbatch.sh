#!/bin/bash
set -e && echo "$0 $*" >&2

set -vx

sbatch_params=$1
script=$2
dependencies=$3
flag=tmp/slurm-signal/dependencies-$(basename $script)-finish

mkdir -p tmp/slurm-signal

echo "#!/bin/bash"                      > $flag.sh
echo ''                                 >> $flag.sh
echo "for i in ${dependencies//,/ }"    >> $flag.sh
echo 'do'                               >> $flag.sh
echo '    while [ ! -f $i ]'            >> $flag.sh
echo '    do'                           >> $flag.sh
echo '        sleep 3'                  >> $flag.sh
echo '    done'                         >> $flag.sh
echo 'done'                             >> $flag.sh
echo ''                                 >> $flag.sh
cat $script                             >> $flag.sh
echo ''                                 >> $flag.sh
echo "touch $flag"                      >> $flag.sh

nohup bash $flag.sh > $sbatch_params.log 2>&1 &

echo $flag
