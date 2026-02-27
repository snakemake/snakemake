#!/bin/bash
echo `date` >> qsub.log
tail -n1 $1 >> qsub.log
# simulate printing of job id by a random number
echo $RANDOM
cat $1 >> qsub.log
sh $1
