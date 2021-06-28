#!/bin/bash
# Simple virus model simulation. By Assaf Amitai, June 2021, amitaiassaf@gmail.com

for SpikeNum in 40;  do
for ((cond=1; cond<=7; cond++)); do
for ((idumNum=402; idumNum<=402; idumNum++)); do
for ((epitope=7; epitope<=234; epitope++)); do
    echo ${idumNum}
    echo SpikeNum ${SpikeNum}
    echo ${epitope}
    echo $cond
    JOB=$(qsub -q forever -v NUMBERARG=$idumNum,EPITOPEID=$epitope,COND=$cond,SN=$SpikeNum  lammps_flu.pbs)

#STRIPSTRING='.lc2'
#N=${PBS_JOBID%%$STRIPSTRING}
#echo NN
#echo ${JOB}
#echo M
#"$PBS_JOBID"

#"$PBS_O_WORKDIR"
done
done
done
done