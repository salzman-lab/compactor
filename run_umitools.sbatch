#!/bin/bash
#
#SBATCH --job-name=run
#SBATCH --output=job_output/run.%j.out
#SBATCH --error=job_output/run.%j.err
#SBATCH --time=24:00:00
#SBATCH -p horence,quake,owners
#SBATCH --nodes=1
#SBATCH --mem=60Gb
#SBATCH --requeue

dataset=$(basename ${PWD})
srp=$(echo ${dataset} | cut -f2 -d_)

ids=$(pysradb srp-to-srr ${srp} | tail -n +2 | cut -f2 -d' ')

for id in ${ids}
do 
    sbatch run_one_umitools.sbatch ${id}
done

