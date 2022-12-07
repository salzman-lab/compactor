#!/bin/bash
#
#SBATCH --job-name=umi
#SBATCH --output=job_output/run.%j.out
#SBATCH --error=job_output/run.%j.err
#SBATCH --time=24:00:00
#SBATCH -p horence,quake,owners
#SBATCH --nodes=1
#SBATCH --mem=60Gb
#SBATCH --requeue

## For Sherlock users, try activating George Henderson's virtual environment, or otherwise be sure umi_tools is available. 
source /oak/stanford/groups/horence/george/dog/bin/activate

id=$1

R1=${id}_1.fastq.gz
R2=${id}_2.fastq.gz
extracted_R1=${id}_1_extracted.fastq.gz
extracted_R2=${id}_2_extracted.fastq.gz

# run UMItools on FASTQ files
umi_tools whitelist \
    --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
    --stdin ${R1} \
    -L whitelist.${id}.log \
    > whitelist_${id}.txt
umi_tools extract \
    --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
    --stdin ${R1} \
    --stdout ${extracted_R1} \
    --read2-in ${R2} \
    --read2-out ${extracted_R2} \
    --whitelist=whitelist_${id}.txt

echo -e "Done with ${id}"

