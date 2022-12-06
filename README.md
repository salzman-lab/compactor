# Compactor Generation
The implementation of compactors used in Tabula Sapiens Smart-Seq 2 analysis. 

Compactors are seed-based contigs constrained to FASTQ read-length. NOMAD-called anchors are taken as input and >= 1 compactors are generated to represent all sequence diversity downstream of the anchor.


## Usage

This script takes as input the results of the Salzman Lab's NextFlow implementation of NOMAD. In aspen.sh, the user specifies:
1. A path to the samplesheet used input to NOMAD.
2. A path to the NOMAD results directory.
Pull this repository to your machine, and call 'aspen.sh' with a single argument to specify the run name: we might do
sbatch aspen.sh 'test_compactors'.

## Dependencies

This code requires the following: 
1. Python 3.9
2. NumPy 1.20.3
3. Pandas 1.3.1 
4. BioPython 1.79


## Parallelization tuning

submit_parse.py is used to parallize jobs parsing samplesheet FASTQs for reads containing NOMAD-called anchors. In cases where the user has a few very large FASTQs (10s of GB), it is recommended to set fastqs_to_process_in_parallel to 1. In cases where the user has many FASTQs (100s, as in the case with Tabula Sapiens SS2 data), it is recommended to set fastqs_to_process_in_parallel to 10. This parameter can be found on line 10 of submit_parse.py. 

parallel_generate.py is used to parallelize jobs generating compactors from anchor-specific intermediate files (sets of reads containing an anchor). We recommend that the user set anchors_to_process_in_parallel to a value such that  {total anchors input to compactor generation} /  anchors_to_process_in_parallel < 1500. This parameter can be found on line 10 of parallel_generate.py. 

## FASTQ downsampling

make_intermediates.py uses 2 parameters to control the amount of reads we collect per anchor in each FASTQ. By default, we take 200 anchor-reads per FASTQ; this default can be modified on line 10 of make_intermediates.py. By default, we also permit 1 million reads to be selected per FASTQ; this can be modified on line 30 of make_intermediates.py. 
