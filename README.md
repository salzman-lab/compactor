# Compactor Generation
Implementation of compactors used in Tabula Sapiens Smart-Seq 2 analysis. 

Compactors are seed-based contigs constrained to FASTQ read-length. NOMAD-called anchors are taken as input and >= 1 compactors are generated to represent all sequence diversity downstream of the anchor.


### Usage

This script takes as input the results of the Salzman Lab's NextFlow implementation of NOMAD. In aspen.sh, the user specifies:
1. A path to the samplesheet used input to NOMAD.
2. A path to the NOMAD results directory.

### Tuning

submit_parse.py is used to parallize jobs parsing samplesheet FASTQs for reads containing NOMAD-called anchors. In cases where the user has few, but very large FASTQs (10s of GB), it is recommended to set fastqs_to_process_in_parallel to 1. In cases where the user has many FASTQs (100s, as in the case with Tabula Sapiens SS2 data), we recommend setting fastqs_to_process_in_parallel to 10. 

parallel_generate.py is used to parallelize jobs generating compactors from anchor-specific intermediate files (sets of reads containing an anchor). We recommend that the user set anchors_to_process_in_parallel to a value such that  {total anchors input to compactor generation} /  anchors_to_process_in_parallel < 1500.

