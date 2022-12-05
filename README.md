# Compactor Generation
Implementation of compactors used in Tabula Sapiens Smart-Seq 2 analysis. 

Compactors are seed-based contigs constrained to FASTQ read-length. NOMAD-called anchors are taken as input and >= 1 compactors are generated to represent all sequence diversity downstream of the anchor.


### Usage

This script takes as input the results of the Salzman Lab's NextFlow implementation of NOMAD. In aspen.sh, the user specifies:
1. A path to the samplesheet used input to NOMAD.
2. A path to the NOMAD results directory.

### In-depth usage.

