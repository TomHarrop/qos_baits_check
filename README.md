# qos_baits_check

Check the k-mer content of your bait files, and cluster them at various
identity levels.

## Overview

This is an informational workflow for bait design. The aim is to get a rough
idea of the number of baits that will be synthesised from the target file.

It takes a target file in FASTA format, counts the unique k-mers (e.g. 80-mers
for 80 nt baits), and then clusters the unique k-mers at various identity
levels.


## Workflow

![`snakemake --rulegraph`](assets/graph.svg)
