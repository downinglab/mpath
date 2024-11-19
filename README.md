# MPATH
Methylation pseudotime analysis for deciphering epigenetic cross-talk across sub-cell-cycle timescales

## Installation

For Python metrics code, run using a Python 3 environment with joblib, pandas, numpy, scipy, and tqdm installed

```bash
git clone https://github.com/downinglab/mpath.git
```

## Usage

MPATH calculates pseudotime for long read data from Oxford Nanopore Technologies sequencing platforms.

### Constructing input bed file

A bed file serves as the input to longread_metrics_parallel.py or longread_metrics.py. The bed file rows should represent individual CpGs and have the following columns: 

| column | name                  | description                                        | type  |
|--------|-----------------------|----------------------------------------------------|-------|
| 1      | chrom                 | name of the chromosome                             | str   |
| 2      | start position        | start position of cpg                              | int   |
| 3      | end position          | end position of cpg                                | int   |
| 4      | strand                | '+' for positive strand '-' for negative strand    | str   |
| 5      | read id               | read id of where the cpg came from                 | str   |
| 6      | long read methylation | methylated (1) or unmethylated (0) state           | int   |
| 7      | wgbs methylation      | methylation ratio of cpg in wgbs data              | float |

Example script of how to generate input bed file:
```bash
TODO
```

## Description of longread_metrics.py output CSV

### Output CSV column descriptions

| column | name                  | description                                                                    | type  |
|--------|-----------------------|--------------------------------------------------------------------------------|-------|
| 1      | read_id               | read id in longread bam                                                        | str   |
| 2      | read_wgbs_distance    | difference of longread and wgbs methylation ratio                              | float |
| 3      | read_meth_ratio       | total methylation ratio of read (methylated:total cpgs)                        | float |
| 4      | distance_bin          | which distance bin this row is for                                             | str   |
| 5      | num_pairs             | number of cpg pairs in this read                                               | int   |
| 6      | smc                   | simple matching coefficient                                                    | float |
| 7      | uniformity            | uniformity score                                                               | float |
| 8      | pearson_r             | pearson rho                                                                    | float |
| 9      | pearson_p             | pearson p value                                                                | float |
