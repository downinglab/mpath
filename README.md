# MPATH
Methylation pseudotime analysis for deciphering epigenetic cross-talk across sub-cell-cycle timescales

## Installation & General Usage

1. For Python metrics code, prepare a Python 3 environment and install joblib, pandas, numpy, scipy, and tqdm. For development, Python 3.10 was used.

2. Clone the project using:

```bash
git clone https://github.com/downinglab/mpath.git
```

3. Run Python code from the the command line.

4. Run Matlab PCA code using output CSV from Python.

## Python Usage

MPATH calculates pseudotime for long read data from Oxford Nanopore Technologies sequencing platforms. The flow of information is as follows:

Nanopore BAM file with MM/ML tags -> modkit -> CpG methylation bed file -> Python metrics code -> MATLAB PCA code -> Pseudotime score

In the future, the PCA Matlab code may be incorporated as Python, and MPATH distributed via pip.

### Python Input BED File Structure

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

### Python Command Line Arguments

The other command line arguments to the code are as follows:

| name                  | description                                         | type  | required | default value         |
|-----------------------|-----------------------------------------------------|-------|----------|-----------------------|
| path_input_bed        | path of input bed file                              | str   | yes      | na                    |
| path_output_csv       | path of output csv file                             | str   | yes      | na                    |
| min_cpgs              | minimum number of cpgs on a read to calculate stats | int   | no       | 3                     |
| bin_limits            | comma-separated limits of distance bins             | str   | no       | 0,100,1000,5000,10000 | 
| batch_size            | number of cpgs to process in a batch (parallel)     | int   | no       | 1e8                   |
| p                     | number of processes (parallel)                      | int   | no       | 1                     |
| use_full_matrix       | use the full matrix to make calculations            | bool  | no       | False                 |

An example call might look like this:

```bash
python ./metrics/longread_metrics_parallel.py -path_input_bed ./data/input.bed -path_output_csv ./data/output.csv -p 8 -min_cpgs 3 -bin_limits 0,100,1000,5000,10000 --use_full_matrix
```

### Command Line Notes

For large input bed files, it is recommended to use the parallel version, as the other version tries to read the entire file into memory.

bin_limits - A comma-separated string of the limits of distance bins, which limits which CpGs are used to calculate metrics. An input of 0,100,1000,5000,10000 will generate the following distance bins: [[0,100], [100,1000], [1000,5000], [5000,10000]]
batch_size - The number of reads to use in a batch. Adjust as needed to avoid using too much RAM.
p - The number of processes to run in parallel. More processes in parallel will use more RAM (multiplicative with batch_size). Adjust as needed.
use_full_matrix - The Python script works by generating matrices representing all CpG pairs in the read. By default, only one corner of the matrix is needed, but the full matrix can also be used. It may slightly change values of output metrics.

### Output CSV column descriptions

The Python code outputs a CSV in long format, where each row corresponds to a particular read. Since multiple distance bins are generated

| column | name                  | description                                               | type  |
|--------|-----------------------|-----------------------------------------------------------|-------|
| 1      | read_id               | read id in longread bam                                   | str   |
| 2      | read_wgbs_distance    | difference of longread and wgbs methylation ratio         | float |
| 3      | read_meth_ratio       | total methylation ratio of read (methylated:total cpgs)   | float |
| 4      | distance_bin          | which distance bin this row is for                        | str   |
| 5      | num_pairs             | number of cpg pairs in this read                          | int   |
| 6      | smc                   | simple matching coefficient                               | float |
| 7      | uniformity            | uniformity score                                          | float |
| 8      | pearson_r             | pearson rho                                               | float |
| 9      | pearson_p             | pearson p value                                           | float |

## MATLAB PCA Code Instructions

After generation of longread_metrics.py output CSVs, the files may be used in NascMatur_PCA_scatterhistogram.m for PCA and visualization. 

### Requirements

An updated version of MATLAB with the Statistics and Machine Learning Toolbox installed.

### Run instructions

1. Prepare two data tables as described in input file guidelines.
1. Open the script NascMatur_PCA_scatterhistogram.m in the MATLAB environment.
1. Run the script.
1. Select the primary ("nascent") data table.
1. Select the secondary ("mature") data table.

### Input file guidelines

The NascMatur_PCA_scatterhistogram.m script calls for one nascent and one mature long read data tables with consistent formatting between the two datasets. If the input files are directly from the Python metric script, no further action is needed. Otherwise, the first column must be of string type and each row must be one uniquely labeled read ID. The remaining N columns must be of numerical types; the data in these N columns will be used for PCA. A specific naming scheme for column variable names is not required, but the metric categories (i.e., the type of calculation) of each column must be the same between the two datasets. Therefore, the number of columns must be the same between the two data tables while the number of rows or unique reads need not be equal. An example of two data tables may be found below:

#### Nascent data table

|read_id|metric_1|metric_2|...|metric_N|
|-------|--------|--------|---|--------|
|read_1 |num_type|num_type|  .|num_type|
|read_2 |num_type|num_type|  .|num_type|
|    ...|     ...|     ...|  .|	    ...|
|read_M1|num_type|num_type|  .|num_type|

#### Mature data table

|read_id|metric_1|metric_2|...|metric_N|
|-------|--------|--------|---|--------|
|read_x |num_type|num_type|  .|num_type|
|read_y |num_type|num_type|  .|num_type|
|    ...|     ...|     ...|  .|	    ...|
|read_M2|num_type|num_type|  .|num_type|