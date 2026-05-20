# epimut-load

`epimut-load` is a lightweight Python command-line tool to detect Stochastic Epigenetic Mutations (SEMs) and compute Epigenetic Mutation Load (EML) from DNA methylation matrices.

The tool is designed to be simple, reproducible, and easy to configure through a YAML file.

## Features

- Detect stochastic epigenetic mutations using an IQR-based approach
- Compute sample-level Epigenetic Mutation Load
- Support internal reference cohorts
- Flexible metadata-based sample filtering
- Compute SEM direction and SEM strength
- Export SEM matrices, SEM annotations, and sample-level metrics
- Simple command-line interface

## Installation

Clone the repository:

```bash
git clone https://github.com/Jackcava/epimut-load.git
cd epimut-load
pip install -e .
```

## Usage
Run the analysis with:

```
epimut-load run --config configs/example_config.yaml
```

## Input files
### DNA methylation matrix
Rows should represent CpG sites and columns should represent samples.
Example:

CpG        |   S01 |   S02 |  S03
---------- | ----: | ----: | ---:
cg00000001 |  0.12 |  0.18 | 0.15
cg00000002 | -1.20 | -1.10 | 2.90
.......... | ..... | ..... | ....


### Metadata file
The metadata file must contain one column matching the sample IDs in the methylation matrix.
Example:

sampleID | group   |  y
-------- | ------- | -:
S01      | control |  0
S02      | control |  0
S03      | case    |  1
........ | ....... | ..


## Configuration
Example configuration:
```yaml
input:
  path: data/input
  methylation_matrix: example_methylation.csv
  metadata: example_metadata.csv
  sample_id_col: sampleID


output:
  path: data/output
  name_exp: example


data:
  matrix_type: M_values
  cpg_id_col: CpG


# reference
reference:
  mode: internal # [internal | external]

  internal:
    include_filter: "group == 'control'"
    exclude_filter: "" # e.g. "y.isna()"

  external:
    path: ""
    methylation_matrix: ""
    metadata: ""
    include_filter: "" # e.g.  "age < 10"
    exclude_filter: "" # e.g. "y.isna()"


# target
target_samples:
  include_filter: ""
  exclude_filter: ""


# SEMs calculation settings
sem_calling:
  iqr_multiplier: 3
  direction: true
  compute_strength: true
```

## Method

For each CpG site, SEMs are detected using the distribution of the selected reference samples.

A CpG value is classified as a SEM if it falls outside:
```
Q1 - k * IQR
Q3 + k * IQR
```
where k is controlled by iqr_multiplier, with default value 3.

Values below the lower threshold are classified as hypomethylated SEMs.
Values above the upper threshold are classified as hypermethylated SEMs.

The Epigenetic Mutation Load is computed as the total number of SEMs detected in each sample.

## Output files

The analysis produces:
```
outliers_<name_exp>.tsv
sem_strength_<name_exp>.tsv
sems_<name_exp>.tsv
sample_file_metrics_<name_exp>.tsv
```

| File                        | Description                                        |
| --------------------------- | -------------------------------------------------- |
| `outliers_*.tsv`            | CpG-by-sample matrix with values `-1`, `0`, `1`    |
| `sem_strength_*.tsv`        | CpG-by-sample matrix containing SEM strength       |
| `sems_*.tsv`                | List of SEMs, hyper-SEMs, and hypo-SEMs per sample |
| `sample_file_metrics_*.tsv` | Metadata annotated with EML and SEM-level metrics  |


## Example data

A small simulated dataset is provided in:
```
data/input
```

It can be used to test the package:
```bash
epimut-load run --config configs/example_config.yaml
```


## References

- Gentilini D, Garagnani P, Pisoni S, et al. Stochastic epigenetic mutations (DNA methylation) increase exponentially in human aging and correlate with X chromosome inactivation skewing in females. Aging (Albany NY). 2015;7(8):568-578. doi:10.18632/aging.100792

- Yan Q, Paul KC, Lu AT, et al. Epigenetic mutation load is weakly correlated with epigenetic age acceleration. Aging (Albany NY). 2020;12(18):17863-17894. doi:10.18632/aging.103950

## Example application

This approach has also been applied in:

- Cavalca G, Vergani M, Cangelosi D, et al. Stochastic epigenetic mutation profiles as biomarkers of clinical activity in juvenile idiopathic arthritis: a multi-omic machine learning approach for gene prioritization. Molecular Medicine. 2025;31:289. doi:10.1186/s10020-025-01348-6



## License

This project is released under the MIT License.