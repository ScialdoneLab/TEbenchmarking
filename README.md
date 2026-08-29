# TEbenchmarking

Benchmarking single-cell, locus-specific TE analysis tools in scRNA-seq

## Preprocessing and TE quantification

Snakemake workflow in folders:

- **config** (config file)
- **workflow** (Snakefile, rules and conda environments)

Dockerfile: Dockerfile_containerize

## Configure input and output paths

Machine-specific paths are defined in `config/config.yaml`.

Before running the workflow, edit the `paths` section to identify:

- the repository directory;
- the input-data directory;
- the output directory used by each tool.

Relative paths are resolved from the directory where Snakemake is executed.
Absolute paths can be used for mounted volumes on virtual machines.

Example:

```yaml
paths:
  main: "/mnt/project/TEbenchmarking"
  data: "/mnt/data/snakemake_data"

  results:
    star: "/mnt/results/star"
    starsolo_te: "/mnt/results/starsolo_te"
    solote: "/mnt/results/solote"
    stellarscope: "/mnt/results/stellarscope"

## Data simulation

Scripts for data simulation using Minnow and splatter in folder **simulations**

## Evaluation

Notebooks with analysis used for tool evaluation available in folder **evaluation**


*Preprint published on Biorxiv (10.64898/2026.02.26.708244)*

*Data is be available on Zenodo (10.5281/zenodo.18772673)*
