# TEbenchmarking

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18772673.svg)](https://doi.org/10.5281/zenodo.18772673)

TEbenchmarking contains a Snakemake workflow for benchmarking locus-specific transposable-element (TE) quantification from single-cell RNA-seq data. It includes workflows for STARsolo, SoloTE and Stellarscope, together with simulation scripts and R notebooks used to evaluate their results.

- Preprint: [10.64898/2026.02.26.708244](https://doi.org/10.64898/2026.02.26.708244)
- Reproducibility data: [Zenodo record 18772673](https://doi.org/10.5281/zenodo.18772673)


Large FASTQ files, reference genomes, annotations and generated results are not stored in GitHub. This README explains how to obtain and arrange those files, configure storage locations, run the workflow, and add datasets or tools.

## Citation

Please cite the preprint, this repository, and the Zenodo dataset when reusing the workflow or benchmark data. Citation metadata are also provided in `CITATION.cff`.

## License

The repository is distributed under the MIT License. External tools, datasets, annotations, and reference genomes remain subject to their own licenses and terms of use.

## Repository contents

```text
TEbenchmarking/
├── config/
│   └── config.yml                 # Paths, datasets, genomes and parameters
├── workflow/
│   ├── Snakefile                  # Workflow entry point and requested targets
│   ├── rules/                     # Rules for alignment and TE quantification
│   └── envs/                      # Per-rule Conda environments
├── code/                          # External tools included as Git submodules
├── simulations/                   # Scripts used to generate simulated data
├── evaluation/                    # Benchmarking and evaluation notebooks
├── Dockerfile
└── README.md
```

## Data simulation and evaluation

Simulation scripts based on Minnow and Splatter are in `simulations/`. Evaluation notebooks and supporting scripts are in `evaluation/`. These components are separate from the main preprocessing and TE-quantification run described below.


## Workflow overview

The rules are divided across:

| File | Purpose |
| --- | --- |
| `workflow/rules/common.smk` | Configuration lookup and shared helper functions |
| `workflow/rules/star_align.smk` | Reference preparation, STAR indexing and STARsolo alignment |
| `workflow/rules/soloTE.smk` | SoloTE setup and quantification |
| `workflow/rules/starSolo_TE.smk` | STARsolo-based TE quantification |
| `workflow/rules/stellarscope.smk` | Stellarscope preprocessing and assignment |

he rule-level dependency graph is available at
[`workflow/schemas/rulegraph.svg`](workflow/schemas/rulegraph.svg).

## Requirements

The workflow is designed for Linux and has been configured with Snakemake 8.11.6. The minimum version declared in the Snakefile is 5.3.0, but using the tested version is recommended.

The rules declare substantial compute requirements. The largest current rules request up to 24 threads and 250 GB RAM. 

## 1. Clone the repository

Clone recursively because SoloTE and modified SoloTE code are Git submodules:

```bash
git clone --recurse-submodules https://github.com/ScialdoneLab/TEbenchmarking.git
cd TEbenchmarking
```

If the repository was cloned without submodules, initialize them afterward:

```bash
git submodule update --init --recursive
```

## 2. Choose data and result locations

Edit the `paths` section of `config/config.yml`:

```yaml
paths:
  main: "."
  data: "snakemake_data"
  results:
    star: "."
    starsolo_te: "."
    solote: "."
    stellarscope: "."
```

`paths.data` is the root containing all external inputs: FASTQs, reference genomes, gene annotations, TE annotations, barcode whitelists and optional cell annotations.

Each entry under `paths.results` is the root used by the corresponding workflow component:

| Key | Contents |
| --- | --- |
| `star` | STAR indexes and standard STARsolo alignments |
| `starsolo_te` | STARsolo TE and gene-plus-TE results |
| `solote` | SoloTE results |
| `stellarscope` | Stellarscope results |

With the default values, inputs are read from `./snakemake_data`, while indexes and results are written inside the repository.

To use separate mounted volumes, for example:

```yaml
paths:
  main: "/mnt/project/TEbenchmarking"
  data: "/mnt/data/TEbenchmarking-inputs"
  results:
    star: "/mnt/results/star"
    starsolo_te: "/mnt/results/starsolo-te"
    solote: "/mnt/results/solote"
    stellarscope: "/mnt/results/stellarscope"
```

The workflow appends `indexes/` or `results/` to these result roots.

Relative paths are interpreted from the directory in which Snakemake runs. Run all commands below from the repository root.

## 3. Obtain the data from Zenodo

The [Zenodo dataset](https://doi.org/10.5281/zenodo.18772673) contains four archives:

| Archive | Contents |
| --- | --- |
| `data.tar.gz` | Simulated count matrices and corresponding FASTQ files |
| `TEannotations.tar.gz` | Mouse and human TE annotations |
| `whitelists.tar.gz` | Barcode whitelists and cell annotations |
| `jc_mm10.txt.tar.gz` | Mouse TE annotation with Jukes–Cantor distances |

Extract the archives into the directory configured as `paths.data`.

## 4. Prepare the input directory

The complete input layout is:

```text
<data-root>/
├── data/
│   └── <dataset>/
│       ├── fastqs/
│       │   └── <sample_id>/
│       │       ├── ...R1...fastq.gz
│       │       └── ...R2...fastq.gz
│       └── annotation/
│           └── <sample_id>/
│               └── celltype_tsv.tsv      # Only for Stellarscope cell-type pooling
├── references/
│   ├── hg38/
│   │   ├── <human-genome>.fa.gz
│   │   └── <human-gene-annotation>.gtf.gz
│   └── mm10/
│       ├── <mouse-genome>.fa.gz
│       └── <mouse-gene-annotation>.gtf.gz
├── TEannotation/
│   ├── hg38_GENCODE_rmsk_TE.gtf.gz
│   └── mm10_GENCODE_rmsk_TE.gtf.gz
└── whitelists/
    ├── 10xPBMC_sub_whitelist.tsv.gz
    ├── 3M-february-2018.txt.gz
    └── 737K-august-2016.txt.gz
```

FASTQ filenames must contain `R1` or `R2`, because the current alignment rules use those strings to identify barcode/UMI and cDNA reads.


### Reference genomes and gene annotations

Reference FASTA and GTF files are not stored in Git. Obtain matching FASTA/GTF pairs from GENCODE or another documented source.

The distributed configuration uses:

```yaml
genomes:
  mm10:
    fasta: "references/mm10/GRCm38.p4.genome.fa.gz"
    gtf: "references/mm10/gencode.vM10.annotation.gtf.gz"
    rmsk_gtf: "TEannotation/mm10_GENCODE_rmsk_TE.gtf.gz"
  hg38:
    fasta: "references/hg38/GRCh38.p14.genome.fa.gz"
    gtf: "references/hg38/gencode.v30.annotation.gtf.gz"
    rmsk_gtf: "TEannotation/hg38_GENCODE_rmsk_TE.gtf.gz"
```

### Real datasets

Raw data for real datasets are not redistributed through GitHub. Place downloaded or locally available paired FASTQs in:

```text
<data-root>/data/<dataset>/fastqs/<sample_id>/
```

Our config file records external download information where available, including the 10x Genomics URL or SRA accession. 


## 5. Configure datasets

Each dataset is defined under `datasets` in `config/config.yml`:

```yaml
datasets:
  my_dataset:
    sample_ids:
      - sample_1
      - sample_2
    genome: "hg38"
    read_length: 98
    UMI_length: "10"
    whitelist: "whitelists/my_whitelist.txt.gz"
    strandedness: "Forward"
    simulated: false
```

The fields are:

| Field | Meaning |
| --- | --- |
| `sample_ids` | Sample directory names used in input and output paths |
| `genome` | Key under the top-level `genomes` section |
| `read_length` | cDNA read length used to set the STAR index overhang |
| `UMI_length` | Value passed to STARsolo as `--soloUMIlen`; some existing entries also append platform-specific STARsolo options |
| `whitelist` | Barcode whitelist used by STARsolo and Stellarscope |
| `strandedness` | `Forward`, `Reverse`, or `Unstranded` |
| `simulated` | YAML Boolean: `true` for simulated data and `false` for real data |

## 6. Run the workflow

Run all requested benchmark targets with:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/config.yml \
  --software-deployment-method conda \
  --cores 24 \
  --printshellcmds \
  --rerun-incomplete 
```

Change `--cores` to match the available machine. 

To run only one result, request its final path rather than `allout`. For example:

```bash
snakemake \
  -s workflow/Snakefile \
  --configfile config/config.yml \
  --software-deployment-method conda \
  --cores 12 \
  results/SoloTEout/10xPBMC/pbmc8k
```

## 7. Expected outputs

With all result roots set to `.`, the main outputs are organized as follows:

```text
TEbenchmarking/
├── indexes/
│   └── genomeIndexes_<genome>_<read_length>/
├── results/
│   ├── STARoutdir/
│   │   └── <dataset>/<sample_id>/
│   │       ├── best_Aligned.sortedByCoord.out.bam
│   │       ├── best_Solo.out/
│   │       ├── score_Aligned.sortedByCoord.out.bam
│   │       └── score_Solo.out/
│   ├── SoloTEout/
│   │   └── <dataset>/<sample_id>/
│   ├── SoloTEout_thr0/
│   ├── SoloTEout_thr1/
│   ├── SoloTEout_thr2/
│   ├── stellarscope_out/
│   │   └── <dataset>/<sample_id>/
│   │       ├── stload/
│   │       ├── pseudobulk/
│   │       ├── pseudobulk_thr0.9/
│   │       └── celltype/
│   ├── STARsolo_EM_TE/
│   │   └── <dataset>/<sample_id>/TE_Solo.out/
│   └── STARsolo_EM_TEwGenes/
│       └── simulated_mm_RA_wGenes/all/TE_Solo.out/
├── logs/
├── benchmarks/
└── tmp/
```

If different result roots are configured, each tool creates its relevant part of this structure under its own root. Logs, benchmarks, temporary files, and workflow diagrams currently remain relative to the repository working directory.

## 8. Add a new dataset

1. Create one FASTQ directory per sample:

   ```text
   <data-root>/data/<dataset>/fastqs/<sample_id>/
   ```

2. Add the dataset and all sample IDs to `config/config.yml`.
3. Select an existing genome key or add a new entry under `genomes`.
4. Provide a compatible barcode whitelist.
5. Set `read_length`, `UMI_length`, `strandedness`, and `simulated`.
6. For Stellarscope cell-type pooling, add:

   ```text
   <data-root>/data/<dataset>/annotation/<sample_id>/celltype_tsv.tsv
   ```

7. Dry-run one final target before adding the dataset to a large aggregation run.

Example:

```bash
snakemake \
  -s workflow/Snakefile \
  --software-deployment-method conda \
  --cores 1 \
  --dry-run \
  results/SoloTEout/my_dataset/sample_1
```

## 9. Add another TE-analysis tool

To extend the benchmark with a new tool:

1. add a rule module such as `workflow/rules/new_tool.smk`;
2. add a pinned environment under `workflow/envs/`;
3. include the rule module from `workflow/Snakefile`;
4. add a configurable output root under `paths.results`;
5. use `{dataset}` and `{sample_id}` wildcards so the tool works with existing datasets;
6. declare threads, memory, logs, and benchmark files;
7. add only the final outputs to the aggregation rule;


