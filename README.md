# NanoQMRA

NanoQMRA is a bioinformatics workflow for species-resolved microbial risk scoring from nanopore long reads. For each read it combines:

- pathogen / WHO-priority association
- plasmid mobility (conjugative / mobile / non-mobile)
- antibiotic resistance gene (ARG) carriage

Outputs:

- a per-read scored table
- a ranked per-read risk list
- a sample-level `risk_score`

Per-read score:

```text
risk_score = WHO_priority_score × ARG_rank_score × mobile_factor × conjugative_plasmid_factor
```

Sample score:

```text
sample risk_score = sum(read risk_score) / number of FASTA records
```

This is a heuristic resistome-risk index, not a classical dose–response QMRA (no exposure or infection probability).

## Installation

`environment.yml` lists the conda environment (Python and pandas). After the env exists, install this repository as the `nanoqmra` command.

```bash
git clone https://github.com/sustc-xylab/NanoQMRA.git
cd NanoQMRA
conda env create -f environment.yml
conda activate nanoqmra
pip install .
nanoqmra check
```

Scoring does **not** need ARGpore2 databases.

### Full workflow dependencies

`nanoqmra run` also needs two external tools in **their own conda envs**. The runner calls them with `conda run` (you do not activate those envs yourself):

| Dependency | Conda env | Role | How to provide |
| --- | --- | --- | --- |
| [ARGpore2](https://github.com/sustc-xylab/ARGpore2) | `argpore2` + `argpore2-r` | ARG + taxonomy on nanopore reads | Clone and set up databases. `argpore2`: seqkit, blast, centrifuge, parallel, perl-bioperl. R packages cannot share that env (solver conflict), so they live in `argpore2-r` (`r-base`, plyr, data.table, dplyr, foreach, doParallel). `run` stitches both via `conda run`. |
| [Plascad](https://github.com/pianpianyouche/plascad) 1.17 | `Plascad` | Plasmid mobility classification | `conda create -n Plascad -y -c bioconda -c pianpianyouche plascad` |

Check what is available:

```bash
nanoqmra check --argpore-home /path/to/ARGpore2 --argpore-env argpore2 --plascad-env Plascad
```

## Usage

One nanopore FASTA per run. `.fa` is accepted; the runner copies it to `.fasta` because Plascad requires that suffix. Avoid spaces in file names.

### Score only

Use this when ARGpore2 and Plascad tables already exist.

```bash
nanoqmra score \
  --arg-taxa sample.arg.w.taxa.tab \
  --plasmid-sum sample_plasmids_classification_sum.txt \
  -o results \
  --prefix sample \
  --total-reads 2407
```

### Full workflow

Writes **only** under `-o`. Nothing is copied back into the git repository.

```bash
nanoqmra run \
  -f examples/example.fa \
  -t 40 \
  -o /path/to/outdir \
  --argpore-home /path/to/ARGpore2 \
  --argpore-env argpore2 \
  --plascad-env Plascad
```

Shell wrapper (same options):

```bash
bash NanoQMRA.sh -f examples/example.fa -t 40 -o /path/to/outdir \
  --argpore-home /path/to/ARGpore2 -p Plascad
```

## Command-line arguments

### `nanoqmra score`

| Argument | Required | Description |
| --- | :---: | --- |
| `--arg-taxa` | yes | ARGpore2 `*.arg.w.taxa.tab` |
| `--plasmid-sum` | yes | Plascad `*_plasmids_classification_sum.txt` |
| `-o, --outdir` | yes | Output directory |
| `--prefix` | no | Output file prefix (default: from `--arg-taxa` name) |
| `--rules` | no | Rules directory (default: packaged `src/nanoqmra/data/rules`) |
| `--arg-csv` | no | ARG rank CSV (default: packaged `ARG.csv`) |
| `--total-reads` | no | FASTA record count; if set, appends sample `risk_score` |

### `nanoqmra run`

| Argument | Required | Description |
| --- | :---: | --- |
| `-f, --fasta` | yes | Input nanopore FASTA (`.fa` or `.fasta`) |
| `-t, --threads` | yes | Threads for ARGpore2 |
| `-o, --outdir` | yes | Output directory |
| `--argpore-home` | yes* | ARGpore2 install directory (`argpore.sh` lives here). Also reads `ARGPORE_HOME` |
| `--argpore-env` | no | Conda env for ARGpore2 tools (default: `argpore2`, or `ARGPORE_ENV`) |
| `--plascad-env` | no | Conda env that provides `Plascad` (default: `Plascad`) |

\* Required unless `ARGPORE_HOME` is set.

### `nanoqmra check`

| Argument | Required | Description |
| --- | :---: | --- |
| `--argpore-home` | no | Directory to inspect for ARGpore2 + databases |
| `--argpore-env` | no | Conda env to test for perl / parallel / R (default: `argpore2`) |
| `--plascad-env` | no | Conda env name to test (default: `Plascad`) |

`check` exits 0 if Python scoring data are present. Missing ARGpore2 databases are reported as warnings; they only block `run`.

## Output layout

```text
outdir/
├── input/                      # normalized *.fasta
├── logs/
│   ├── *_ARGpore.log
│   └── *_Plascad.log
├── intermediate/
│   ├── argpore2/
│   └── plascad/
└── results/
    ├── *_mother_table.scored.tsv
    └── *_risk_score_list.tsv
```

`*_mother_table.scored.tsv` — ARG, WHO priority, mobility factors, and per-read `risk_score`.

`*_risk_score_list.tsv` — ranked reads; sample-level `risk_score` is appended when `--total-reads` is set (always for `run`).

Packaged reference data (do not duplicate at repo root):

- `src/nanoqmra/data/rules/` — WHO matching rules
- `src/nanoqmra/data/reference/ARG.csv` — ARG rank table

## License

No license file has been added yet. The GitHub default is “no assertion”; choose MIT / Apache-2.0 / GPL-3.0 before a public release.

## Citation

See `CITATION.cff`. Related papers:

Y. CHE, et al., Conjugative plasmids interact with insertion sequences to shape the horizontal transfer of antimicrobial resistance genes. PNAS. 118(6):e2008731118. doi: 10.1073/pnas.2008731118

Ziqi Wu, You Che, Chenyuan Dang, Miao Zhang, Xuyang Zhang, Yuhong Sun, Xiang Li, Tong Zhang, Yu Xia*. 2022. Nanopore-based Long-read Metagenomics Uncover the Resistome Intrusion by Antibiotic Resistant Bacteria from Treated Wastewater in Receiving Water Body. Water Research. 226:119282. (https://doi.org/10.1016/j.watres.2022.119282)

Yuxi Yan, Jiao Yang, Miao Zhang, Bixi Zhao, Cailong Nie, Baobo Liu, Chen Wang, Yongsheng Chen, Yu Xia*. A Species-resolved Quantitative Microbial Risk (QMRA) Framework for Hospital Wastewater Based on Nanopore Long-Read Sequencing.
