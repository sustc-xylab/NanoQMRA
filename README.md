# NanoAMRRank

NanoAMRRank is an assembly-free bioinformatics workflow for antimicrobial resistance (AMR) risk scoring and prioritization from nanopore long reads. It links antibiotic resistance genes (ARGs), mobile genetic elements (MGEs), and pathogen-associated taxonomic signals (PATSs) on individual reads, providing read-level evidence for species-resolved AMR risk assessment. For each read it combines:

- pathogen-associated taxonomy and WHO priority
- plasmid mobility (conjugative / mobile / non-mobile)
- antibiotic resistance gene (ARG) carriage

Outputs:

- a per-read scored table
- a ranked per-read risk list
- a sample-level relative AMR risk index (`risk_score`)

Read-level relative risk score:

```text
risk_score = WHO_priority_score × ARG_rank_score × mobile_factor × conjugative_plasmid_factor
```

Sample-level relative risk index:

```text
sample risk_score = sum(read risk_score) / number of FASTA records
```

These scores support comparative AMR risk prioritization. The workflow does not calculate exposure doses or infection probabilities and does not implement pathogen-specific dose–response modeling. It can provide inputs for subsequent quantitative microbial risk assessment (QMRA).

The outputs listed here rank individual reads and summarize the sample-level index. Species-level ranking requires downstream aggregation of read scores by assigned host species. The scoring commands described below do not incorporate iDNA/eDNA weighting, absolute 16S rRNA gene abundance, or questionnaire-derived exposure information.

## Installation

`environment.yml` lists the conda environment (Python and pandas). After the env exists, install this repository as the `nanoamrrank` command.

```bash
git clone https://github.com/sustc-xylab/NanoAMRRank.git
cd NanoAMRRank
conda env create -f environment.yml
conda activate nanoamrrank
pip install .
nanoamrrank check
```

Scoring does **not** need ARGpore2 databases.

### Full workflow dependencies

`nanoamrrank run` also needs two external tools in **their own conda envs**. The runner calls them with `conda run` (you do not activate those envs yourself):

| Dependency | Conda env | Role | How to provide |
| --- | --- | --- | --- |
| [ARGpore2](https://github.com/sustc-xylab/ARGpore2) 2.1.1 | `argpore2` + `argpore2-r` | ARG + taxonomy on nanopore reads | `conda env create -f environment-tools.yml` and `environment-r.yml`, then `pip install -e .`, `argpore2 setup-tools`, and `argpore2 setup-db`. `run` calls `argpore2 -f ... -t ... -o ...` in the `argpore2` environment. |
| [Plascad](https://github.com/pianpianyouche/plascad) 1.17 | `Plascad` | Plasmid mobility classification | `conda create -n Plascad -y -c bioconda -c pianpianyouche plascad` |

Check what is available:

```bash
nanoamrrank check --argpore-home /path/to/ARGpore2 --argpore-env argpore2 --plascad-env Plascad
```

## Usage

One nanopore FASTA per run. Accepted suffixes are `.fa`, `.fasta`, `.fas`, `.fna`, and `.ffn`. Accepting the suffix does not rename the file. Plascad later opens `<stem>.fasta` in its working directory (`Plasmids_plot.py` uses `args.i + ".fasta"`), so the runner hard-links the input to `<prefix>.fasta` under `-o` when the filesystem allows it, and copies only when a hard link cannot be created. Avoid spaces in file names.

### Score only

Use this when ARGpore2 and Plascad tables already exist.

```bash
nanoamrrank score \
  --arg-taxa sample.arg.w.taxa.tab \
  --plasmid-sum sample_plasmids_classification_sum.txt \
  -o results \
  --prefix sample \
  --total-reads 2407
```

### Full workflow

Writes **only** under `-o`. Nothing is copied back into the git repository.

```bash
nanoamrrank run \
  -f examples/example.fa \
  -t 40 \
  -o /path/to/outdir \
  --argpore-home /path/to/ARGpore2 \
  --argpore-env argpore2 \
  --plascad-env Plascad
```

On the reference databases shipped with ARGpore2 (PLSDB `2020_03_04`, SARG `20211207`, Centrifuge `p_compressed+h+v`), `examples/example.fa` (2,407 records) yields a sample-level relative AMR risk index of **0.154965**.

## Command-line arguments

### `nanoamrrank score`

| Argument | Required | Description |
| --- | :---: | --- |
| `--arg-taxa` | yes | ARGpore2 `*.arg.w.taxa.tab` |
| `--plasmid-sum` | yes | Plascad `*_plasmids_classification_sum.txt` |
| `-o, --outdir` | yes | Output directory |
| `--prefix` | no | Output file prefix (default: from `--arg-taxa` name) |
| `--rules` | no | Rules directory (default: packaged `src/nanoamrrank/data/rules`) |
| `--arg-csv` | no | ARG rank CSV (default: packaged `ARG.csv`) |
| `--total-reads` | no | FASTA record count; if set, appends sample `risk_score` |

### `nanoamrrank run`

| Argument | Required | Description |
| --- | :---: | --- |
| `-f, --fasta` | yes | Input nanopore FASTA (`.fa`, `.fasta`, `.fas`, `.fna`, `.ffn`) |
| `-t, --threads` | yes | Threads for ARGpore2 |
| `-o, --outdir` | yes | Output directory |
| `--argpore-home` | yes* | ARGpore2 install directory (`argpore.sh` lives here). Also reads `ARGPORE_HOME` |
| `--argpore-env` | no | Conda env for ARGpore2 tools (default: `argpore2`, or `ARGPORE_ENV`) |
| `--plascad-env` | no | Conda env that provides `Plascad` (default: `Plascad`) |

\* Required unless `ARGPORE_HOME` is set.

### `nanoamrrank check`

| Argument | Required | Description |
| --- | :---: | --- |
| `--argpore-home` | no | Directory to inspect for ARGpore2 + databases |
| `--argpore-env` | no | Conda env to test for the `argpore2` command, perl, and parallel (default: `argpore2`) |
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

`*_mother_table.scored.tsv` — ARG annotations, WHO priority, mobility factors, and read-level relative `risk_score`.

`*_risk_score_list.tsv` — ranked reads; the sample-level relative AMR risk index (`risk_score`) is appended when `--total-reads` is set (always for `run`).

Packaged reference data (do not duplicate at repo root):

- `src/nanoamrrank/data/rules/` — WHO matching rules
- `src/nanoamrrank/data/reference/ARG.csv` — ARG rank table

## License

A license has not yet been specified. Consult the repository maintainers regarding permission to use, modify, or redistribute the software.

## Citation

See `CITATION.cff`. Related papers:

Y. CHE, et al., Conjugative plasmids interact with insertion sequences to shape the horizontal transfer of antimicrobial resistance genes. PNAS. 118(6):e2008731118. doi: 10.1073/pnas.2008731118

Ziqi Wu, You Che, Chenyuan Dang, Miao Zhang, Xuyang Zhang, Yuhong Sun, Xiang Li, Tong Zhang, Yu Xia*. 2022. Nanopore-based Long-read Metagenomics Uncover the Resistome Intrusion by Antibiotic Resistant Bacteria from Treated Wastewater in Receiving Water Body. Water Research. 226:119282. (https://doi.org/10.1016/j.watres.2022.119282)

Yuxi Yan, Jiao Yang, Miao Zhang, Bixi Zhao, Cailong Nie, Baobo Liu, Chen Wang, Yongsheng Chen, Yu Xia*. A Species-resolved Quantitative Microbial Risk (QMRA) Framework for Hospital Wastewater Based on Nanopore Long-Read Sequencing.
