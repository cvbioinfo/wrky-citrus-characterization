# WRKY characterization workflow for Citrus

This workflow accepts a FASTA of candidate WRKY proteins and produces one final TSV per species. The species/sample prefix is derived from the FASTA filename, so the scripts do not need to be edited for every Citrus species.

## 1. Filename convention

Name the input FASTA `<prefix>_proteins.fasta`, for example:

- `data/RE_proteins.fasta` → prefix `RE`
- `data/Csinensis_proteins.fasta` → prefix `Csinensis`
- `data/Citrus_sinensis_proteins.fasta` → prefix `Citrus_sinensis`

The part before `_proteins` becomes the prefix. If the filename does not end in a recognized suffix, the whole filename stem is used. You can override the detected prefix with `--prefix`.

Keep protein IDs in the FASTA unique and unchanged when submitting the sequences to InterPro and WoLF PSORT; the scripts use those IDs to join results from the different tools.

## 2. Requirements

Activate the project environment:

```bash
micromamba activate wrky-analysis
```

Python 3.10, Biopython and SignalP 6 are required for the local steps. Create the environment and install the Python packages:

```bash
micromamba create -n wrky-analysis python=3.10 pip
micromamba activate wrky-analysis
pip install -r requirements.txt
```

Download the **portable SignalP 6 fast package** from the official DTU website. The download is large and requires following the site's access/licence instructions. Extract it into the project. For example, if the downloaded archive is `~/Downloads/signalp-6.0i.fast.tar.gz`:

```bash
mkdir -p ~/wrky_analysis/software/signalp6
tar -xzf ~/Downloads/signalp-6.0i.fast.tar.gz \
  -C ~/wrky_analysis/software/signalp6
cd ~/wrky_analysis/software/signalp6/signalp6_fast/signalp-6-package
pip install "numpy<2"
pip install .
cd ~/wrky_analysis
signalp6 --version
```

If the downloaded archive has a different name, replace it in the `tar` command. SignalP and its model files are large and are intentionally **not included in this Git repository**. The workflow expects the models at:

```text
software/signalp6/signalp6_fast/signalp-6-package/models
```

If they are elsewhere, provide `--model-dir` when running the pipeline.

## 3. Run the local steps

From the repository root, place your FASTA in `data/` and run, for example:

```bash
python scripts/run_pipeline.py --fasta data/Csinensis_proteins.fasta
```

The workflow derives the prefix (`Csinensis`) and creates `results/Csinensis/`. It runs SignalP and calculates protein length, molecular weight and pI. On this first run, it will tell you the exact paths where the web results need to be saved.

## 4. Run the web analyses

The InterPro and WoLF PSORT submissions are manual because this workflow does not automate their websites.

### InterPro

1. Submit the same FASTA at <https://www.ebi.ac.uk/interpro/search/sequence/>.
2. Download **JSON output** and **TSV output**.
3. Save them using the prefix-derived names in `results/Csinensis/`:

```text
results/Csinensis/Csinensis_interpro.json
results/Csinensis/Csinensis_interpro.tsv
```

### WoLF PSORT

1. Submit the same FASTA at <https://wolfpsort.hgc.jp/>.
2. Select **Plant**.
3. Copy the complete text result to:

```text
results/Csinensis/Csinensis_wolfpsort.txt
```

Replace `Csinensis` with the prefix for your input FASTA.

## 5. Finish the analysis

Run the same command again:

```bash
python scripts/run_pipeline.py --fasta data/Csinensis_proteins.fasta
```

The script reuses up-to-date local results, detects the saved InterPro and WoLF PSORT outputs, parses the InterPro JSON, and builds:

```text
results/Csinensis/Csinensis_final_characterization.tsv
```

You can force the local analyses to run again with:

```bash
python scripts/run_pipeline.py --fasta data/Csinensis_proteins.fasta --force
```

If SignalP's model files are stored elsewhere:

```bash
python scripts/run_pipeline.py \
  --fasta data/Csinensis_proteins.fasta \
  --model-dir /path/to/signalp-6-package/models
```

## Outputs

Each species gets its own folder, preventing results from being overwritten by another species:

```text
results/<prefix>/
├── <prefix>_results.tsv
├── <prefix>_interpro.json
├── <prefix>_interpro.tsv
├── <prefix>_wolfpsort.txt
├── signalp/
│   └── prediction_results.txt
├── interpro_parsed/
│   ├── <prefix>_interpro_metadata.txt
│   ├── <prefix>_interpro_matches.tsv
│   ├── <prefix>_domains.tsv
│   └── <prefix>_domains_summary.tsv
└── <prefix>_final_characterization.tsv
```

The final TSV combines physicochemical properties, SignalP prediction, WoLF PSORT localization and InterPro domains. It does not discover WRKY proteins in a complete proteome; it characterizes the candidate proteins provided in the FASTA. Proteins without a detected WRKY domain should be reviewed rather than automatically discarded.
