# WRKY protein characterization

Fluxo para caracterizar as proteínas do arquivo `data/RE_proteins.fasta`.

## Ordem de execução

### 1. SignalP

Rodar o SignalP 6:

```bash
signalp6 \
  --fastafile ~/wrky_analysis/data/RE_proteins.fasta \
  --output_dir ~/wrky_analysis/results/signalp \
  --organism eukarya \
  --mode fast \
  --format none \
  --bsize 10 \
  --write_procs 2 \
  --torch_num_threads 4 \
  --model_dir ~/wrky_analysis/software/signalp6/signalp6_fast/signalp-6-package/models
```

Isso gera:

```text
results/signalp/prediction_results.txt
```

### 2. Propriedades físico-químicas + SignalP

```bash
python scripts/caracterize_proteins.py \
  -i data/RE_proteins.fasta \
  --signalp results/signalp/prediction_results.txt \
  -o results/RE_results.tsv
```

Calcula tamanho, peso molecular e pI e adiciona os resultados do SignalP.

### 3. InterPro

Enviar `data/RE_proteins.fasta` ao site do InterPro.

Baixar os resultados em:

```text
results/RE_interpro.tsv
results/RE_interpro.json
```

O JSON é usado no passo seguinte.

### 4. Processar o InterPro

```bash
python scripts/parse_interpro.py \
  -i results/RE_interpro.json \
  -o results/interpro_parsed
```

O arquivo principal usado depois é:

```text
results/interpro_parsed/RE_domains.tsv
```

### 5. WoLF PSORT

Enviar o mesmo FASTA ao WoLF PSORT usando a opção **Plant**.

Salvar o resultado em:

```text
results/wolfpsort.txt
```

### 6. Gerar a tabela final

```bash
python scripts/build_final_table.py \
  --base results/RE_results.tsv \
  --wolfpsort results/wolfpsort.txt \
  --domains results/interpro_parsed/RE_domains.tsv \
  -o results/RE_final_characterization.tsv
```

Resultado final:

```text
results/RE_final_characterization.tsv
```

## Resumo

```text
RE_proteins.fasta
├── SignalP → prediction_results.txt
│   └── caracterize_proteins.py → RE_results.tsv
├── InterPro → RE_interpro.json
│   └── parse_interpro.py → RE_domains.tsv
└── WoLF PSORT → wolfpsort.txt

RE_results.tsv + RE_domains.tsv + wolfpsort.txt
└── build_final_table.py
    └── RE_final_characterization.tsv
```
