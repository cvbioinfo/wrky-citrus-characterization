# Caracterização de proteínas WRKY em Citrus

Este fluxo recebe um FASTA de proteínas candidatas WRKY e gera uma tabela final com:

- tamanho da proteína;
- peso molecular;
- ponto isoelétrico (pI);
- peptídeo sinal — SignalP 6;
- localização subcelular — WoLF PSORT;
- domínios proteicos — InterPro;
- número e posição dos domínios WRKY;
- outros domínios encontrados.

Ao final, todas as informações são reunidas em um único arquivo TSV.

---

# 1. Instalação

## Criar o ambiente

```bash
micromamba create -n wrky-analysis python=3.10 pip
micromamba activate wrky-analysis
```

Instalar Biopython:

```bash
pip install biopython
```

O `caracterize_proteins.py` usa Biopython para calcular tamanho, peso molecular e pI.

Os scripts `parse_interpro.py` e `build_final_table.py` usam apenas módulos padrão do Python.

## Instalar SignalP 6

Baixe a versão portátil **SignalP 6 — fast** no site oficial da DTU e extraia o arquivo.

Exemplo de organização:

```text
wrky_analysis/
└── software/
    └── signalp6/
        └── signalp6_fast/
            └── signalp-6-package/
                ├── models/
                ├── signalp/
                ├── setup.py
                └── pyproject.toml
```

Ative o ambiente:

```bash
micromamba activate wrky-analysis
```

Instale uma versão compatível do NumPy:

```bash
pip install "numpy<2"
```

Entre na pasta do pacote:

```bash
cd ~/wrky_analysis/software/signalp6/signalp6_fast/signalp-6-package
```

Instale:

```bash
pip install .
```

Confira:

```bash
signalp6 --version
```

Depois volte à raiz do projeto:

```bash
cd ~/wrky_analysis
```

---

# 2. Preparar o FASTA da espécie

Cada pessoa deve usar um FASTA contendo as proteínas candidatas WRKY da sua espécie.

Exemplo:

```text
data/Csinensis_proteins.fasta
```

Neste README, `Csinensis` é apenas um exemplo. Troque esse nome pelo identificador da espécie analisada em todos os comandos.

Crie uma pasta de resultados:

```bash
mkdir -p results/Csinensis
```

---

# 3. Rodar SignalP

Para proteínas vegetais, usar `eukarya`.

```bash
signalp6 \
  --fastafile ~/wrky_analysis/data/Csinensis_proteins.fasta \
  --output_dir ~/wrky_analysis/results/Csinensis/signalp \
  --organism eukarya \
  --mode fast \
  --format none \
  --bsize 10 \
  --write_procs 2 \
  --torch_num_threads 4 \
  --model_dir ~/wrky_analysis/software/signalp6/signalp6_fast/signalp-6-package/models
```

O arquivo necessário para a próxima etapa será:

```text
results/Csinensis/signalp/prediction_results.txt
```

---

# 4. Calcular tamanho, peso molecular e pI

Execute:

```bash
python scripts/caracterize_proteins.py \
  -i data/Csinensis_proteins.fasta \
  --signalp results/Csinensis/signalp/prediction_results.txt \
  -o results/Csinensis/Csinensis_results.tsv
```

Esse script calcula:

- comprimento;
- peso molecular;
- pI;

e adiciona os resultados do SignalP.

Saída:

```text
results/Csinensis/Csinensis_results.tsv
```

---

# 5. Rodar InterPro

Acesse:

```text
https://www.ebi.ac.uk/interpro/search/sequence/
```

No site:

1. envie `data/Csinensis_proteins.fasta`;
2. submeta todas as sequências;
3. aguarde a análise;
4. clique em **Download**;
5. baixe `TSV output`;
6. baixe `JSON output`.

Salve como:

```text
results/Csinensis/Csinensis_interpro.tsv
results/Csinensis/Csinensis_interpro.json
```

O JSON é usado pelo próximo script.

---

# 6. Processar os domínios do InterPro

Execute:

```bash
python scripts/parse_interpro.py \
  -i results/Csinensis/Csinensis_interpro.json \
  -o results/Csinensis/interpro_parsed
```

O arquivo principal para a etapa final será:

```text
results/Csinensis/interpro_parsed/RE_domains.tsv
```

O parser consolida diferentes bancos que detectam a mesma região proteica. Assim, por exemplo, Pfam, SMART e PROSITE reconhecendo a mesma região WRKY não são contados como três domínios diferentes.

---

# 7. Rodar WoLF PSORT

Acesse:

```text
https://wolfpsort.hgc.jp/
```

No site:

1. selecione **Plant**;
2. cole todas as sequências do FASTA no campo de texto;
3. execute;
4. copie todo o resultado.

Para mostrar o FASTA no terminal:

```bash
cat data/Csinensis_proteins.fasta
```

Salve a saída do WoLF PSORT em:

```text
results/Csinensis/wolfpsort.txt
```

O arquivo deve manter o formato original, por exemplo:

```text
CITRE_001t000750 details nucl: 14
CITRE_001t001510 details nucl: 10, chlo: 1, cyto: 1, mito: 1, extr: 1
```

Não edite manualmente os resultados.

---

# 8. Gerar a tabela final

Quando existirem:

```text
results/Csinensis/Csinensis_results.tsv
results/Csinensis/wolfpsort.txt
results/Csinensis/interpro_parsed/RE_domains.tsv
```

execute:

```bash
python scripts/build_final_table.py \
  --base results/Csinensis/Csinensis_results.tsv \
  --wolfpsort results/Csinensis/wolfpsort.txt \
  --domains results/Csinensis/interpro_parsed/RE_domains.tsv \
  -o results/Csinensis/Csinensis_final_characterization.tsv
```

Resultado:

```text
results/Csinensis/Csinensis_final_characterization.tsv
```

---

# 9. Conferir o resultado

Conte as proteínas do FASTA:

```bash
grep -c '^>' data/Csinensis_proteins.fasta
```

Conte as linhas da tabela final:

```bash
wc -l results/Csinensis/Csinensis_final_characterization.tsv
```

A tabela final deve ter uma linha a mais que o número de proteínas, por causa do cabeçalho.

Exemplo:

```text
56 proteínas no FASTA
57 linhas na tabela final
```

---

# Ordem resumida

```text
FASTA
 |
 +-- SignalP
 |    |
 |    +-- prediction_results.txt
 |           |
 |           +-- caracterize_proteins.py
 |                  |
 |                  +-- Csinensis_results.tsv
 |
 +-- InterPro
 |    |
 |    +-- Csinensis_interpro.json
 |           |
 |           +-- parse_interpro.py
 |                  |
 |                  +-- RE_domains.tsv
 |
 +-- WoLF PSORT
      |
      +-- wolfpsort.txt

Csinensis_results.tsv
+ RE_domains.tsv
+ wolfpsort.txt
        |
        +-- build_final_table.py
                |
                +-- Csinensis_final_characterization.tsv
```

---

# Arquivos principais a manter

```text
data/Csinensis_proteins.fasta

results/Csinensis/
├── Csinensis_results.tsv
├── Csinensis_interpro.tsv
├── Csinensis_interpro.json
├── wolfpsort.txt
├── signalp/
│   └── prediction_results.txt
├── interpro_parsed/
│   └── RE_domains.tsv
└── Csinensis_final_characterization.tsv
```

O principal arquivo para comparação entre espécies é:

```text
Csinensis_final_characterization.tsv
```

---

# Observação

Este fluxo não procura proteínas WRKY em um proteoma completo.

Ele recebe um FASTA de proteínas candidatas previamente selecionadas e caracteriza essas sequências.

Proteínas sem domínio WRKY detectado pelo InterPro não devem ser excluídas automaticamente. Esses casos devem ser revisados separadamente.
