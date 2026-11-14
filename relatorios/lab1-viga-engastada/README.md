# Laboratório 1 — Viga engastada (PCMEC0147 Vibrações Mecânicas)

Relatório reconstruído a partir do projeto do Overleaf (introdução teórica) e dos dados
dos ensaios. PDF compilado: `Relatorio_Lab1_Viga_Engastada.pdf`.

## Estrutura

```
main.tex                 documento principal
content/                 01 introdução, 02 metodologia, 03 resultados, 04 conclusão, apêndice
settings/                pacotes, comandos e estilos (igual ao projeto original + \trecho e \linhas)
codigo/                  scripts simbólicos da introdução (SymPy) e análise numérica
  parametros.py          dados da viga e conversão lambda_r <-> f_r
  dados.py               leitura dos arquivos e FRF do martelo a partir da série temporal
  identificacao.py       picos e método da meia potência
  modelo.py              equações características, modos normalizados e FRF analítica
  ajuste_molas.py        ajuste de kappa_R e kappa_L aos picos 1 e 2, validação, sensibilidade
  gera_resultados.py     gera figures/*.pdf e tabelas/*.tex a partir de dados/
dados/                   arquivos dos ensaios (renomeados; nomes originais em codigo/dados.py)
figures/                 figuras (SciencePlots, estilo "science" + "grid")
tabelas/                 linhas das tabelas, valores citados no texto e resumo.json
errata_notas.md          correções feitas nas notas manuscritas
```

## Reproduzir

```bash
cd codigo
pip install -r requirements.txt     # precisa de LaTeX instalado (SciencePlots usa usetex)
python gera_resultados.py           # ~1 min; reescreve figures/ e tabelas/
cd .. && latexmk -pdf main.tex
```

O preâmbulo usa `tblr-extras` (TeX Live 2024+, como no Overleaf).
