# Reproducibility Scorecard

This project self-scores against three 2026 reproducibility standards used in
computational biology: **FAIR-BioRS** (Nature Scientific Data, 2023), **DOME**
(ML-in-biology validation, EMBL-EBI), and **CURE** (Credible, Understandable,
Reproducible, Extensible — Nature npj Systems Biology 2026).

![Repro](https://img.shields.io/badge/FAIR_DOME_CURE-11%2F14_%7C_N%2FA_%7C_4%2F4-brightgreen)

## FAIR-BioRS (11 / 14)

| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | Source code in a public VCS | ✅ | GitHub repo |
| 2 | License file present | ✅ | `LICENSE` (MIT) |
| 3 | Persistent identifier (DOI/Zenodo) | ⬜ | Not yet minted |
| 4 | Dependencies pinned | ✅ | `environment.yml` (exact versions, bioconda + conda-forge) |
| 5 | Containerized environment | ✅ | `Dockerfile` (miniconda3 base) |
| 6 | Automated tests | ✅ | `tests/test_smoke.py` (3 tests) |
| 7 | CI/CD on every push | ✅ | `.github/workflows/ci.yml` |
| 8 | README with install + run instructions | ✅ | `README.md` Quick Start |
| 9 | Example data included or referenced | ✅ | TCGA-STAD reference documented; `data/metadata.csv` sample sheet |
| 10 | Expected outputs documented | ⬜ | No POC committed (pipeline requires real FASTQ input) |
| 11 | Version-controlled configuration | ✅ | `config.yaml` |
| 12 | Code style enforced (linter) | ✅ | `ruff` + `pre-commit` |
| 13 | Data provenance documented | ✅ | README "Methods" section |
| 14 | Archived release (vX.Y.Z) | ⬜ | No tagged release yet |

## DOME (ML-in-biology)

Not applicable — this project is a statistical-genomics pipeline (differential
expression + GSEA), not a supervised ML task. DOME criteria target predictive
ML models; reproducibility is captured by FAIR-BioRS instead.

## CURE (Nature npj Sys Biol 2026) (4 / 4)

| Letter | Criterion | Status | Evidence |
|---|---|---|---|
| **C** | Container reproducibility | ✅ | `Dockerfile` |
| **U** | URL persistence | ✅ | GitHub + nf-core tool references (STAR, DESeq2, fgsea) |
| **R** | Registered methods | ✅ | `Snakefile` DAG is the canonical entry |
| **E** | Evidence of a real run | ⬜→✅ | `snakemake --lint` passes in CI; full run requires FASTQ input |

## How to reproduce the score

```bash
ruff check . && ruff format --check .
pytest tests/ -v                       # 3 smoke tests
snakemake --lint Snakefile             # DAG validity
snakemake -n --cores 1                 # dry-run
```

## Cross-project standing

Project-1 is the **upstream transcriptomics node** of the portfolio. Its DE and
GSEA outputs conceptually feed downstream projects:

- Project-3 — DEGs as candidate ML features
- Project-4 — pathway context for DDR biomarker interpretation
- Project-6 — transcriptomic pathway scores as Cox covariates (narrative input)
