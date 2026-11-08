# Diarrhoeagenic *Escherichia coli* Attribution Paradox in Africa & Globally: A Two-Part Systematic Review, Meta-Analysis, and Virulome Case-Control Synthesis (1970–2026)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-green.svg)](https://www.python.org/)
[![Target Journal](https://img.shields.io/badge/Target%20Journal-PLOS%20Neglected%20Tropical%20Diseases-orange.svg)](https://journals.plos.org/plosntds/)

This repository hosts the complete, reproducible computational codebase, curated datasets, statistical pooling models, and high-resolution figures for the dual-track systematic review and meta-analysis resolving the diarrhoeagenic *Escherichia coli* (DEC) disease attribution paradox in African paediatric populations and international cohorts.

---

## 🔬 Study Overview & The Attribution Paradox

Diarrhoeagenic *E. coli* represents a leading bacterial cause of paediatric morbidity in low- and middle-income countries. However, clinical attribution has remained intensely controversial:
- **Part A (Continental Macro-Epidemiology):** In single-arm hospital prevalence series across Africa ($k = 28$ primary cohorts, $N = 15,750$), **Enteroaggregative *E. coli* (EAEC)** is the most prevalent pathotype recovered from children with acute diarrhoea (**$28.3\%$**, $95\%$ CI: $23.2\text{--}33.7\%$). Yet in matched case-control surveillance, healthy asymptomatic community controls harbor EAEC at nearly identical baseline rates ($25.1\%$), yielding a null unstratified attribution odds ratio ($\text{OR} = 1.38$, Knapp–Hartung $95\%$ CI: $0.95\text{--}2.00$). Identical carriage anomalies affect atypical EPEC (`aEPEC`, $\text{OR} = 0.89$) and DAEC ($\text{OR} = 1.21$).
- **Part B (Virulome & Host Ontogeny Resolution):** Meta-analysis of $154$ case-control contingency strata across 29 distinct virulence loci in five continents resolves the paradox:
  1. **EAEC:** Canonical master regulon markers (`aggR`, `aatA`) are ubiquitous colonizers, but specific tissue-damaging autotransporters (**`sepA`**, $\text{OR} = 4.68$, $P = 1.15 \times 10^{-9}$; **`pet`**, $\text{OR} = 2.19$) and high-affinity adhesins (**`aafA`**, $\text{OR} = 3.08$; **`CS22`**, $\text{OR} = 3.59$) drive acute secretory disease ($I^2 = 0.0\%$).
  2. **DAEC Age-Ontogeny Switch:** DAEC is completely non-pathogenic in infants $<12$ months ($\text{OR} = 0.91$, $P = 0.509$) due to maternal breastmilk IgA protection and immature brush-border DAF/CEACAM receptor density, but switches abruptly to a potent diarrhoeal pathogen in toddlers $\ge 12\text{--}59$ months (**$\text{OR} = 2.21$**, $P = 1.84 \times 10^{-10}$).
  3. **Atypical EPEC:** Unstratified `eae+ bfpA-` isolates show null attribution, but strains carrying the **OI-122 pathogenicity island** (**`lifA/efa1`**, $\text{OR} = 6.01$; **`nleB`**, $\text{OR} = 4.20$; **`nleE`**, $\text{OR} = 3.80$) are among the most virulent enteropathogens evaluated.

---

## 📁 Repository Structure

```
dec-attribution-meta-analysis/
├── README.md                      # Comprehensive study guide and documentation
├── LICENSE                        # MIT Open Source License
├── requirements.txt               # Exact Python package dependencies
├── environment.yml                # Conda environment definition
├── data/
│   ├── raw/
│   │   ├── part_a_african_cohorts.csv      # Primary African surveillance cohorts (k=28, N=15,750)
│   │   └── part_b_global_virulence_loci.csv# Global & African 2x2 contingency tables (154 strata)
│   └── processed/
│       ├── meta_analysis_summary_statistics.json # Complete pooled estimates, REML/DL & I2 metrics
│       └── meta_regression_models.json           # Multivariate meta-regression model parameters
├── scripts/
│   ├── 01_prevalence_and_odds_ratio_pooling.py   # Primary Freeman-Tukey & odds ratio meta-analysis
│   ├── 02_multivariate_meta_regression.py        # Mixed-effects REML meta-regression
│   ├── 03_diagnostic_presence_absence_matrix.py  # 68-cohort diagnostic tile matrix (Fig 2A / Fig S18)
│   ├── 04_methodological_mindmap.py              # Visual evidence architecture (Fig 1B)
│   ├── 05_generate_supplementary_figures.py      # Complete S1-S22 supplementary figure suite
│   └── 06_format_landscape_tables.py             # OpenXML landscape formatting engine
├── figures/
│   ├── main/                      # Publication-ready main manuscript figures (PNG & PDF)
│   └── supplementary/             # High-resolution supplementary figures (S1 Fig to S22 Fig)
└── manuscript/
    ├── manuscript.md              # Full submission-ready manuscript (PLOS NTD format)
    ├── manuscript.pdf             # Formatted PDF with 18 landscape data tables
    ├── supplementary_information.md # Complete supplementary text, tables, and notes
    └── supplementary_information.pdf# Formatted Supplementary PDF with 25 landscape pages
```

---

## 📊 Data Dictionaries

### `data/raw/part_a_african_cohorts.csv`
Contains study-level metadata, sample sizes, and pathotype recovery counts across $k = 28$ primary cohorts in 18 African nations:
- `study_id`: Unique identifier (e.g., `GEMS_2013`, `VIDA_2023`, `AKINLABI_2023A`).
- `author_year`: APA citation and surveillance year.
- `country`: Study nation and multi-site locations.
- `un_subregion`: United Nations African Subregion (Western, Eastern, Southern, Central, Northern, Pan-African).
- `study_design`: Study classification (`Matched Case-Control` vs. `Cross-Sectional / Genomic Surveillance`).
- `modality`: Primary diagnostic platform deployed (`Conventional Culture`, `Multiplex PCR`, `Quantitative TAC`, `WGS`).
- `cases_n` / `controls_n`: Total diarrhoeal case and asymptomatic control counts.
- `eaec_pos`, `etec_pos`, `epec_pos`, `tepec_pos`, `aepec_pos`, `daec_pos`, `stec_pos`, `eiec_pos`, `hybrid_pos`: Pathotype-positive counts in clinical cases.
- `*_ctrl_pos`: Corresponding pathotype counts in matched healthy community controls.

### `data/raw/part_b_global_virulence_loci.csv`
Contains 154 granular $2 \times 2$ contingency tables across 29 individual virulence genes:
- `study_id`: Cohort identifier.
- `country` / `region`: Geographic origin across Africa, Latin America, Europe, Asia, and Oceania.
- `pathotype`: Evaluated pathotype (`EAEC`, `DAEC`, or `EPEC`).
- `locus`: Virulence determinant (e.g., `sepA`, `aafA`, `CS22`, `pet`, `sat`, `lifA_efa1`, `nleB`, `nleE`, `daaC`, `bfpA`).
- `cases_pos`, `cases_total`: Positive and total counts in diarrhoeal cases.
- `ctrl_pos`, `ctrl_total`: Positive and total counts in matched healthy controls.
- `age_strata`: Host age group (`<12m`, `12-59m`, or `All`).

---

## 🚀 Reproduction & Execution

### 1. Environment Setup
```bash
# Clone the repository
git clone https://github.com/TochiAgboeze/dec-attribution-meta-analysis.git
cd dec-attribution-meta-analysis

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install requirements
pip install -r requirements.txt
```

### 2. Execute Statistical Models & Generate Figures
```bash
# 1. Run complete prevalence and case-control odds ratio pooling
python3 scripts/01_prevalence_and_odds_ratio_pooling.py

# 2. Run multivariate mixed-effects meta-regression
python3 scripts/02_multivariate_meta_regression.py

# 3. Generate 68-cohort diagnostic presence/absence matrix (Fig 2A / Fig S18)
python3 scripts/03_diagnostic_presence_absence_matrix.py

# 4. Generate the full supplementary figure suite (S1–S22 Fig)
python3 scripts/05_generate_supplementary_figures.py

# 5. Format Word tables into landscape publication format
python3 scripts/06_format_landscape_tables.py
```

---

## 📜 Citation & Attribution

If you utilize the datasets, statistical pooling models, or diagnostic methodologies from this repository, please cite:

```bibtex
@article{Agboeze2026DECAttribution,
  title={Diarrhoeagenic Escherichia coli Attribution Paradox in Africa and Globally: A Two-Part Systematic Review, Freeman–Tukey Random-Effects Meta-Analysis, Diagnostic Sensitivity Evaluation, and Virulome Case-Control Synthesis (1970–2026)},
  author={Agboeze, Tochukwu and Okeke, Iruka N. and {African DEC Genomics Surveillance Consortium}},
  journal={PLOS Neglected Tropical Diseases},
  year={2026},
  note={Under Review}
}
```

---

## 📄 License
This repository is licensed under the [MIT License](LICENSE). Curated datasets and figures are licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).
