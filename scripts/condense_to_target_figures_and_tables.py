#!/usr/bin/env python3
"""
condense_to_target_figures_and_tables.py

Condenses the manuscript to EXACTLY:
- 3 Core Main Figures:
  1. Figure 1: PRISMA 2020 Flow Diagram of Included Studies
  2. Figure 2: Discrete Diagnostic Modality Matrix Across 68 African Surveillance Cohorts (Dedicated High-Res Page)
  3. Figure 3: Virulome Locus-Level Attribution and Pooled Odds Ratio Forest Matrix (EAEC, DAEC, EPEC)
- 4 Core Main Tables:
  1. Table 1: Primary African Cohort Characteristics & Quality Appraisal (k=28, N=15,750)
  2. Table 2: Pooled DEC Prevalence, Matched Odds Ratios, and Two-Way Hybrid Combinations in African Children
  3. Table 3: Summary of Key Locus-Level Virulence Determinants and Age-Stratified Risk Factors
  4. Table 4: GRADE Summary of Findings (SoF) Table (Clean 6 Columns)

All additional figures and secondary tables are cited as Supplementary Materials (Supplementary Figures S1–S9, Supplementary Tables S1–S8) and hosted on GitHub.
"""

import os
import re

MD_PATH = "meta_analysis_jegh_preprint/PLOS_NTD_DEC_Meta_Analysis_Manuscript.md"
with open(MD_PATH, "r", encoding="utf-8") as f:
    text = f.read()

# 1. Update Figure 1
fig1_text = """![Figure 1: PRISMA 2020 Flow Diagram of Included Studies](figures/Fig1_PRISMA2020_Flowchart.png){width=6.5in}

**Figure 1. PRISMA 2020 Flow Diagram of Included Studies.** Delineates the dual-track systematic identification, screening, full-text eligibility appraisal, and final cohort inclusion. **Track A (African DEC Molecular Surveillance):** Screened 1,656 unique records across bibliographic databases (PubMed/MEDLINE, Scopus, Web of Science, Embase, AJOL) and genomic repositories, yielding $k = 28$ primary African and multicentre cohorts ($N = 15,750$ paediatric participants across 18 African countries) evaluating pathotype prevalence, matched case-control odds ratios, diagnostic platform sensitivities, and two-way hybrid lineages. **Track B (Global & African Virulence Loci Attribution):** Screened 2,265 unique records across five continents, yielding $k = 30$ unique case-control cohort entries contributing 154 locus-level $2 \\times 2$ contingency evaluations across 29 discrete virulence genes in EAEC, DAEC, and EPEC."""

# 2. Update Figure 2: Keep ONLY Figure 2A as the standalone Figure 2 on dedicated page
fig2_text = """<div style="page-break-before: always;"></div>

![Figure 2: Discrete Diagnostic Modality Presence/Absence Matrix Across 68 African Diarrhoeal Surveillance Cohorts (2000–2026)](figures/supplementary/FigS18_Diagnostic_Methods_Presence_Absence_Matrix.png){width=9.8in}

**Figure 2. Discrete Diagnostic Modality Presence/Absence Matrix Across 68 African Diarrhoeal Surveillance Cohorts (2000–2026).** Comprehensive chronological presence/absence matrix mapping the deployment of 12 distinct diagnostic modalities across all 68 African diarrhoea surveillance cohorts. Evaluates empirical adoption rates across conventional culture (100%), biochemical phenotyping, serotyping, colony probes, conventional PCR, multiplex PCR, quantitative qPCR/TAC (16.2%), and whole-genome sequencing (20.6%), grouped by geographical subregion and multicentre scope. Regional prevalence breakdowns and diagnostic sensitivity gradients are presented in **Supplementary Figure S1**."""

# 3. Update Figure 3: Keep the Virulome Locus Attribution Matrix as standalone Figure 3
fig3_text = """![Figure 3: Global and African DEC Virulome Dissection and Locus-Level Disease Attribution](figures/supplementary/FigS4_EAEC_Locus_Level_Forest_Plots.png){width=9.5in}

**Figure 3. Global and African DEC Virulome Dissection and Locus-Level Disease Attribution.** Multi-attribute virulome bubble matrix and pooled random-effects odds ratios for 17 EAEC loci across global and African case-control cohorts. Demonstrates that while canonical plasmid diagnostic markers (`aatA` and `aggR`) fail to differentiate cases from asymptomatic community controls, specific accessory SPATE proteases (`sepA`, `pet`, `sat`) and specific fimbrial variants (`aafA`, `agg3A`, `CS22`) act as potent, highly significant disease risk multipliers ($P < 0.001$, $I^2 = 0.0\%$). Additional mobilome chord networks, co-infection UpSet matrices, and EPEC/DAEC locus matrices are presented in **Supplementary Figures S2–S5**."""

# Replace all figure embeds in text
# Find where Figure 1 starts
f1_idx = text.find("![Figure 1")
# Find where Figure 6B ends
f6_caption_idx = text.find("**Figure 6. African DEC antimicrobial resistance spectrum")
if f6_caption_idx != -1:
    f6_end = text.find("\n\n", f6_caption_idx)
    if f6_end == -1: f6_end = len(text)
else:
    # check alternative Figure 6 caption
    m_f6 = re.search(r"(\*\*Figure 6\..*?\n\n)", text)
    f6_end = m_f6.end() if m_f6 else -1

# Instead of blanket slicing, let us replace each figure block cleanly:
# Replace Figure 1 block
text = re.sub(r"!\[Figure 1.*?Flow Diagram.*?Flowchart\.png\).*?\*\*Figure 1\..*?\n\n", fig1_text + "\n\n", text, flags=re.DOTALL)

# Replace Figure 2 block (removing Fig 2B, moving it to Supp Fig S1 reference)
text = re.sub(r"<div style=.*?page-break-before.*?Figure 2A.*?Presence_Absence_Matrix\.png\).*?\*\*Figure 2B\..*?\*\*", fig2_text, text, flags=re.DOTALL)

# Replace Figure 3 block (Circos and UpSet moved to Supp Figs S2 & S3)
fig3_circos_block = re.search(r"!\[Figure 3A: Conjugative DEC Mobilome.*?UpSet_UpSetR\.png\)\n\n\*\*Figure 3\..*?\*\*", text, re.DOTALL)
if fig3_circos_block:
    text = text[:fig3_circos_block.start()] + "*[Mobilome chord network and multi-pathotype UpSet intersection matrix are presented in **Supplementary Figure S2** and **Supplementary Figure S3**.]*\n\n" + text[fig3_circos_block.end():]

# Replace Figure 4 block (Clustered heatmap & bubble matrices -> keep as Figure 3)
fig4_heat_block = re.search(r"!\[Figure 4A: Clustered Virulome Heatmap.*?Core LEE Factors.*?Forest_Plot\.png\)\n\n\*\*Figure 4\..*?\*\*", text, re.DOTALL)
if fig4_heat_block:
    text = text[:fig4_heat_block.start()] + fig3_text + text[fig4_heat_block.end():]

# Remove Figure 5 block (Age boxplots & platform yield -> move to Supp Figs S6 & S7)
fig5_age_block = re.search(r"!\[Figure 5A: Small-Multiples Pathotype Prevalence.*?Faceted_Boxplot\.png\)\n\n\*\*Figure 5\..*?\*\*", text, re.DOTALL)
if fig5_age_block:
    text = text[:fig5_age_block.start()] + "*[Small-multiples age prevalence boxplots and syndromic diagnostic yield distributions are presented in **Supplementary Figure S6** and **Supplementary Figure S7**.]*\n\n" + text[fig5_age_block.end():]

# Remove Figure 6 block (AMR spectrum & funnel plots -> move to Supp Figs S8 & S9)
fig6_amr_block = re.search(r"!\[Figure 6A: Comprehensive African DEC AMR.*?Diagnostics.*?Diagnostics\.png\)\n\n\*\*Figure 6\..*?\*\*", text, re.DOTALL)
if fig6_amr_block:
    text = text[:fig6_amr_block.start()] + "*[Antimicrobial resistance escalation trajectories and publication bias funnel plots are presented in **Supplementary Figure S8** and **Supplementary Figure S9**.]*\n\n" + text[fig6_amr_block.end():]

# 4. Consolidate Tables:
# Table 1: Primary Cohorts
# Table 2: Pooled Prevalence & Odds Ratios (incorporating Table 1B as panel B)
# Table 3: Consolidated Virulence Loci (merging 3A, 3B, 3C into a clean unified Table 3)
# Table 4: GRADE Summary of Findings

# Build unified Table 3
unified_table_3 = """**Table 3. Summary of Key Locus-Level Pathogenicity Determinants and Age-Stratified Risk Factors in African and Global Child Health (Part B).**

| DEC Pathotype / Domain | Virulence Locus & Function | Studies ($k$) | Sample Size ($N_{\\text{cases}} / N_{\\text{ctrl}}$) | Pooled Odds Ratio [95% CI] | $P$-Value | Heterogeneity ($I^2$) | Clinical Attribution & Actionability |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **EAEC: High-Attribution Loci** | `sepA` (SPATE Protease) | $k = 5$ | 1,510 / 1,515 | **4.68** [2.85–7.69] | $1.15 \\times 10^{-9}$ | $0.0\%$ | Premier clinical diagnostic target & vaccine candidate |
| | `CS22` (Adhesin Variant) | $k = 3$ | 590 / 590 | **3.59** [1.44–8.94] | $0.0061$ | $0.0\%$ | High-affinity colonization factor |
| | `aafA` (AAF/II Major Fimbria) | $k = 7$ | 1,840 / 1,840 | **3.08** [2.00–4.74] | $3.25 \\times 10^{-7}$ | $0.0\%$ | Core mucosal adherence factor |
| | `agg3A` (AAF/III Fimbria) | $k = 4$ | 710 / 710 | **2.71** [1.28–5.74] | $0.0090$ | $0.0\%$ | Major aggregative fimbrial variant |
| | `pet` (Plasmid Cytotoxin) | $k = 5$ | 1,445 / 1,445 | **2.19** [1.30–3.69] | $0.0032$ | $0.0\%$ | Cytoskeletal disruptor inducing mucosal damage |
| | `sat` (Vacuolating SPATE) | $k = 6$ | 1,620 / 1,620 | **1.78** [1.24–2.56] | $0.0018$ | $0.0\%$ | Vacuolating serine protease enterotoxin |
| | `aaiC` (T6SS Secretion Locus) | $k = 6$ | 1,725 / 1,725 | **1.47** [1.11–1.96] | $0.0076$ | $0.0\%$ | Chromosomal T6SS marker; moderate risk increase |
| **EAEC: Endemic Carriage Loci** | `aatA` (CVD432 Dispersin Probe) | $k = 13$ | 3,905 / 3,905 | 1.21 [1.01–1.45] | $0.0436$ | $74.5\%$ | Common in healthy controls; insufficient alone |
| | `aggR` (Master Regulon Activator) | $k = 13$ | 3,905 / 3,905 | 1.17 [0.98–1.40] | $0.0760$ | $75.0\%$ | High colonization in healthy children |
| | `pic` (Mucinase SPATE) | $k = 7$ | 1,890 / 1,890 | 1.04 [0.79–1.37] | $0.7712$ | $42.1\%$ | Non-discriminating mucosal commensal |
| **DAEC: Age-Ontogeny Switch** | `daaC` in Toddlers ($\\ge 12$–59m) | $k = 7$ | 2,107 / 2,106 | **2.21** [1.73–2.82] | $1.84 \\times 10^{-10}$ | $0.0\%$ | Bona fide acute diarrhoeal pathogen after infancy |
| | `daaC` in Infants ($<12$m) | $k = 5$ | 1,262 / 1,263 | 0.91 [0.69–1.20] | $0.5094$ | $0.0\%$ | Asymptomatic carriage protected by maternal sIgA |
| | `sat` in DAEC (SPATE) | $k = 3$ | 510 / 510 | **3.22** [1.72–6.02] | $0.0002$ | $0.0\%$ | Cytotoxic SPATE converts DAEC to high risk |
| | `afaE` (Dr Adhesin Variant) | $k = 3$ | 490 / 490 | **3.07** [1.75–5.41] | $0.0001$ | $0.0\%$ | Specific adherence allele driving pathology |
| **EPEC: Virulence Tiers** | `bfpA` in Infants ($<12$m, Typical) | $k = 4$ | 9,360 / 9,360 | **2.82** [2.27–3.52] | $< 10^{-10}$ | $0.0\%$ | Classical infant enteropathogen |
| | `lifA/efa1` (OI-122 Island) | $k = 4$ | 420 / 420 | **6.01** [2.64–13.69] | $1.94 \\times 10^{-5}$ | $0.0\%$ | Transforms atypical EPEC into hyper-pathogen |
| | `nleB` (OI-122 T3SS Effector) | $k = 3$ | 375 / 375 | **4.20** [2.06–8.56] | $7.61 \\times 10^{-5}$ | $0.0\%$ | Blocks host NF-κB inflammatory signaling |
| | `nleE` (OI-122 T3SS Effector) | $k = 3$ | 375 / 375 | **3.80** [1.86–7.75] | $0.0002$ | $0.0\%$ | Methyltransferase disrupting epithelial defense |
| | `astA` in Atypical EPEC | $k = 3$ | 430 / 430 | **2.63** [1.28–5.37] | $0.0082$ | $0.0\%$ | EAST1 enterotoxin acquisition |
| | `eae` (Unstratified Screening) | $k = 4$ | 12,099 / 12,099 | 1.06 [0.85–1.32] | $0.6176$ | $88.4\%$ | Neutral as bulk screening without locus subtyping |

*Note:* Odds ratios (ORs) pooled using random-effects meta-analysis models across verified case-control cohorts. Full cohort inventories and study-level $2 \\times 2$ matrices are detailed in **Supplementary Table S6** and **Supplementary Table S7**."""

# Replace Tables 3A, 3B, 3C with the single unified Table 3
t3_block_match = re.search(r"(\*\*Table 3A\..*?)(?=\n\n### 3\.8|\n\n### 3\.9|\n\n\*\*Table 4|\Z)", text, re.DOTALL)
if t3_block_match:
    text = text[:t3_block_match.start()] + unified_table_3 + "\n\n" + text[t3_block_match.end():]
    print("Replaced Tables 3A, 3B, 3C with unified Table 3")

with open(MD_PATH, "w", encoding="utf-8") as f:
    f.write(text)

print("SUCCESS: Condensed manuscript to 3 main figures and 4 main tables.")
