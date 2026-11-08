#!/usr/bin/env python3
"""
generate_figs18_discrete_tile_matrix.py

Generates the Discrete Presence/Absence Diagnostic Methods Matrix for Supplementary Figure S18.
Enhanced with a dedicated UN African Subregion metadata layer:
1. Top Panel: Integer count of distinct diagnostic modalities deployed per study (1 to 6).
2. Main Heatmap Grid: 68 African cohorts (chronological 2000–2026) x 12 Diagnostic Modalities.
   Each present method is filled with its distinctive thematic color and hatch pattern.
   Absent cells are subtle crisp slate (#F8FAFC).
3. Regional Layer Strip: Dedicated horizontal metadata track displaying the UN African Subregion
   for each of the 68 cohorts (Western, Eastern, Southern, Central, Northern, Pan-African).
4. Right Margin Bar: Cumulative number of cohorts deploying each modality (1 to 68).
5. Comprehensive Legend Cards: Both for the 12 diagnostic modalities and the 6 UN African Subregions.

Outputs saved to both figures/supplementary/ and meta_analysis_jegh_preprint/figures/supplementary/.
"""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_PNG1 = os.path.join(BASE_DIR, "figures/supplementary/FigS18_Diagnostic_Methods_Presence_Absence_Matrix.png")
OUT_PNG2 = os.path.join(BASE_DIR, "figures/main/fig2a_diagnostic_methods_matrix.png")
OUT_PDF1 = os.path.join(BASE_DIR, "figures/supplementary/FigS18_Diagnostic_Methods_Presence_Absence_Matrix.pdf")
OUT_PDF2 = os.path.join(BASE_DIR, "figures/main/fig2a_diagnostic_methods_matrix.pdf")

os.makedirs(os.path.dirname(OUT_PNG1), exist_ok=True)
os.makedirs(os.path.dirname(OUT_PNG2), exist_ok=True)

# 12 Verified Diagnostic Modalities
METHODS = [
    "Conventional culture",
    "Biochemical tests",
    "Serotyping",
    "Cell adherence assays",
    "Colony DNA hybridization probes",
    "Immuno-assays",
    "Conventional PCR",
    "Multiplex PCR",
    "Quantitative qPCR / TAC",
    "PFGE",
    "MLST",
    "WGS / NGS"
]

COLORS = [
    "#0284C7",  # Culture: Bright Cerulean Blue
    "#10B981",  # Biochem: Emerald Green
    "#F97316",  # Serotyping: Vibrant Coral Orange
    "#8B5CF6",  # Adherence: Royal Violet / Purple
    "#06B6D4",  # Colony probes: Crisp Turquoise
    "#EC4899",  # Immuno-assays: Deep Fuchsia Pink
    "#64748B",  # Conv PCR: Slate Grey
    "#E11D48",  # Multiplex PCR: Rich Crimson
    "#EAB308",  # qPCR / TAC: Gold / Yellow
    "#4338CA",  # PFGE: Indigo
    "#D97706",  # MLST: Amber
    "#9333EA"   # WGS / NGS: Purple / Violet
]

HATCHES = [
    "",       # Culture
    "///",    # Biochem
    "\\\\\\", # Serotyping
    "xx",     # Adherence
    "..",     # Colony probes
    "|||",    # Immuno-assays
    "---",    # Conv PCR
    "+++",    # Multiplex PCR
    "***",    # qPCR / TAC
    "oo",     # PFGE
    "//..",   # MLST
    "XX"      # WGS / NGS
]

# UN African Subregion Color Palette
REGION_COLORS = {
    "Western Africa": "#1D4ED8",           # Deep Royal Blue
    "Eastern Africa": "#10B981",           # Vibrant Emerald Green
    "Southern Africa": "#F59E0B",          # Warm Amber Gold
    "Central Africa": "#8B5CF6",           # Royal Violet Purple
    "Northern Africa": "#EC4899",          # Vibrant Rose Magenta
    "Pan-African Multi-region": "#06B6D4"  # Bright Cyan
}

# 68 African cohorts ordered chronologically (2000–2026) with Country, ISO code, and Subregion
STUDIES = [
    # Era 1: Classical Phenotyping & Colony Probes (2000–2008, k = 10)
    ('Okeke et al. (2000)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Cell adherence assays', 'Colony DNA hybridization probes', 'Conventional PCR']),
    ('Obi et al. (2004)', 'South Africa', 'ZA', 'Southern Africa',
     ['Conventional culture', 'Biochemical tests', 'Colony DNA hybridization probes', 'Conventional PCR']),
    ('Dow et al. (2006)', 'Botswana', 'BW', 'Southern Africa',
     ['Conventional culture', 'Serotyping', 'Conventional PCR', 'Multiplex PCR']),
    ('Al-Gallas et al. (2007)', 'Tunisia', 'TN', 'Northern Africa',
     ['Conventional culture', 'Biochemical tests', 'Serotyping', 'Cell adherence assays', 'Conventional PCR', 'PFGE']),
    ('Mandomando et al. (2007)', 'Mozambique', 'MZ', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Metwally et al. (2007)', 'Egypt', 'EG', 'Northern Africa',
     ['Conventional culture', 'Conventional PCR', 'Multiplex PCR']),
    ('Moyo et al. (2007)', 'Tanzania', 'TZ', 'Eastern Africa',
     ['Conventional culture', 'Biochemical tests', 'Multiplex PCR']),
    ('Ateba et al. (2008)', 'South Africa', 'ZA', 'Southern Africa',
     ['Conventional culture', 'Biochemical tests', 'Serotyping', 'Conventional PCR']),
    ('Oundo et al. (2008)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Conventional PCR', 'PFGE']),
    ('Raji et al. (2008)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Serotyping', 'Conventional PCR']),

    # Era 2: Molecular Gel Multiplex PCR & PFGE Era (2009–2015, k = 19)
    ('Onyango et al. (2009)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Biochemical tests', 'Conventional PCR', 'Multiplex PCR']),
    ('Oundo et al. (2009)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Conventional PCR', 'PFGE']),
    ('Nweze et al. (2010)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Serotyping', 'Multiplex PCR']),
    ('Opintan et al. (2010)', 'Ghana', 'GH', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Multiplex PCR']),
    ('Rahouma et al. (2011)', 'Libya', 'LY', 'Northern Africa',
     ['Conventional culture', 'Biochemical tests', 'Serotyping', 'Conventional PCR']),
    ('Boisen et al. (2012)', 'Mali', 'ML', 'Western Africa',
     ['Conventional culture', 'Multiplex PCR', 'MLST', 'WGS / NGS']),
    ('Osman et al. (2012)', 'Sudan', 'SD', 'Northern Africa',
     ['Conventional culture', 'Biochemical tests', 'Serotyping', 'Conventional PCR']),
    ('Sang et al. (2012)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Tau et al. (2012)', 'South Africa', 'ZA', 'Southern Africa',
     ['Conventional culture', 'Conventional PCR', 'PFGE']),
    ('Kotloff et al. (2013) [GEMS]', 'Multi-site', 'Pan-Afr', 'Pan-African Multi-region',
     ['Conventional culture', 'Biochemical tests', 'Cell adherence assays', 'Multiplex PCR', 'WGS / NGS']),
    ('Swierczewski et al. (2013)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Ali et al. (2014)', 'Egypt', 'EG', 'Northern Africa',
     ['Conventional culture', 'Biochemical tests', 'Conventional PCR', 'Multiplex PCR']),
    ('Kabayiza et al. (2014)', 'Rwanda', 'RW', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC']),
    ('Lääveri et al. (2014)', 'Benin', 'BJ', 'Western Africa',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC']),
    ('Onanuga et al. (2014)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Conventional PCR']),
    ('Stensvold et al. (2014)', 'Guinea-Bissau', 'GW', 'Western Africa',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC']),
    ('Diarra et al. (2015)', 'Mali', 'ML', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Multiplex PCR']),
    ('Ifebunandu et al. (2015)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Conventional PCR', 'Multiplex PCR']),
    ('Platts-Mills et al. (2015) [MAL-ED]', 'Multi-site', 'Pan-Afr', 'Pan-African Multi-region',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC']),

    # Era 3: Genomic & Multi-Plex Expansion (2016–2022, k = 31)
    ('Alabi et al. (2016)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Garrine et al. (2016)', 'Mozambique', 'MZ', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Hegazy et al. (2016)', 'Egypt', 'EG', 'Northern Africa',
     ['Conventional culture', 'Biochemical tests', 'Multiplex PCR']),
    ('Mbae et al. (2016)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('El-Sharkawy et al. (2017)', 'Egypt', 'EG', 'Northern Africa',
     ['Conventional culture', 'Biochemical tests', 'Multiplex PCR']),
    ('Hebbelstrup Jensen et al. (2017)', 'Guinea-Bissau', 'GW', 'Western Africa',
     ['Conventional culture', 'Multiplex PCR', 'WGS / NGS']),
    ('Ifeanyi et al. (2017)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Conventional PCR', 'Multiplex PCR']),
    ('Loadman et al. (2017)', 'South Africa', 'ZA', 'Southern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Mshana et al. (2017)', 'Tanzania', 'TZ', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Rajendram et al. (2017)', 'Madagascar', 'MG', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Svennerholm et al. (2017)', 'Multi-site', 'Pan-Afr', 'Pan-African Multi-region',
     ['Conventional culture', 'Multiplex PCR', 'Immuno-assays']),
    ('Vilchez et al. (2017)', 'Multi-site', 'Pan-Afr', 'Pan-African Multi-region',
     ['Conventional culture', 'Multiplex PCR']),
    ('Abba et al. (2018)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Conventional PCR']),
    ('Duan et al. (2018)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC']),
    ('Komba et al. (2018)', 'Tanzania', 'TZ', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Krumkamp et al. (2018)', 'Ghana', 'GH', 'Western Africa',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC']),
    ('Manning et al. (2018)', 'Mali', 'ML', 'Western Africa',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC']),
    ('Ochieng et al. (2018)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Dallman et al. (2019)', 'Multi-site', 'Pan-Afr', 'Pan-African Multi-region',
     ['Conventional culture', 'WGS / NGS']),
    ('Kopitopoulos et al. (2019)', 'DRC', 'CD', 'Central Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Levine et al. (2020) [GEMS-1A]', 'Multi-site', 'Pan-Afr', 'Pan-African Multi-region',
     ['Conventional culture', 'Cell adherence assays', 'Multiplex PCR', 'Quantitative qPCR / TAC', 'WGS / NGS']),
    ('Matias et al. (2020)', 'Angola', 'AO', 'Central Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Modgil et al. (2020)', 'Multi-site', 'Pan-Afr', 'Pan-African Multi-region',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC']),
    ('Siddik et al. (2020)', 'Ethiopia', 'ET', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Suleiman et al. (2020)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Biochemical tests', 'Conventional PCR']),
    ('Amézquita-López et al. (2021)', 'Multi-site', 'Pan-Afr', 'Pan-African Multi-region',
     ['Conventional culture', 'Multiplex PCR']),
    ('El-Sayed et al. (2021)', 'Egypt', 'EG', 'Northern Africa',
     ['Conventional culture', 'Biochemical tests', 'Multiplex PCR']),
    ('Juma et al. (2021)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR']),
    ('Kongsbak-Wismann et al. (2021)', 'Guinea-Bissau', 'GW', 'Western Africa',
     ['Conventional culture', 'WGS / NGS']),
    ('Okeke et al. (2021)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Cell adherence assays', 'Multiplex PCR', 'WGS / NGS']),
    ('Tchatchouang et al. (2021)', 'Cameroon', 'CM', 'Central Africa',
     ['Conventional culture', 'Biochemical tests', 'Multiplex PCR']),

    # Era 4: Modern High-Resolution WGS & Hybrid Surveillance (2023–2026, k = 8)
    ('Akinlabi et al. (2023a)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Cell adherence assays', 'Multiplex PCR', 'MLST', 'WGS / NGS']),
    ('Akinlabi et al. (2023b)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Cell adherence assays', 'Multiplex PCR', 'PFGE', 'MLST', 'WGS / NGS']),
    ('Kiiru et al. (2024)', 'Kenya', 'KE', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR', 'MLST', 'WGS / NGS']),
    ('Wolde et al. (2024)', 'Ethiopia', 'ET', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR', 'MLST', 'WGS / NGS']),
    ('Mahamat et al. (2024)', 'Cameroon', 'CM', 'Central Africa',
     ['Conventional culture', 'Multiplex PCR', 'MLST', 'WGS / NGS']),
    ('Medugu et al. (2025)', 'Nigeria', 'NG', 'Western Africa',
     ['Conventional culture', 'Multiplex PCR', 'MLST', 'WGS / NGS']),
    ('Garrine et al. (2025)', 'Mozambique', 'MZ', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR', 'Quantitative qPCR / TAC', 'WGS / NGS']),
    ('Kuchaka et al. (2026)', 'Tanzania', 'TZ', 'Eastern Africa',
     ['Conventional culture', 'Multiplex PCR', 'MLST', 'WGS / NGS'])
]

num_studies = len(STUDIES)
num_methods = len(METHODS)

# Create binary presence/absence matrix (num_methods rows x num_studies cols)
matrix = np.zeros((num_methods, num_studies), dtype=int)
for col_idx, s in enumerate(STUDIES):
    m_list = s[4]
    for m in m_list:
        if m in METHODS:
            row_idx = METHODS.index(m)
            matrix[row_idx, col_idx] = 1

study_counts = np.sum(matrix, axis=0) # Methods per study
method_counts = np.sum(matrix, axis=1) # Studies per method

# Regional counts
subregion_list = [s[3] for s in STUDIES]
from collections import Counter
reg_counter = Counter(subregion_list)

# Set up figure with GridSpec:
# Row 0: Top count bar (Methods / Study)
# Row 1: Main tile matrix (12 Diagnostic Modalities)
# Row 2: Regional Layer Strip (UN African Subregion)
fig = plt.figure(figsize=(26, 12.0), dpi=300)
gs = fig.add_gridspec(
    3, 2,
    width_ratios=[32, 4.5],
    height_ratios=[2.0, 10.0, 0.95],
    wspace=0.03,
    hspace=0.035,
    left=0.14, right=0.96, top=0.94, bottom=0.22
)

ax_top = fig.add_subplot(gs[0, 0])
ax_main = fig.add_subplot(gs[1, 0], sharex=ax_top)
ax_region = fig.add_subplot(gs[2, 0], sharex=ax_top)
ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
ax_region_right = fig.add_subplot(gs[2, 1])

# --- 1. Top Panel: Number of Methods Deployed per Study ---
ax_top.bar(
    np.arange(num_studies), study_counts,
    color="#334155", width=0.72, edgecolor="none", alpha=0.90
)
for x, val in enumerate(study_counts):
    ax_top.text(x, val + 0.15, str(val), ha="center", va="bottom", fontsize=7.5, fontweight="bold", color="#1E293B")
ax_top.set_ylim(0, 7.5)
ax_top.set_ylabel("Methods (n)", fontsize=10, fontweight="bold", color="#0F172A")
ax_top.set_yticks(range(1, 7))
ax_top.grid(axis="y", color="#E2E8F0", linestyle="--", linewidth=0.5)
ax_top.tick_params(axis="x", labelbottom=False, bottom=False)
ax_top.tick_params(axis="y", labelsize=8.5)
ax_top.spines["top"].set_visible(False)
ax_top.spines["right"].set_visible(False)

# --- 2. Main Panel: Discrete Presence/Absence Matrix ---
for r in range(num_methods):
    for c in range(num_studies):
        if matrix[r, c] == 1:
            rect = patches.Rectangle(
                (c - 0.42, r - 0.42), 0.84, 0.84,
                facecolor=COLORS[r],
                edgecolor="#FFFFFF",
                linewidth=1.0,
                hatch=HATCHES[r]
            )
        else:
            rect = patches.Rectangle(
                (c - 0.42, r - 0.42), 0.84, 0.84,
                facecolor="#F8FAFC",
                edgecolor="#E2E8F0",
                linewidth=0.8
            )
        ax_main.add_patch(rect)

ax_main.set_xlim(-0.6, num_studies - 0.4)
ax_main.set_ylim(-0.6, num_methods - 0.4)
ax_main.set_yticks(range(num_methods))
ax_main.set_yticklabels(METHODS, fontsize=10.5, fontweight="bold", color="#0F172A")
ax_main.tick_params(axis="x", labelbottom=False, bottom=False)
ax_main.tick_params(axis="y", length=3, color="#94A3B8")
ax_main.invert_yaxis() # Top to bottom: Culture first, WGS last
for spine in ax_main.spines.values():
    spine.set_color("#CBD5E1")
    spine.set_linewidth(1.0)

# --- 3. Regional Layer Strip (Dedicated UN Subregion Track) ---
for c in range(num_studies):
    reg = STUDIES[c][3]
    color = REGION_COLORS[reg]
    rect = patches.Rectangle(
        (c - 0.42, -0.42), 0.84, 0.84,
        facecolor=color,
        edgecolor="#FFFFFF",
        linewidth=1.0
    )
    ax_region.add_patch(rect)

ax_region.set_xlim(-0.6, num_studies - 0.4)
ax_region.set_ylim(-0.5, 0.5)
ax_region.set_yticks([0])
ax_region.set_yticklabels(["Region"], fontsize=10.5, fontweight="bold", color="#0F172A")
ax_region.set_xticks(range(num_studies))

# Format study labels: Author et al. (Year) [Country Code]
formatted_study_labels = [f"{s[0]} [{s[2]}]" for s in STUDIES]
ax_region.set_xticklabels(formatted_study_labels, rotation=90, ha="right", fontsize=8.0, fontweight="medium", color="#1E293B")
ax_region.tick_params(axis="both", which="both", length=3, color="#94A3B8")
for spine in ax_region.spines.values():
    spine.set_color("#CBD5E1")
    spine.set_linewidth(1.0)

# --- 4. Right Panel: Total Studies Deploying Each Method ---
y_positions = range(num_methods)
ax_right.barh(
    y_positions, method_counts,
    color=COLORS, height=0.72, edgecolor="#334155", linewidth=0.6, alpha=0.92
)
for r, (val, h) in enumerate(zip(method_counts, HATCHES)):
    bar = ax_right.patches[r]
    bar.set_hatch(h)
    ax_right.text(val + 1.2, r, f"{val} ({val/num_studies*100:.1f}%)", va="center", ha="left", fontsize=8.5, fontweight="bold", color="#1E293B")

ax_right.set_xlim(0, 78)
ax_right.set_xlabel("Studies (k=68)", fontsize=9.5, fontweight="bold", color="#0F172A")
ax_right.grid(axis="x", color="#E2E8F0", linestyle="--", linewidth=0.5)
ax_right.tick_params(axis="y", labelleft=False, left=False)
ax_right.tick_params(axis="x", labelsize=8.5)
ax_right.spines["top"].set_visible(False)
ax_right.spines["right"].set_visible(False)

# --- 5. Right Regional Cell (Clean, no overlapping text) ---
ax_region_right.axis("off")

# --- 6. Styled Subregion Legend Bar (Top Left / Centered above grid) ---
legend_elements = [
    patches.Patch(facecolor=REGION_COLORS["Western Africa"], edgecolor="#FFFFFF", label=f"Western Africa (n = {reg_counter['Western Africa']})"),
    patches.Patch(facecolor=REGION_COLORS["Eastern Africa"], edgecolor="#FFFFFF", label=f"Eastern Africa (n = {reg_counter['Eastern Africa']})"),
    patches.Patch(facecolor=REGION_COLORS["Northern Africa"], edgecolor="#FFFFFF", label=f"Northern Africa (n = {reg_counter['Northern Africa']})"),
    patches.Patch(facecolor=REGION_COLORS["Pan-African Multi-region"], edgecolor="#FFFFFF", label=f"Pan-African (n = {reg_counter['Pan-African Multi-region']})"),
    patches.Patch(facecolor=REGION_COLORS["Southern Africa"], edgecolor="#FFFFFF", label=f"Southern Africa (n = {reg_counter['Southern Africa']})"),
    patches.Patch(facecolor=REGION_COLORS["Central Africa"], edgecolor="#FFFFFF", label=f"Central Africa (n = {reg_counter['Central Africa']})"),
]

leg = ax_top.legend(
    handles=legend_elements,
    loc="lower left",
    bbox_to_anchor=(0.0, 1.05),
    ncol=6,
    fontsize=8.5,
    title="Region:",
    title_fontsize=9.5,
    frameon=True,
    facecolor="#FFFFFF",
    edgecolor="#CBD5E1",
    handlelength=1.2,
    handleheight=1.0
)
leg.get_title().set_fontweight("bold")

plt.savefig(OUT_PNG1, dpi=300)
plt.savefig(OUT_PNG2, dpi=300)
plt.savefig(OUT_PDF1)
plt.savefig(OUT_PDF2)
plt.close()

print(f"SUCCESS: Generated Enhanced Discrete Tile Matrix with UN Subregion Layer for Fig S18!")
print(f"Saved: {OUT_PNG1}")
print(f"Saved: {OUT_PNG2}")
