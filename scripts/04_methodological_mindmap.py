#!/usr/bin/env python3
"""
generate_fig1_mindmap.py
========================
Publication-Grade Minimalist Scientific Evidence & Workflow Mindmap for Figure 1.
Replaces legacy PRISMA flowchart with an exquisite Nature/PLOS/Cell-style architecture.
Resolution: 300 DPI (4,800 x 3,200 px).
"""

import os
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

# Typography setup
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['mathtext.fontset'] = 'dejavusans'

FIG_W = 16.0
FIG_H = 10.6667
DPI = 300

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI, facecolor='#FFFFFF')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

# ------------------------------------------------------------------------------
# PALETTE DEFINITION (Nature / PLOS Minimalist Academic Palette)
# ------------------------------------------------------------------------------
NAVY_DARK   = '#0F172A'  # Slate 900
NAVY_PRIMARY= '#1E3A8A'  # Blue 900
NAVY_LIGHT  = '#EFF6FF'  # Blue 50
NAVY_BORDER = '#93C5FD'  # Blue 300
NAVY_ACCENT = '#2563EB'  # Blue 600

GREEN_DARK  = '#064E3B'  # Emerald 900
GREEN_PRIM  = '#065F46'  # Emerald 800
GREEN_LIGHT = '#F0FDF4'  # Emerald 50
GREEN_BORDER= '#86EFAC'  # Emerald 300
GREEN_ACCENT= '#059669'  # Emerald 600

RED_DARK    = '#7F1D1D'  # Red 900
RED_PRIM    = '#991B1B'  # Red 800
RED_LIGHT   = '#FEF2F2'  # Red 50
RED_BORDER  = '#FCA5A5'  # Red 300
RED_ACCENT  = '#DC2626'  # Red 600

AMBER_DARK  = '#78350F'  # Amber 900
AMBER_PRIM  = '#92400E'  # Amber 800
AMBER_LIGHT = '#FFFBEB'  # Amber 50
AMBER_BORDER= '#FCD34D'  # Amber 300
AMBER_ACCENT= '#D97706'  # Amber 600

SLATE_DARK  = '#1E293B'  # Slate 800
SLATE_MED   = '#475569'  # Slate 600
SLATE_LIGHT = '#94A3B8'  # Slate 400
BG_SUBTLE   = '#F8FAFC'  # Slate 50
BORDER_BASE = '#E2E8F0'  # Slate 200

# ------------------------------------------------------------------------------
# DRAWING HELPERS
# ------------------------------------------------------------------------------
def draw_rounded_card(ax, x, y, w, h, bg_color, border_color, border_width=1.1, radius=0.008, zorder=2):
    box = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=bg_color,
        edgecolor=border_color,
        linewidth=border_width,
        zorder=zorder
    )
    ax.add_patch(box)
    return box

def draw_header_pill(ax, x, y, w, h, bg_color, border_color, title, subtitle=None, title_color='#FFFFFF', sub_color='#E2E8F0'):
    pill = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={h*0.48}",
        facecolor=bg_color,
        edgecolor=border_color,
        linewidth=1.3,
        zorder=5
    )
    ax.add_patch(pill)
    if subtitle:
        ax.text(x + w/2, y + h*0.62, title, ha='center', va='center',
                fontsize=10.8, fontweight='bold', color=title_color, zorder=6)
        ax.text(x + w/2, y + h*0.26, subtitle, ha='center', va='center',
                fontsize=8.0, fontweight='normal', color=sub_color, zorder=6)
    else:
        ax.text(x + w/2, y + h*0.5, title, ha='center', va='center',
                fontsize=10.8, fontweight='bold', color=title_color, zorder=6)

def draw_minimalist_card(ax, x, y, w, h, kicker, title, bullets, theme_color, border_color,
                         wrap_width=32, font_size=7.2, kicker_size=6.3, title_size=8.6):
    """
    Draws a clean, non-colliding scientific card with adaptive line spacing
    guaranteeing that text never crosses the bottom card border.
    """
    # Card Background
    draw_rounded_card(ax, x, y, w, h, bg_color='#FFFFFF', border_color=border_color, border_width=1.1, radius=0.007, zorder=3)
    
    # Top accent bar
    accent_h = 0.004
    ax.add_patch(patches.FancyBboxPatch(
        (x, y + h - accent_h), w, accent_h,
        boxstyle="round,pad=0,rounding_size=0.001",
        facecolor=theme_color, edgecolor='none', zorder=4
    ))
    
    # Kicker (uppercase category)
    ax.text(x + 0.008, y + h - 0.0125, kicker.upper(), ha='left', va='center',
            fontsize=kicker_size, fontweight='bold', color=theme_color, zorder=5)
    
    # Title
    ax.text(x + 0.008, y + h - 0.0235, title, ha='left', va='center',
            fontsize=title_size, fontweight='bold', color=SLATE_DARK, zorder=5)
    
    # Thin divider under title
    ax.plot([x + 0.008, x + w - 0.008], [y + h - 0.030, y + h - 0.030],
            color=BORDER_BASE, lw=0.7, zorder=4)
    
    # Wrap all bullets and gather line objects
    all_lines = []
    for b in bullets:
        wlines = textwrap.wrap(b, width=wrap_width)
        for idx, wl in enumerate(wlines):
            all_lines.append((idx == 0, wl))
            
    # Calculate available vertical height
    y_start = y + h - 0.040
    y_end = y + 0.008
    avail_h = y_start - y_end
    n_lines = len(all_lines)
    
    # Dynamic line spacing bounded safely
    line_spacing = min(0.0118, avail_h / max(n_lines, 1))
    
    curr_y = y_start
    for is_bullet, line_text in all_lines:
        if is_bullet:
            ax.plot(x + 0.010, curr_y, marker='o', markersize=2.4, color=theme_color, zorder=5)
            ax.text(x + 0.016, curr_y, line_text, ha='left', va='center',
                    fontsize=font_size, color=SLATE_DARK, zorder=5)
        else:
            ax.text(x + 0.016, curr_y, line_text, ha='left', va='center',
                    fontsize=font_size, color=SLATE_DARK, zorder=5)
        curr_y -= line_spacing

def draw_spline(ax, p0, p1, p2, p3, color='#94A3B8', lw=1.5, ls='-', alpha=0.75, zorder=1):
    verts = [p0, p1, p2, p3]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', edgecolor=color, lw=lw, linestyle=ls, alpha=alpha, zorder=zorder)
    ax.add_patch(patch)

# ------------------------------------------------------------------------------
# 1. TOP HEADER BANNER
# ------------------------------------------------------------------------------
ax.text(0.50, 0.980, "S Y S T E M A T I C   E V I D E N C E   &   M E T H O D O L O G I C A L   W O R K F L O W   A R C H I T E C T U R E",
        ha='center', va='center', fontsize=9.2, fontweight='bold', color=NAVY_ACCENT)

ax.text(0.50, 0.955, "Global & African Diarrhoeagenic Escherichia coli (DEC) Synthesis: 55-Year Molecular Surveillance (1970–2026)",
        ha='center', va='center', fontsize=16.5, fontweight='bold', color=SLATE_DARK)

ax.text(0.50, 0.933, "Dual-track continental epidemiology (k = 28, N = 15,750), global virulome meta-analysis (k = 47, 29 loci), biostatistical modeling, and policy translation",
        ha='center', va='center', fontsize=10.0, fontstyle='italic', color=SLATE_MED)

# Header divider rule
ax.plot([0.030, 0.970], [0.918, 0.918], color=BORDER_BASE, lw=1.2, zorder=2)

# ------------------------------------------------------------------------------
# 2. CENTRAL CORE (THE SYNTHESIS HUB)
# ------------------------------------------------------------------------------
hub_x, hub_y, hub_w, hub_h = 0.405, 0.415, 0.190, 0.160

# Double card border for subtle depth
draw_rounded_card(ax, hub_x - 0.004, hub_y - 0.004, hub_w + 0.008, hub_h + 0.008,
                  bg_color=NAVY_LIGHT, border_color=NAVY_BORDER, border_width=1.0, radius=0.014, zorder=2)
draw_rounded_card(ax, hub_x, hub_y, hub_w, hub_h,
                  bg_color='#FFFFFF', border_color=NAVY_PRIMARY, border_width=2.0, radius=0.012, zorder=4)

# Hub Header Tag
ax.text(0.50, hub_y + hub_h - 0.018, "S Y N T H E S I S   C O R E", ha='center', va='center',
        fontsize=8.0, fontweight='bold', color=NAVY_ACCENT, zorder=6)

ax.text(0.50, hub_y + hub_h*0.64, "Global & African DEC\nEvidence Architecture", ha='center', va='center',
        fontsize=12.5, fontweight='bold', color=NAVY_DARK, linespacing=1.15, zorder=6)

core_pills = [
    "55-Year Horizon (1970–2026)",
    "k = 75 Studies · N > 25,000 Isolates",
    "REML Random-Effects & GRADE"
]
for i, pill_txt in enumerate(core_pills):
    py = hub_y + 0.048 - i * 0.017
    pill_box = patches.FancyBboxPatch(
        (0.50 - 0.082, py - 0.0065), 0.164, 0.0135,
        boxstyle="round,pad=0,rounding_size=0.006",
        facecolor=NAVY_LIGHT, edgecolor=NAVY_BORDER, lw=0.8, zorder=5
    )
    ax.add_patch(pill_box)
    ax.text(0.50, py, pill_txt, ha='center', va='center',
            fontsize=7.6, fontweight='semibold', color=NAVY_PRIMARY, zorder=6)

# ------------------------------------------------------------------------------
# 3. CLUSTER 1 (TOP-LEFT): PART A — CONTINENTAL SURVEILLANCE & AMR
# ------------------------------------------------------------------------------
c1_header_x, c1_header_y, c1_header_w, c1_header_h = 0.030, 0.865, 0.350, 0.038
draw_header_pill(ax, c1_header_x, c1_header_y, c1_header_w, c1_header_h,
                 bg_color=NAVY_PRIMARY, border_color=NAVY_DARK,
                 title="PART A: CONTINENTAL AFRICAN SURVEILLANCE & AMR",
                 subtitle="k = 28 Primary Cohorts · N = 15,750 Subjects · 18 African Nations (1970–2026)")

col1_x = 0.030
col2_x = 0.212
card_w = 0.168
row1_y = 0.710
row2_y = 0.550
card_h1 = 0.142
card_h2 = 0.148

# Card 1A: Canonical Pathotypes & Hybrids
draw_minimalist_card(
    ax, col1_x, row1_y, card_w, card_h1,
    kicker="Pathotype Diversity",
    title="Canonical Types & Hybrids",
    bullets=[
        "6 canonical classes: EAEC, ETEC, EPEC, DAEC, STEC, EIEC",
        "Emergent hybrids: EAEC/ETEC, EAEC/STEC, DEC/ExPEC",
        "EAEC dominant in persistent infant diarrhea (38.3%)",
        "High asymptomatic carriage challenges attribution"
    ],
    theme_color=NAVY_PRIMARY, border_color=NAVY_BORDER,
    wrap_width=32, font_size=7.2
)

# Card 1B: Subregional Endemic Gradient
draw_minimalist_card(
    ax, col2_x, row1_y, card_w, card_h1,
    kicker="Spatial Distribution",
    title="Subregional Endemic Gradient",
    bullets=[
        "Western Africa (k=12): Peak burden (Nigeria, Ghana, Mali)",
        "Eastern Africa (k=8): Balanced EAEC/EPEC co-circulation",
        "Central (k=2), Southern (k=2), & Northern Africa (k=1)",
        "Multicentre cohorts: GEMS (N=4,110), MAL-ED, VIDA"
    ],
    theme_color=NAVY_PRIMARY, border_color=NAVY_BORDER,
    wrap_width=32, font_size=7.2
)

# Card 1C: Diagnostic Modality Evolution
draw_minimalist_card(
    ax, col1_x, row2_y, card_w, card_h2,
    kicker="Technology Transition",
    title="Diagnostic Modality Evolution",
    bullets=[
        "Phenotypic / HEp-2 Adherence: Historical benchmark (1970–90s)",
        "Conventional & Multiplex PCR: Genotypic shift (2000–2015)",
        "Quantitative TAC (qPCR): High-sensitivity burden (2015+)",
        "Whole-Genome Sequencing: High-resolution phylogenomics"
    ],
    theme_color=NAVY_PRIMARY, border_color=NAVY_BORDER,
    wrap_width=32, font_size=7.2
)

# Card 1D: High-Acuity AMR & Mobilome
draw_minimalist_card(
    ax, col2_x, row2_y, card_w, card_h2,
    kicker="Antimicrobial Resistance",
    title="High-Acuity AMR & Mobilome",
    bullets=[
        "MDR phenotypes: Ampicillin (88%), SXT (79%), Tet (74%)",
        "Cephalosporin resistance driven by mobile blaCTX-M-15",
        "Mobilome core: class 1 integrons (intI1) & IncFII/I1 plasmids",
        "Severe empirical risk for fluoroquinolones & β-lactams"
    ],
    theme_color=NAVY_PRIMARY, border_color=NAVY_BORDER,
    wrap_width=32, font_size=7.2
)

# ------------------------------------------------------------------------------
# 4. CLUSTER 2 (BOTTOM-LEFT): PART B — GLOBAL VIRULOME META-ANALYSIS
# ------------------------------------------------------------------------------
c2_header_x, c2_header_y, c2_header_w, c2_header_h = 0.030, 0.445, 0.350, 0.038
draw_header_pill(ax, c2_header_x, c2_header_y, c2_header_w, c2_header_h,
                 bg_color=GREEN_PRIM, border_color=GREEN_DARK,
                 title="PART B: GLOBAL VIRULOME CASE-CONTROL META-ANALYSIS",
                 subtitle="k = 47 Cohorts · 29 Virulence Loci · 154 Locus Strata · 5 Continents")

row_b1_y = 0.320
row_b2_y = 0.200
row_b3_y = 0.080
card_hb = 0.110

# Card 2A: SPATE Proteases
draw_minimalist_card(
    ax, col1_x, row_b1_y, card_w, card_hb,
    kicker="Serine Proteases",
    title="SPATE Autotransporters",
    bullets=[
        "Class 1 Cytotoxic: sat (mOR=1.82), pet enterotoxin",
        "Class 2 Mucinases: pic (mOR=1.68), sepA protease",
        "Mediates mucosal disruption & immune evasion"
    ],
    theme_color=GREEN_PRIM, border_color=GREEN_BORDER,
    wrap_width=32, font_size=7.1
)

# Card 2B: Aggregative Adherence Fimbriae
draw_minimalist_card(
    ax, col2_x, row_b1_y, card_w, card_hb,
    kicker="Fimbrial Adhesins",
    title="Aggregative Adherence (AAF)",
    bullets=[
        "Fimbrial subunits: aafA (AAF/II), agg3A, agg4A, aggA",
        "Mediates stacked-brick adherence to enterocytes",
        "AAF/II strongly linked to severe dehydration"
    ],
    theme_color=GREEN_PRIM, border_color=GREEN_BORDER,
    wrap_width=32, font_size=7.1
)

# Card 2C: Canonical Master Regulon
draw_minimalist_card(
    ax, col1_x, row_b2_y, card_w, card_hb,
    kicker="Regulon Core",
    title="Master Regulon (pAA)",
    bullets=[
        "aggR master AraC transcriptional activator (mOR=1.56)",
        "aatA ABC transporter (CVD432) & aap dispersin",
        "aaiC Type VI secretion system & astA EAST-1 toxin"
    ],
    theme_color=GREEN_PRIM, border_color=GREEN_BORDER,
    wrap_width=32, font_size=7.1
)

# Card 2D: EPEC Island OI-122 Effectors
draw_minimalist_card(
    ax, col2_x, row_b2_y, card_w, card_hb,
    kicker="Non-LEE Island",
    title="EPEC Island OI-122 Effectors",
    bullets=[
        "lifA / efa1: Lymphocyte inhibition & mucosal adherence",
        "nleB & nleE: NF-κB suppression & tissue necrosis",
        "Highest clinical severity among atypical EPEC"
    ],
    theme_color=GREEN_PRIM, border_color=GREEN_BORDER,
    wrap_width=32, font_size=7.1
)

# Card 2E: DAEC Adhesins & Age Inversion
draw_minimalist_card(
    ax, col1_x, row_b3_y, card_w, card_hb,
    kicker="Adhesin Genetics",
    title="DAEC Adhesins & Age Inversion",
    bullets=[
        "Surface adhesins: daaC, afaE, and sat_DAEC variants",
        "Age-stratified inversion: significant risk >24m",
        "Asymptomatic mucosal colonization in early infancy"
    ],
    theme_color=GREEN_PRIM, border_color=GREEN_BORDER,
    wrap_width=32, font_size=7.1
)

# Card 2F: Siderophore Iron Acquisition
draw_minimalist_card(
    ax, col2_x, row_b3_y, card_w, card_hb,
    kicker="Iron Acquisition",
    title="Siderophore Iron Scavenging",
    bullets=[
        "High-pathogenicity island (HPI): irp2 & fyuA locus",
        "Promotes microaerophilic luminal iron uptake",
        "Correlates with intra-host persistence & fitness"
    ],
    theme_color=GREEN_PRIM, border_color=GREEN_BORDER,
    wrap_width=32, font_size=7.1
)

# ------------------------------------------------------------------------------
# 5. CLUSTER 3 (TOP-RIGHT): BIOSTATISTICAL MODELING & VARIANCE PARTITIONING
# ------------------------------------------------------------------------------
c3_header_x, c3_header_y, c3_header_w, c3_header_h = 0.620, 0.865, 0.350, 0.038
draw_header_pill(ax, c3_header_x, c3_header_y, c3_header_w, c3_header_h,
                 bg_color=RED_PRIM, border_color=RED_DARK,
                 title="BIOSTATISTICAL MODELING & VARIANCE PARTITIONING",
                 subtitle="Freeman–Tukey REML · Continuity Corrections · Mixed Meta-Regression")

col3_x = 0.620
col4_x = 0.802

# Card 3A: Freeman-Tukey REML Pooling
draw_minimalist_card(
    ax, col3_x, row1_y, card_w, card_h1,
    kicker="Random-Effects Pooling",
    title="Freeman–Tukey REML Model",
    bullets=[
        "Double-arcsine transformation stabilizes binomial variance",
        "Restricted Maximum Likelihood (REML) random effects",
        "Heterogeneity diagnostics: Higgins I², Cochran Q, tau²",
        "Subgroup pooling across all 6 DEC pathotypes"
    ],
    theme_color=RED_PRIM, border_color=RED_BORDER,
    wrap_width=32, font_size=7.2
)

# Card 3B: Continuity Corrections
draw_minimalist_card(
    ax, col4_x, row1_y, card_w, card_h1,
    kicker="Sparse Event Handling",
    title="Empirical Continuity (+0.5)",
    bullets=[
        "Haldane–Anscombe correction for zero-event strata",
        "Eliminates undefined log-ORs in sparse 2x2 contingency tables",
        "Preserves study weights without artificial inflation",
        "Benchmarked against profile-likelihood estimators"
    ],
    theme_color=RED_PRIM, border_color=RED_BORDER,
    wrap_width=32, font_size=7.2
)

# Card 3C: Mixed-Effects Meta-Regression
draw_minimalist_card(
    ax, col3_x, row2_y, card_w, card_h2,
    kicker="Variance Partitioning",
    title="Mixed-Effects Meta-Regression",
    bullets=[
        "R² variance partitioning: marginal vs conditional variance",
        "Moderators: Diagnostic Modality, Historical Era, Age",
        "Modality accounts for 64.2% of between-cohort variance",
        "Resolves historical controversies in DEC epidemiology"
    ],
    theme_color=RED_PRIM, border_color=RED_BORDER,
    wrap_width=32, font_size=7.2
)

# Card 3D: Disease Attribution Odds Ratios
draw_minimalist_card(
    ax, col4_x, row2_y, card_w, card_h2,
    kicker="Aetiological Fraction",
    title="Disease Attribution Odds Ratios",
    bullets=[
        "Mantel–Haenszel matched odds ratios (mOR) with 95% CIs",
        "Firth penalized likelihood corrects small-sample bias",
        "Separates causal aetiological agents from bystanders",
        "Stratified attributions delineate infant risk profiles"
    ],
    theme_color=RED_PRIM, border_color=RED_BORDER,
    wrap_width=32, font_size=7.2
)

# ------------------------------------------------------------------------------
# 6. CLUSTER 4 (BOTTOM-RIGHT): EVIDENCE CERTAINTY & POLICY TRANSLATION
# ------------------------------------------------------------------------------
c4_header_x, c4_header_y, c4_header_w, c4_header_h = 0.620, 0.445, 0.350, 0.038
draw_header_pill(ax, c4_header_x, c4_header_y, c4_header_w, c4_header_h,
                 bg_color=AMBER_PRIM, border_color=AMBER_DARK,
                 title="EVIDENCE CERTAINTY & POLICY TRANSLATION",
                 subtitle="GRADE Profiling · Africa CDC PGI Strategy · WHO AFRO IDSR Revision")

# Card 4A: GRADE Multi-Domain Profiling (Full width of cluster)
grade_w = 0.350
grade_h = 0.175
grade_y = 0.255
draw_minimalist_card(
    ax, c4_header_x, grade_y, grade_w, grade_h,
    kicker="Systematic Evidence Quality",
    title="GRADE Multi-Domain Certainty Profiling across 5 Domains",
    bullets=[
        "Domain 1 — Risk of Bias: Newcastle–Ottawa Scale (NOS ≥ 7/9) & JBI critical appraisal",
        "Domain 2 — Inconsistency: Heterogeneity assessment (I² > 75%, Cochran Q p < 0.05)",
        "Domain 3 — Indirectness: Generalizability from tertiary centers to rural primary health",
        "Domain 4 — Imprecision: 95% CI breadth and statistical power sample-size thresholds",
        "Domain 5 — Publication Bias: Egger's regression, Begg's rank test & funnel asymmetry",
        "Overall Evidence Tier: High Certainty (Grade A) for EAEC/ETEC; Moderate (Grade B) for specific loci"
    ],
    theme_color=AMBER_PRIM, border_color=AMBER_BORDER,
    wrap_width=74, font_size=7.1
)

# Bottom Row: Africa CDC PGI and WHO AFRO IDSR side by side
c4_row2_y = 0.080
c4_card_h2 = 0.160

# Card 4B: Africa CDC PGI Architecture
draw_minimalist_card(
    ax, col3_x, c4_row2_y, card_w, c4_card_h2,
    kicker="Continental Genomics",
    title="Africa CDC PGI Architecture",
    bullets=[
        "Standardized enteric genomic surveillance network",
        "Distributed biobanking & bioinformatics pipelines",
        "Genomic tracking of epidemic clones (ST10, ST38)",
        "Plasmid mobilome tracking for blaCTX-M-15",
        "Decentralized sequencing capacity across NPHIs"
    ],
    theme_color=AMBER_PRIM, border_color=AMBER_BORDER,
    wrap_width=32, font_size=7.1
)

# Card 4C: WHO AFRO IDSR Protocol Revision
draw_minimalist_card(
    ax, col4_x, c4_row2_y, card_w, c4_card_h2,
    kicker="Public Health Response",
    title="WHO AFRO IDSR Revision",
    bullets=[
        "Revision of 3rd Edition IDSR guidelines",
        "Mandating multiplex molecular syndromic testing",
        "Differential management for invasive DEC hybrids",
        "Antimicrobial stewardship: antibiotic de-escalation",
        "Early outbreak triggers for community diarrhoeal clusters"
    ],
    theme_color=AMBER_PRIM, border_color=AMBER_BORDER,
    wrap_width=32, font_size=7.1
)

# ------------------------------------------------------------------------------
# 7. ROUTED VECTOR CONNECTORS & METHODOLOGICAL FLOW
# ------------------------------------------------------------------------------
# 1. Bridge to Cluster 2 (Bottom-Left Header):
ax.plot([0.405, 0.380], [c2_header_y + c2_header_h*0.5, c2_header_y + c2_header_h*0.5],
        color=GREEN_ACCENT, lw=2.4, zorder=3)
ax.plot(0.380, c2_header_y + c2_header_h*0.5, marker='<', markersize=6.0, color=GREEN_ACCENT, zorder=4)

# 2. Bridge to Cluster 4 (Bottom-Right Header):
ax.plot([0.595, 0.620], [c4_header_y + c4_header_h*0.5, c4_header_y + c4_header_h*0.5],
        color=AMBER_ACCENT, lw=2.4, zorder=3)
ax.plot(0.620, c4_header_y + c4_header_h*0.5, marker='>', markersize=6.0, color=AMBER_ACCENT, zorder=4)

# 3. Routed Spline to Cluster 1 (Top-Left Header):
draw_spline(ax, (0.425, hub_y + hub_h), (0.405, hub_y + hub_h + 0.10),
            (0.395, 0.82), (c1_header_x + c1_header_w, c1_header_y + c1_header_h*0.5),
            color=NAVY_ACCENT, lw=2.2, alpha=0.9, zorder=3)
ax.plot(c1_header_x + c1_header_w, c1_header_y + c1_header_h*0.5, marker='<', markersize=6.0, color=NAVY_ACCENT, zorder=4)

# 4. Routed Spline to Cluster 3 (Top-Right Header):
draw_spline(ax, (0.575, hub_y + hub_h), (0.595, hub_y + hub_h + 0.10),
            (0.605, 0.82), (c3_header_x, c3_header_y + c3_header_h*0.5),
            color=RED_ACCENT, lw=2.2, alpha=0.9, zorder=3)
ax.plot(c3_header_x, c3_header_y + c3_header_h*0.5, marker='>', markersize=6.0, color=RED_ACCENT, zorder=4)

# Methodological progression dashed connectors:
# Horizontal data integration feed from C1 to C3:
ax.annotate("", xy=(c3_header_x - 0.008, c3_header_y + c3_header_h*0.5),
            xytext=(c1_header_x + c1_header_w + 0.008, c1_header_y + c1_header_h*0.5),
            arrowprops=dict(arrowstyle="->,head_width=0.35,head_length=0.45",
                            color=SLATE_LIGHT, lw=1.5, linestyle='--'),
            zorder=3)
ax.text(0.50, 0.884, "Harmonized Data Integration Feed", ha='center', va='center',
        fontsize=7.8, fontstyle='italic', color=SLATE_MED,
        bbox=dict(boxstyle="round,pad=0.25", facecolor='#FFFFFF', edgecolor=BORDER_BASE, lw=0.8),
        zorder=5)

# Left vertical cross-talk feed (Part A -> Part B virulome alignment):
ax.annotate("", xy=(c2_header_x + c2_header_w*0.5, c2_header_y + c2_header_h + 0.006),
            xytext=(c1_header_x + c1_header_w*0.5, row2_y - 0.006),
            arrowprops=dict(arrowstyle="->,head_width=0.35,head_length=0.45",
                            color=SLATE_LIGHT, lw=1.5, linestyle='--'),
            zorder=3)
ax.text(c1_header_x + c1_header_w*0.5, (c2_header_y + c2_header_h + row2_y)/2,
        "Molecular Typing & Virulome Alignment", ha='center', va='center',
        fontsize=7.2, fontstyle='italic', color=SLATE_MED,
        bbox=dict(boxstyle="round,pad=0.2", facecolor='#FFFFFF', edgecolor=BORDER_BASE, lw=0.8),
        zorder=5)

# Right vertical statistical synthesis feed (Biostatistics -> GRADE & Policy):
ax.annotate("", xy=(c4_header_x + c4_header_w*0.5, c4_header_y + c4_header_h + 0.006),
            xytext=(c3_header_x + c3_header_w*0.5, row2_y - 0.006),
            arrowprops=dict(arrowstyle="->,head_width=0.35,head_length=0.45",
                            color=SLATE_LIGHT, lw=1.5, linestyle='--'),
            zorder=3)
ax.text(c4_header_x + c4_header_w*0.5, (c4_header_y + c4_header_h + row2_y)/2,
        "Statistical Synthesis Feed", ha='center', va='center',
        fontsize=7.2, fontstyle='italic', color=SLATE_MED,
        bbox=dict(boxstyle="round,pad=0.2", facecolor='#FFFFFF', edgecolor=BORDER_BASE, lw=0.8),
        zorder=5)

# ------------------------------------------------------------------------------
# 8. FOOTER / INTEGRATED EVIDENCE METHODOLOGY STRIP
# ------------------------------------------------------------------------------
ax.plot([0.030, 0.970], [0.062, 0.062], color=BORDER_BASE, lw=1.2, zorder=2)

steps = [
    ("Step 1: Systematic Harmonization", "PRISMA PICOS Search (PubMed, Scopus, AJOL, ENA)", NAVY_PRIMARY),
    ("Step 2: Dual-Track Data Stratification", "Part A: African DEC | Part B: Global Virulome", NAVY_ACCENT),
    ("Step 3: Freeman–Tukey REML Modeling", "Double-Arcsine Pooling & Mixed Meta-Regression", RED_PRIM),
    ("Step 4: Multi-Domain GRADE Profiling", "Certainty Appraisal across 5 Rigorous Domains", AMBER_PRIM),
    ("Step 5: Public Health Translation", "Africa CDC Genomics & WHO AFRO IDSR Revision", GREEN_PRIM),
]

n_steps = len(steps)
step_w = 0.176
step_gap = (0.940 - n_steps * step_w) / (n_steps - 1)

for s_idx, (stitle, ssub, scolor) in enumerate(steps):
    sx = 0.030 + s_idx * (step_w + step_gap)
    sy = 0.016
    sh = 0.036
    
    # Step card
    sbox = patches.FancyBboxPatch(
        (sx, sy), step_w, sh,
        boxstyle="round,pad=0,rounding_size=0.006",
        facecolor=BG_SUBTLE, edgecolor=BORDER_BASE, lw=1.0, zorder=3
    )
    ax.add_patch(sbox)
    
    # Colored accent dot
    ax.plot(sx + 0.010, sy + sh*0.5, marker='o', markersize=3.8, color=scolor, zorder=4)
    
    # Text
    ax.text(sx + 0.020, sy + sh*0.68, stitle, ha='left', va='center',
            fontsize=7.8, fontweight='bold', color=SLATE_DARK, zorder=4)
    ax.text(sx + 0.020, sy + sh*0.28, ssub, ha='left', va='center',
            fontsize=6.8, color=SLATE_MED, zorder=4)
    
    # Arrow between steps
    if s_idx < n_steps - 1:
        arr_x1 = sx + step_w + 0.002
        arr_x2 = sx + step_w + step_gap - 0.002
        ax.annotate("", xy=(arr_x2, sy + sh*0.5), xytext=(arr_x1, sy + sh*0.5),
                    arrowprops=dict(arrowstyle="->,head_width=0.25,head_length=0.35",
                                    color=SLATE_LIGHT, lw=1.2),
                    zorder=4)

# ------------------------------------------------------------------------------
# SAVE OUTPUTS
# ------------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_paths = [
    os.path.join(BASE_DIR, "figures/Fig1_PRISMA2020_Flowchart.png"),
    os.path.join(BASE_DIR, "meta_analysis_jegh_preprint/figures/Fig1_PRISMA2020_Flowchart.png"),
    os.path.join(BASE_DIR, "figures/Fig1_Minimalist_Workflow_Mindmap.png"),
    os.path.join(BASE_DIR, "meta_analysis_jegh_preprint/figures/Fig1_Minimalist_Workflow_Mindmap.png"),
    os.path.join(BASE_DIR, "figures/Fig1_Minimalist_Workflow_Mindmap.pdf"),
    os.path.join(BASE_DIR, "meta_analysis_jegh_preprint/figures/Fig1_Minimalist_Workflow_Mindmap.pdf"),
    os.path.join(BASE_DIR, "figures/Fig1_PRISMA2020_Flowchart.pdf"),
    os.path.join(BASE_DIR, "meta_analysis_jegh_preprint/figures/Fig1_PRISMA2020_Flowchart.pdf"),
    os.path.join(BASE_DIR, "figures/Fig1_Minimalist_Workflow_Mindmap.svg"),
    os.path.join(BASE_DIR, "meta_analysis_jegh_preprint/figures/Fig1_Minimalist_Workflow_Mindmap.svg")
]

for p in output_paths:
    os.makedirs(os.path.dirname(p), exist_ok=True)
    fmt = 'svg' if p.endswith('.svg') else ('pdf' if p.endswith('.pdf') else 'png')
    plt.savefig(p, format=fmt, dpi=DPI, bbox_inches=None, facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Saved: {p}")

plt.close()
print("Figure 1 generation completed successfully.")
