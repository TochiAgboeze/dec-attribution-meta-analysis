#!/usr/bin/env python3
"""
clean_and_rebuild_manuscript.py

Rebuilds PLOS_NTD_DEC_Meta_Analysis_Manuscript.md from COMPLETE_JEGH_PREPRINT_MANUSCRIPT.md:
1. Retains strictly:
   - Brief title: # Diarrhoeagenic *Escherichia coli* Attribution Paradox in Africa and Globally: A Two-Part Systematic Review and Meta-Analysis
   - Body from ## Abstract through ## 5. Conclusions
   - # References
2. Strips:
   - Author metadata, affiliations, journal placeholders
   - Abbreviations section
   - Supporting Information file lists
   - Declarations
3. Transforms Table 9 into a publication-grade 6-column GRADE Summary of Findings (SoF) table.
4. Cleans Table 1 pathotype columns (removes redundant denominators).
5. Strips all inline LaTeX math ($) and formatting tags from all table cells, captions, and table notes across all 9 tables.
"""

import os
import re

SOURCE_PATH = "meta_analysis_jegh_preprint/COMPLETE_JEGH_PREPRINT_MANUSCRIPT.md"
DEST_PATH = "meta_analysis_jegh_preprint/PLOS_NTD_DEC_Meta_Analysis_Manuscript.md"

with open(SOURCE_PATH, "r", encoding="utf-8") as f:
    text = f.read()

# 1. Extract body: from ## Abstract to ## 6. Abbreviations
abstract_idx = text.find("## Abstract")
abbrev_idx = text.find("## 6. Abbreviations")
ref_idx = text.find("# References")

if abstract_idx == -1 or abbrev_idx == -1 or ref_idx == -1:
    raise ValueError(f"Indices missing: abstract={abstract_idx}, abbrev={abbrev_idx}, ref={ref_idx}")

body = text[abstract_idx:abbrev_idx].strip()
references = text[ref_idx:].strip()

# Clean table math helper
def clean_math_in_str(s):
    s = s.replace(r"\text{cases}", "cases").replace(r"\text{ctrl}", "ctrl")
    s = s.replace(r"\text{res}", "res").replace(r"\text{mOR}", "mOR")
    s = s.replace(r"\text{OR}", "OR").replace(r"\text{lnOR}", "lnOR")
    s = s.replace(r"\text{--}", "–")
    s = s.replace(r"\tau^2_{\text{res}}", "τ²_res")
    s = s.replace(r"\tau^2", "τ²").replace(r"\tau", "τ")
    s = s.replace(r"I^2", "I²").replace(r"R^2", "R²")
    s = s.replace(r"Q_M", "QM").replace(r"Q_E", "QE")
    s = s.replace(r"\beta", "β").replace(r"\oplus", "⊕").replace(r"\bigcirc", "○")
    s = s.replace(r"\%", "%")
    s = s.replace(r"\ge", "≥").replace(r"\le", "≤").replace(r"\pm", "±")
    s = s.replace("$", "")
    return s

# 2. Refactor Table 9 to 6 columns
t9_match = re.search(r"(\*\*Table 9\..*?\n\s*\n\|.*?)(?=\n\n##|\Z)", body, re.DOTALL)
if not t9_match:
    raise ValueError("Table 9 not found in body")

t9_raw = t9_match.group(1)
t9_lines = t9_raw.split("\n")

t9_caption = clean_math_in_str(t9_lines[0])
new_t9_lines = [
    t9_caption,
    "",
    "| Domain / Target Evaluated | Studies (k) & Design | Sample Size (Cases / Controls) | Pooled Effect Estimate [95% CI] | GRADE Certainty & Domain Appraisal | Clinical Interpretation & Actionability |",
    "| :--- | :---: | :---: | :---: | :--- | :--- |"
]

# Find header row index
pipe_start = 0
for idx, l in enumerate(t9_lines):
    if l.startswith("| :---"):
        pipe_start = idx + 1
        break

for line in t9_lines[pipe_start:]:
    if not line.strip() or not line.startswith("|"):
        continue
    cols = [clean_math_in_str(c.strip()) for c in line.split("|")[1:-1]]
    if len(cols) < 11:
        continue
    
    domain = cols[0]
    # Section header rows
    if "PART A:" in domain or "PART B:" in domain:
        new_t9_lines.append(f"| {domain} | | | | | |")
        continue
    
    k_design = cols[1]
    sample_sz = cols[2]
    rob = cols[3]
    inconsistency = cols[4]
    indirectness = cols[5]
    imprecision = cols[6]
    other = cols[7]
    effect = cols[8]
    grade = cols[9]
    clin = cols[10]
    
    # Synthesize GRADE appraisal
    appraisal_parts = []
    if "serious" in inconsistency.lower():
        appraisal_parts.append(f"Inconsistency: {inconsistency}")
    if "serious" in indirectness.lower():
        appraisal_parts.append(f"Indirectness: {indirectness}")
    if "serious" in imprecision.lower():
        appraisal_parts.append(f"Imprecision: {imprecision}")
    if other and other.lower() != "none":
        appraisal_parts.append(f"Magnitude: {other}")
    if rob and "low" not in rob.lower():
        appraisal_parts.append(f"Bias: {rob}")
    
    if appraisal_parts:
        grade_combined = f"{grade}; {'; '.join(appraisal_parts)}"
    else:
        grade_combined = f"{grade}; No downgrade (Low risk of bias, no inconsistency, indirectness, or imprecision)"
    
    # Ensure square brackets in effect estimate
    effect_clean = effect.replace("(", "[").replace(")", "]") if "[" not in effect else effect
    
    row_str = f"| {domain} | {k_design} | {sample_sz} | {effect_clean} | {grade_combined} | {clin} |"
    new_t9_lines.append(row_str)

new_t9_text = "\n".join(new_t9_lines)

# Replace Table 9 in body
body = body[:t9_match.start()] + new_t9_text + body[t9_match.end():]

# 3. Clean Table 1 pathotype columns (remove /denominator repetition)
t1_match = re.search(r"(\*\*Table 1\..*?\n\s*\n\|.*?)(?=\n\n\*\*|\n\n###|\n\n##|\Z)", body, re.DOTALL)
if t1_match:
    t1_raw = t1_match.group(1)
    t1_lines = t1_raw.split("\n")
    cleaned_t1_lines = []
    for idx, line in enumerate(t1_lines):
        if line.startswith("|") and not line.startswith("| :---") and "Study (Author" not in line:
            cells = line.split("|")
            for c_idx in range(7, 13): # Pathotype cols EAEC..EIEC
                if c_idx < len(cells) - 1:
                    cell_val = cells[c_idx].strip()
                    cell_val = re.sub(r"(\d+)/\d+\s*(\([0-9.]+%\))", r"\1 \2", cell_val)
                    cells[c_idx] = f" {cell_val} "
            cleaned_t1_lines.append("|".join(cells))
        else:
            cleaned_t1_lines.append(line)
    body = body[:t1_match.start()] + "\n".join(cleaned_t1_lines) + body[t1_match.end():]

# 4. Clean math across ALL tables (cells, captions, footnotes)
def replace_table_block(m):
    block = m.group(0)
    lines = block.split("\n")
    cleaned_lines = []
    for l in lines:
        cleaned_lines.append(clean_math_in_str(l))
    return "\n".join(cleaned_lines)

body = re.sub(r"(\*\*Table \d+\..*?\n\s*\n\|[\s\S]*?)(?=\n\n[^\s\|]|\n\n\Z|\Z)", replace_table_block, body)

# 5. Assemble final manuscript
TITLE = "# Diarrhoeagenic *Escherichia coli* Attribution Paradox in Africa and Globally: A Two-Part Systematic Review and Meta-Analysis\n\n"
final_md = TITLE + body + "\n\n---\n\n" + references + "\n"

with open(DEST_PATH, "w", encoding="utf-8") as f:
    f.write(final_md)

print(f"SUCCESS: Rebuilt clean manuscript at {DEST_PATH}")
print(f"Total characters: {len(final_md)}")
print(f"Total lines: {len(final_md.splitlines())}")
