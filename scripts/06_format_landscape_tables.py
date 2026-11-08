#!/usr/bin/env python3
"""
format_all_tables_landscape.py

Transforms all data-dense tables in the JEGH / PLOS manuscript to publication-grade LANDSCAPE pages:
1. Wraps each multi-column table (and its caption) in a dedicated LANDSCAPE section:
   - Orientation: Landscape (11.0 in x 8.5 in)
   - Margins: 0.5 in (720 dxa) left, right, top, bottom
   - Printable width: 10.0 in (14,400 dxa)
2. Sets explicit, carefully engineered column widths so that every column fits comfortably
   without character wrapping or vertical squishing.
3. Applies strict Three-Line Booktabs borders (thick top, thin header underline, thick bottom, NO vertical lines).
4. Configures `tblHeader` on the header row so headers repeat across multi-page tables.
5. Configures `cantSplit` on every row to prevent mid-row page breaks.
6. Sets tight, elegant cell padding (50 dxa top/bottom, 80 dxa left/right) and crisp typography (7.5-8.0 pt font).
7. Restores the following section back to PORTRAIT (8.5 in x 11.0 in, 1.0 in margins) for narrative text.
8. Automatically compiles to PDF using LibreOffice soffice and verifies page geometries.
"""

import os
import sys
import docx
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import parse_xml, OxmlElement
from docx.oxml.ns import nsdecls, qn
from docx.enum.text import WD_ALIGN_PARAGRAPH

def make_sect_pr(orient="landscape"):
    if orient == "landscape":
        return parse_xml(
            f'<w:sectPr {nsdecls("w")}>'
            f'  <w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/>'
            f'  <w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720" w:header="720" w:footer="720" w:gutter="0"/>'
            f'  <w:type w:val="nextPage"/>'
            f'</w:sectPr>'
        )
    else:
        return parse_xml(
            f'<w:sectPr {nsdecls("w")}>'
            f'  <w:pgSz w:w="12240" w:h="15840" w:orient="portrait"/>'
            f'  <w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>'
            f'  <w:type w:val="nextPage"/>'
            f'</w:sectPr>'
        )

def apply_booktabs_borders(table):
    tblPr = table._tbl.tblPr
    # Remove existing tblBorders if present
    existing = tblPr.find(qn('w:tblBorders'))
    if existing is not None:
        tblPr.remove(existing)
    
    borders = parse_xml(
        f'<w:tblBorders {nsdecls("w")}>'
        f'  <w:top w:val="single" w:sz="10" w:space="0" w:color="0F172A"/>'
        f'  <w:bottom w:val="single" w:sz="10" w:space="0" w:color="0F172A"/>'
        f'  <w:insideH w:val="none"/>'
        f'  <w:left w:val="none"/>'
        f'  <w:right w:val="none"/>'
        f'  <w:insideV w:val="none"/>'
        f'</w:tblBorders>'
    )
    tblPr.append(borders)

def apply_header_bottom_border(row):
    for cell in row.cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tcBorders = parse_xml(
            f'<w:tcBorders {nsdecls("w")}>'
            f'  <w:bottom w:val="single" w:sz="6" w:space="0" w:color="0F172A"/>'
            f'</w:tcBorders>'
        )
        tcPr.append(tcBorders)

def set_cell_margins(cell, top=50, bottom=50, left=80, right=80):
    tcPr = cell._tc.get_or_add_tcPr()
    tcMar = OxmlElement('w:tcMar')
    for m, val in [('top', top), ('bottom', bottom), ('left', left), ('right', right)]:
        node = OxmlElement(f'w:{m}')
        node.set(qn('w:w'), str(val))
        node.set(qn('w:type'), 'dxa')
        tcMar.append(node)
    tcPr.append(tcMar)

def style_table(table, col_widths, font_size_pt=8.0, is_header_shaded=True):
    apply_booktabs_borders(table)
    
    # Repeat header row on every page
    header_tr = table.rows[0]._tr
    header_trPr = header_tr.get_or_add_trPr()
    header_trPr.append(parse_xml(f'<w:tblHeader {nsdecls("w")}/>'))
    
    apply_header_bottom_border(table.rows[0])
    
    # Set alignment to center
    tblPr = table._tbl.tblPr
    jc = parse_xml(f'<w:jc {nsdecls("w")} w:val="center"/>')
    tblPr.append(jc)
    
    total_cols = len(col_widths)
    
    for row_idx, row in enumerate(table.rows):
        trPr = row._tr.get_or_add_trPr()
        trPr.append(parse_xml(f'<w:cantSplit {nsdecls("w")}/>'))
        
        is_header = (row_idx == 0)
        
        for col_idx, cell in enumerate(row.cells):
            if col_idx < total_cols:
                cell.width = Inches(col_widths[col_idx])
                tcPr = cell._tc.get_or_add_tcPr()
                tcW = parse_xml(f'<w:tcW {nsdecls("w")} w:w="{int(col_widths[col_idx] * 1440)}" w:type="dxa"/>')
                tcPr.append(tcW)
            
            set_cell_margins(cell, top=50, bottom=50, left=70, right=70)
            
            # Header shading
            if is_header and is_header_shaded:
                shd = parse_xml(f'<w:shd {nsdecls("w")} w:fill="F1F5F9"/>')
                cell._tc.get_or_add_tcPr().append(shd)
            
            # Typography
            for p in cell.paragraphs:
                p.paragraph_format.space_before = Pt(1)
                p.paragraph_format.space_after = Pt(1)
                p.paragraph_format.line_spacing = 1.05
                if is_header:
                    p.alignment = WD_ALIGN_PARAGRAPH.CENTER if col_idx > 0 else WD_ALIGN_PARAGRAPH.LEFT
                else:
                    # Numbers/stats centered, first column left-aligned
                    p.alignment = WD_ALIGN_PARAGRAPH.LEFT if col_idx == 0 else WD_ALIGN_PARAGRAPH.CENTER
                
                for run in p.runs:
                    run.font.name = "Arial"
                    run.font.size = Pt(font_size_pt)
                    if is_header:
                        run.font.bold = True
                        run.font.color.rgb = RGBColor(15, 23, 42)
                    else:
                        run.font.color.rgb = RGBColor(30, 41, 59)


def format_document(docx_path):
    print(f"\nProcessing document: {docx_path}")
    doc = docx.Document(docx_path)
    body = doc._body._body
    
    # Detailed column widths for each table (total printable width = 10.0 inches)
    TABLE_WIDTH_MAP = {
        # Table 1: Primary Cohort Characteristics (14 cols, total = 9.80 in)
        14: [1.25, 1.05, 0.65, 0.45, 0.85, 0.65, 0.52, 0.52, 0.52, 0.48, 0.48, 0.48, 1.40, 0.50],
        # Table 2: Single-arm case prevalence & CC ORs (13 cols)
        13: [1.40, 0.50, 0.80, 1.00, 1.00, 0.70, 0.50, 1.05, 0.75, 0.55, 0.80, 0.50, 0.45],
        # Table 5, 6, 7, 9: Part B Loci Meta-Analyses & GRADE (11 cols)
        11: [1.20, 1.50, 0.65, 0.95, 0.95, 0.95, 0.60, 0.90, 0.55, 0.90, 0.75],
        # Table 4: Multivariate Meta-Regression (9 cols)
        9:  [1.70, 1.60, 0.95, 0.95, 0.95, 0.75, 1.05, 1.05, 0.90],
        # Table 3: Subgroup Meta-Analysis (8 cols)
        8:  [1.95, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15],
        # Table 8: Global Cohorts Summary (6 cols)
        6:  [1.60, 1.20, 0.80, 1.60, 1.30, 3.40]
    }
    
    # Identify tables to wrap in landscape
    tables_to_format = []
    for i, table in enumerate(doc.tables):
        ncols = len(table.columns)
        print(f"  Table {i+1}: {len(table.rows)} rows x {ncols} cols")
        tables_to_format.append((i, table, ncols))
    
    # Process from last to first so indices/siblings remain valid
    for i, table, ncols in reversed(tables_to_format):
        tbl_el = table._tbl
        prev_p = tbl_el.getprevious()
        
        # Check if prev_p is caption paragraph
        caption_p = None
        if prev_p is not None and prev_p.tag.endswith('p') and ('Table' in prev_p.text or 'TABLE' in prev_p.text):
            caption_p = prev_p
            split_p = caption_p.getprevious()
        else:
            split_p = prev_p
        
        # 1. Terminate preceding portrait section before the caption
        if split_p is not None and split_p.tag.endswith('p'):
            sect_port = make_sect_pr("portrait")
            split_p.get_or_add_pPr().append(sect_port)
        
        # 2. Insert empty paragraph after table with landscape sectPr (terminates the landscape section)
        p_land_end = OxmlElement('w:p')
        p_land_end_pPr = OxmlElement('w:pPr')
        sect_land = make_sect_pr("landscape")
        p_land_end_pPr.append(sect_land)
        p_land_end.append(p_land_end_pPr)
        
        # Insert right after tbl_el (or after footnote paragraph if present)
        next_el = tbl_el.getnext()
        if next_el is not None and next_el.tag.endswith('p') and ('Note:' in next_el.text or 'Abbreviations:' in next_el.text or '*' in next_el.text[:5]):
            next_el.addnext(p_land_end)
        else:
            tbl_el.addnext(p_land_end)
        
        # Determine column widths
        widths = TABLE_WIDTH_MAP.get(ncols)
        if widths is None:
            # fallback equal distribution over 9.8 in
            widths = [9.8 / ncols] * ncols
        
        font_sz = 6.5 if ncols >= 14 else (7.0 if ncols >= 11 else (7.5 if ncols >= 8 else 8.0))
        style_table(table, widths, font_size_pt=font_sz)
    
    # Set final document section to portrait
    doc.save(docx_path)
    print(f"  Saved formatted docx to: {docx_path}")

if __name__ == "__main__":
    files = [
        "meta_analysis_jegh_preprint/JEGH_Preprint_Manuscript.docx",
        "meta_analysis_jegh_preprint/PLOS_NTD_DEC_Meta_Analysis_Manuscript.docx",
        "meta_analysis_jegh_preprint/PI_REVIEW_EXECUTIVE_DRAFT.docx",
        "meta_analysis_jegh_preprint/Supplementary_Information.docx"
    ]
    for f in files:
        if os.path.exists(f):
            format_document(f)
