"""
Step 2 - Stress-granule (SG / punctae) analysis for a folder of cleaned v14 CSVs.

For each *.csv (one image = one Series) it:
  * appends a SUMMARY block to the CSV: cells with punctae (SG+), total cells,
    percent SG+, total punctae in SG+ cells, mean punctae per cell,
    mean punctae per SG+ cell.  (Re-running regenerates the block, no dupes.)

It then aggregates ALL images into a new Excel workbook (SG_analysis_summary.xlsx)
with sheets:
  Summary_per_image, Group_means, Punctae_per_cell, Statistics, Figures

Summary_per_image and Punctae_per_cell carry full provenance (source CSV, imaging
and analysis date, Series, Cell_ID, and which z-slices each count came from), so a
number can always be traced back to the image it came from and per-cell rows from
different experiments can be pooled into one table later.
and writes two figures (also embedded in the workbook):
  Figure1_SGpositive_percent.png       - SG+ cells (%) per treatment group
  Figure2_SG_per_cell.png              - number of SG per cell (all cells)
  Figure3_SG_per_SGpositive_cell.png   - number of SG per SG+ cell only

Statistics (statistical unit = per image / Series):
  * Paired Student's t-test of each treatment vs the 0-min control. Images are
    paired IN SORTED ORDER (Series numbers are not meaningful) and each
    comparison is truncated to the smaller group's size.
  * One-way ANOVA across all timepoints for SG-per-cell.

Grouping is auto-detected from the filenames:
  * TIME mode        -> names with "_<N>_min_Series<M>_" : numeric timepoints,
                        control = 0 min.
  * CATEGORICAL mode -> otherwise, group = first keyword found from CAT_KEYWORDS
                        (e.g. DMSO / BTZ), control = the first keyword.
Edit CAT_KEYWORDS (keep the control first) to add new treatment types.

Usage:
    python sg_analysis.py "PATH\\TO\\v14 analysis"

Run clean_csv.py FIRST if the CSVs still contain Cell_ID=0 rows / tiny slices.
Requires: numpy, scipy, matplotlib, openpyxl.
"""
import csv
import glob
import os
import re
import sys
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpyxl
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

from sg_provenance import cell_rows, parse_dates, slices_used

CURRENT_FOLDER = r"J:\FATHO_LAB\YANG Pei-Tzu\Biological evaluation\Immunostaining\260518_PY_ Sterallis autofluorescent 260508 RPMI 1-45 min As\Extracted_Hyperstacks\RPMI\v14 analysis"

# ---- Grouping is auto-detected per folder from the file names ----
# TIME mode  : name contains "<N> min" (any separator: "_0_min_", " 15 min-",
#              etc.) -> numeric timepoints, control = 0 min, x-axis = As exposure.
# CATEGORICAL: otherwise, group = first matching keyword in CAT_KEYWORDS
#              (control = the FIRST keyword).
# Replicate number = "Series<M>" if present, else the "-<M>_v14" trailing index.
# Add keywords to CAT_KEYWORDS as new treatment types appear; keep the intended
# control first.
TIME_RE = re.compile(r"(\d+)\s*_?\s*min", re.IGNORECASE)   # "0_min", " 15 min"
SERIES_RE = re.compile(r"Series(\d+)")
REP_RE = re.compile(r"-(\d+)_v14")                          # "...-01_v14"
CAT_KEYWORDS = ["DMSO", "BTZ"]   # display order; FIRST entry is the control

TIME_CONTROL, TIME_UNIT, TIME_XLABEL = 0, "min", "Sodium arsenite exposure"
CAT_CONTROL, CAT_UNIT, CAT_XLABEL = CAT_KEYWORDS[0], "", "Treatment"


def _replicate(fname):
    sm = SERIES_RE.search(fname) or REP_RE.search(fname)
    return int(sm.group(1)) if sm else -1


def parse_name(fname):
    """Return (treatment, series); treatment is int (time) or str (categorical)."""
    tm = TIME_RE.search(fname)
    if tm:
        return int(tm.group(1)), _replicate(fname)
    up = fname.upper()
    grp = next((k for k in CAT_KEYWORDS if k.upper() in up), "unknown")
    return grp, _replicate(fname)


def sig_mark(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"


def analyze_folder(folder):
    files = sorted(glob.glob(os.path.join(folder, "*.csv")))
    print(f"Found {len(files)} CSV files")

    per_image = []
    per_cell = []

    for path in files:
        fname = os.path.basename(path)
        treatment, series = parse_name(fname)

        with open(path, "r", newline="") as f:
            rows = list(csv.reader(f))
        header = rows[0]
        cell_idx = header.index("Cell_ID")

        # Detect data rows by an integer Cell_ID (NOT by a ".tif" File name, which
        # would drop .czi and other datasets). Skip files with no parseable data.
        def is_data(r):
            if len(r) < len(header):
                return False
            try:
                int(float(r[cell_idx]))
                return True
            except (ValueError, IndexError):
                return False

        data = [r for r in rows[1:] if is_data(r)]
        if not data:
            print(f"  skipping (no parseable data rows): {fname}")
            continue

        punc_idx = header.index("Punctae_Count")
        slice_idx = header.index("Slice")
        # A cell (Cell_ID) is measured once per z-slice. Aggregate to ONE value
        # per unique Cell_ID, summing Punctae_Count across that cell's slices, so
        # total cells = number of unique Cell_IDs (Cell_ID 0 = background, dropped)
        # and a cell is SG+ if it has any punctae. (Row-based counting would count
        # each cell once per slice and is wrong.)
        cell_punct = {}
        cell_slices = {}          # Cell_ID -> set of slices that cell appears in
        for r in data:
            cid = r[cell_idx].strip()
            if cid == "0":
                continue
            cell_punct[cid] = cell_punct.get(cid, 0) + int(float(r[punc_idx]))
            cell_slices.setdefault(cid, set()).add(r[slice_idx].strip())
        percell = list(cell_punct.values())
        if not percell:
            print(f"  skipping (no cells after excluding Cell_ID 0): {fname}")
            continue
        total_cells = len(percell)
        cells_with = sum(1 for v in percell if v > 0)
        total_punctae = sum(percell)
        pct = (cells_with / total_cells * 100) if total_cells else 0.0
        mean_per_cell = (total_punctae / total_cells) if total_cells else 0.0
        mean_per_pos = (total_punctae / cells_with) if cells_with else 0.0
        image_name = data[0][0] if data else fname.replace("_v14.csv", ".tif")
        ana_date, acq_date = parse_dates(fname)
        img_slices_str, n_img_slices = slices_used(cell_slices)

        out_rows = [header] + data
        out_rows.append([])
        out_rows.append(["SUMMARY"])
        out_rows.append(["Cells_with_punctae_(SG+)", cells_with])
        out_rows.append(["Total_cells", total_cells])
        out_rows.append(["Percent_SG_positive_cells_(%)", round(pct, 4)])
        out_rows.append(["Total_punctae_in_SG+_cells", total_punctae])
        out_rows.append(["Mean_punctae_per_cell", round(mean_per_cell, 4)])
        out_rows.append(["Mean_punctae_per_SG+_cell", round(mean_per_pos, 4)])
        out_rows.append(["Slices_used", img_slices_str])
        out_rows.append(["N_slices", n_img_slices])
        out_rows.append(["Acquisition_date", acq_date])
        out_rows.append(["Analysis_date", ana_date])
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows(out_rows)

        per_image.append({
            "image": image_name, "csv": fname,
            "acq_date": acq_date, "ana_date": ana_date,
            "treatment": treatment, "series": series,
            "slices": img_slices_str, "n_slices": n_img_slices,
            "total_cells": total_cells, "cells_with": cells_with,
            "pct": pct, "total_punctae": total_punctae,
            "mean_per_cell": mean_per_cell, "mean_per_pos": mean_per_pos,
        })
        # One row per counted cell, keeping enough provenance to trace the count
        # back to its image and slices, and to pool cells across experiments.
        for c in cell_rows(cell_punct, cell_slices):
            c.update(image=image_name, csv=fname,
                     acq_date=acq_date, ana_date=ana_date,
                     treatment=treatment, series=series)
            per_cell.append(c)

    numeric = all(isinstance(d["treatment"], int) for d in per_image)
    if numeric:
        CONTROL, GROUP_UNIT = TIME_CONTROL, TIME_UNIT
        # x-axis drug label auto-detected from the file names
        names = " ".join(d["image"] for d in per_image).upper()
        if "BTZ" in names or "BORTEZOMIB" in names:
            X_AXIS_LABEL = "Bortezomib (BTZ) exposure"
        else:
            X_AXIS_LABEL = TIME_XLABEL
        per_image.sort(key=lambda d: (d["treatment"], d["series"]))
        treatments = sorted({d["treatment"] for d in per_image})
    else:
        CONTROL, GROUP_UNIT, X_AXIS_LABEL = CAT_CONTROL, CAT_UNIT, CAT_XLABEL
        order = {k: i for i, k in enumerate(CAT_KEYWORDS)}
        per_image.sort(key=lambda d: (order.get(d["treatment"], 99), d["series"]))
        treatments = sorted({d["treatment"] for d in per_image},
                            key=lambda t: order.get(t, 99))

    def lab(t):
        return f"{t} {GROUP_UNIT}".strip()
    ctrl_label = lab(CONTROL)

    print("Treatments:", treatments)
    for t in treatments:
        n = sum(1 for d in per_image if d["treatment"] == t)
        print(f"  {lab(t)}: {n} images")

    def group_vals(key, t):
        return [d[key] for d in per_image if d["treatment"] == t]

    def paired_vs_control(key):
        res = {}
        ctrl = group_vals(key, CONTROL)
        for t in treatments:
            if t == CONTROL:
                continue
            other = group_vals(key, t)
            n = min(len(ctrl), len(other))
            tstat, p = stats.ttest_rel(np.array(ctrl[:n]), np.array(other[:n]))
            res[t] = (tstat, p, n)
        return res

    pct_tests = paired_vs_control("pct")
    spc_tests = paired_vs_control("mean_per_cell")
    pos_tests = paired_vs_control("mean_per_pos")
    def anova(key):
        groups = [g for g in (group_vals(key, t) for t in treatments) if len(g) > 0]
        if len(groups) < 2:
            return float("nan"), float("nan")
        return stats.f_oneway(*groups)
    F, p_anova = anova("mean_per_cell")
    F2, p_anova2 = anova("mean_per_pos")

    def summarize(key):
        means, sds, ns, allvals = [], [], [], []
        for t in treatments:
            v = np.array(group_vals(key, t))
            means.append(v.mean())
            sds.append(v.std(ddof=1) if len(v) > 1 else 0.0)
            ns.append(len(v)); allvals.append(v)
        return means, sds, ns, allvals

    labels = [lab(t) for t in treatments]
    xpos = np.arange(len(treatments))
    rng = np.random.default_rng(0)

    def make_barfig(key, tests, color, ylabel, title, outname):
        means, sds, ns, allvals = summarize(key)
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        ax.bar(xpos, means, yerr=sds, capsize=5, color=color,
               edgecolor="black", alpha=0.85, zorder=2)
        for i, v in enumerate(allvals):
            jit = rng.uniform(-0.12, 0.12, size=len(v))
            ax.scatter(xpos[i] + jit, v, color="black", s=20, zorder=3)
        ax.set_xticks(xpos); ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel); ax.set_xlabel(X_AXIS_LABEL)
        ax.set_title(title)
        top = max(m + s for m, s in zip(means, sds))
        top = max(top, max((max(v) if len(v) else 0) for v in allvals))
        step = top * 0.10 if top else 1.0
        base = top + step
        ctrl_i = treatments.index(CONTROL) if CONTROL in treatments else 0
        others = [i for i, t in enumerate(treatments) if t != CONTROL]
        for level, i in enumerate(others):
            _, p, _ = tests[treatments[i]]
            y = base + level * step
            ax.plot([xpos[ctrl_i], xpos[ctrl_i], xpos[i], xpos[i]],
                    [y, y + step*0.25, y + step*0.25, y], lw=1.2, c="black")
            ax.text((xpos[ctrl_i] + xpos[i]) / 2, y + step*0.28,
                    sig_mark(p), ha="center", va="bottom", fontsize=11)
        ax.set_ylim(0, base + len(others) * step + step)
        fig.tight_layout()
        out = os.path.join(folder, outname)
        fig.savefig(out, dpi=200); plt.close(fig)
        return out, means, sds, ns

    fig1, means, sds, ns = make_barfig(
        "pct", pct_tests, "#4C72B0", "SG+ cells (%)",
        f"Stress-granule-positive cells (%)\n(paired t-test vs {ctrl_label})",
        "Figure1_SGpositive_percent.png")
    fig2, means2, sds2, ns2 = make_barfig(
        "mean_per_cell", spc_tests, "#C44E52", "Number of SG per cell",
        f"SG per cell (one-way ANOVA p = {p_anova:.3g})\n(paired t-test vs {ctrl_label})",
        "Figure2_SG_per_cell.png")
    fig3, means3, sds3, ns3 = make_barfig(
        "mean_per_pos", pos_tests, "#55A868", "SG per SG+ cell",
        f"SG per SG+ cell (one-way ANOVA p = {p_anova2:.3g})\n(paired t-test vs {ctrl_label})",
        "Figure3_SG_per_SGpositive_cell.png")

    # ---- Excel workbook ----
    wb = openpyxl.Workbook()
    hdr_font = Font(bold=True, color="FFFFFF")
    hdr_fill = PatternFill("solid", fgColor="404040")

    def style_header(ws, ncol, row=1):
        for c in range(1, ncol + 1):
            cell = ws.cell(row=row, column=c)
            cell.font = hdr_font; cell.fill = hdr_fill
            cell.alignment = Alignment(horizontal="center")

    ws = wb.active; ws.title = "Summary_per_image"
    cols = ["Image", "Source_CSV", "Acquisition_date", "Analysis_date",
            "Treatment", "Series", "Slices_used", "N_slices",
            "Total_cells", "Cells_with_punctae_(SG+)",
            "SG_positive_(%)", "Total_punctae_in_SG+_cells", "Mean_punctae_per_cell",
            "Mean_punctae_per_SG+_cell"]
    ws.append(cols)
    for d in per_image:
        ws.append([d["image"], d["csv"], d["acq_date"], d["ana_date"],
                   d["treatment"], d["series"], d["slices"], d["n_slices"],
                   d["total_cells"], d["cells_with"], round(d["pct"], 3),
                   d["total_punctae"],
                   round(d["mean_per_cell"], 4), round(d["mean_per_pos"], 4)])
    style_header(ws, len(cols)); ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    for i, w in enumerate([46, 46, 17, 15, 12, 8, 16, 10, 12, 22, 14, 24, 20, 22],
                          start=1):
        ws.column_dimensions[get_column_letter(i)].width = w

    ws2 = wb.create_sheet("Group_means")
    ws2.append(["Treatment", "n_images", "Mean_SG+_(%)", "SD_SG+_(%)",
                "Mean_SG_per_cell", "SD_SG_per_cell",
                "Mean_SG_per_SG+cell", "SD_SG_per_SG+cell"])
    for i, t in enumerate(treatments):
        ws2.append([t, ns[i], round(means[i], 3), round(sds[i], 3),
                    round(means2[i], 4), round(sds2[i], 4),
                    round(means3[i], 4), round(sds3[i], 4)])
    style_header(ws2, 8)
    for i, w in enumerate([12, 10, 14, 14, 18, 18, 20, 20], start=1):
        ws2.column_dimensions[get_column_letter(i)].width = w

    ws3 = wb.create_sheet("Punctae_per_cell")
    c3 = ["Image", "Source_CSV", "Acquisition_date", "Analysis_date", "Treatment",
          "Series", "Cell_ID", "Slices", "N_slices", "Punctae_per_cell",
          "SG_positive"]
    ws3.append(c3)
    for d in per_cell:
        ws3.append([d["image"], d["csv"], d["acq_date"], d["ana_date"],
                    d["treatment"], d["series"], d["cell_id"], d["slices"],
                    d["n_slices"], d["punctae"], d["sg_pos"]])
    style_header(ws3, len(c3)); ws3.freeze_panes = "A2"
    ws3.auto_filter.ref = ws3.dimensions
    for i, w in enumerate([46, 46, 17, 15, 12, 8, 9, 16, 10, 18, 12], start=1):
        ws3.column_dimensions[get_column_letter(i)].width = w

    ws4 = wb.create_sheet("Statistics")
    ws4.cell(1, 1, f"Paired Student's t-test vs {ctrl_label}  "
                   "(per-image; images paired in order, n = smaller group)").font = Font(bold=True)
    ws4.append(["Metric", "Comparison", "n_pairs", "t", "p_value", "Significance"])
    style_header(ws4, 6, row=2)
    for t in treatments:
        if t == CONTROL: continue
        tstat, p, n = pct_tests[t]
        ws4.append(["SG+ cells (%)", f"{ctrl_label} vs {lab(t)}", n,
                    round(tstat, 4), round(p, 6), sig_mark(p)])
    for t in treatments:
        if t == CONTROL: continue
        tstat, p, n = spc_tests[t]
        ws4.append(["SG per cell", f"{ctrl_label} vs {lab(t)}", n,
                    round(tstat, 4), round(p, 6), sig_mark(p)])
    for t in treatments:
        if t == CONTROL: continue
        tstat, p, n = pos_tests[t]
        ws4.append(["SG per SG+ cell", f"{ctrl_label} vs {lab(t)}", n,
                    round(tstat, 4), round(p, 6), sig_mark(p)])
    ws4.append([])
    ws4.append(["One-way ANOVA (SG per cell)", "F", round(F, 4),
                "p_value", round(p_anova, 6), sig_mark(p_anova)])
    ws4.append(["One-way ANOVA (SG per SG+ cell)", "F", round(F2, 4),
                "p_value", round(p_anova2, 6), sig_mark(p_anova2)])
    ws4.append([])
    ws4.append(["Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant"])
    for i, w in enumerate([28, 18, 10, 12, 12, 14], start=1):
        ws4.column_dimensions[get_column_letter(i)].width = w

    ws5 = wb.create_sheet("Figures")
    ws5.add_image(XLImage(fig1), "A1")
    ws5.add_image(XLImage(fig2), "L1")
    ws5.add_image(XLImage(fig3), "W1")

    out_xlsx = os.path.join(folder, "SG_analysis_summary.xlsx")
    wb.save(out_xlsx)

    print("\n=== Statistics ===")
    print("SG+ (%) paired t vs control:",
          {t: (round(pct_tests[t][1], 6), sig_mark(pct_tests[t][1])) for t in pct_tests})
    print("SG/cell paired t vs control:",
          {t: (round(spc_tests[t][1], 6), sig_mark(spc_tests[t][1])) for t in spc_tests})
    print("SG/SG+cell paired t vs control:",
          {t: (round(pos_tests[t][1], 6), sig_mark(pos_tests[t][1])) for t in pos_tests})
    print(f"ANOVA SG/cell: F={F:.4f}, p={p_anova:.5g} ({sig_mark(p_anova)})")
    print(f"ANOVA SG/SG+cell: F={F2:.4f}, p={p_anova2:.5g} ({sig_mark(p_anova2)})")
    print("\nSaved:")
    for f in (out_xlsx, fig1, fig2, fig3):
        print(" ", f)


if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else CURRENT_FOLDER
    analyze_folder(folder)
