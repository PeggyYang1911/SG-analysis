"""
Per-cell-line treatment comparison: BTZ vs DMSO (DMSO = control) within EACH
RPMI 8226 cell line, on the 5-slice subset.

Unlike the celline_* time-course scripts (x-axis = NaAsO2 minutes), this is a
CATEGORICAL design: for every cell line the BTZ 1 uM / 3 h condition is compared
to that same line's DMSO vehicle control. Grouped bars, one group per cell line,
two bars per group (DMSO, BTZ), mean +/- SD with jittered per-image points.

Counting rule (unchanged): unique Cell_ID (0 = background, excluded); per-cell
punctae = SUM of Punctae_Count over that cell's kept slices; a cell is SG+ if
that sum > 0. Statistical unit = the image (Series), so group n = images.

Stats: within each cell line, BTZ vs DMSO, PAIRED Student's t-test (scipy
ttest_rel). Images are paired by Series number after sorting; the two arms must
have equal n. For the parental line the DMSO arm has more images than BTZ, so a
DMSO_SERIES filter selects which DMSO Series enter the pairing (see CELL_LINES).
Only the paired images feed the bars, points and stats, so the figure and the
p-values always describe the same data; unpaired/dropped images are still listed
in the workbook's Per_image sheet (Used_in_pairing = No) for traceability.

Nothing on disk is modified; reads the already cleaned + 5-slice-subset CSVs and
writes new files only.

Outputs into OUT_DIR:
  Fig_percelline_SGpositive_percent.png / .svg / .pdf  (+ _PPT.svg, Arial)
  Fig_percelline_SG_per_cell.*
  Fig_percelline_SG_per_SGpos_cell.*
  SG_percelline_treatment_summary.xlsx  (Per_image, Per_cell_counts,
                                          Group_means, Stats)
"""
import csv, glob, os, re, sys
import numpy as np
from scipy import stats
import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

from sg_provenance import (PER_CELL_TAIL, cell_rows, parse_dates, slices_used,
                           write_sheet)

# ============================== DATA SOURCE ================================
FOLDER = (r"J:\FATHO_LAB\YANG Pei-Tzu\Biological evaluation\Immunostaining"
          r"\260921_PY_RPMI_BTZ\260921_analysis\V15 analysis\5slices_subset")
OUT_DIR = FOLDER

# (display label, Okabe-Ito colour, DMSO include-substring, BTZ include-substring,
#  DMSO_SERIES, BTZ_SERIES) -- *_SERIES restricts which Series numbers of that arm
# enter the pairing (None = use all). Substrings are unique per (line, treatment);
# the underscore vs hyphen after "R" separates parental / -R / -R-KO / -RY.
#
# *** PARENTAL RPMI 8226 SLIDE SWAP (per PY, corrected acquisition metadata) ***
# The slides were mislabelled at the microscope: for the parental line the file
# group named "RPMI_DMSO_3h" (10 images, HIGH SG+) is actually the BTZ-treated
# sample, and "RPMI_BTZ_1_uM_3hr" (5 images, LOW SG+) is actually the DMSO
# vehicle control. So the substrings are deliberately crossed here: DMSO<-RPMI_BTZ,
# BTZ<-RPMI_DMSO. The 10-image arm (now BTZ) is thinned to Series 6-10 to pair 5v5
# with the 5-image DMSO arm -- the same physical images selected before the swap.
# In the workbook, Source_CSV therefore says "..._DMSO_..." on rows whose
# Treatment is BTZ (and vice versa): that mismatch is the swap, not a bug.
CELL_LINES = [
    ("RPMI 8226", "#0072B2", "RPMI_BTZ",       "RPMI_DMSO",     None, set(range(6, 11))),
    ("RPMI-R",    "#009E73", "RPMI-R_DMSO",    "RPMI-R_BTZ",    None, None),
    ("RPMI-R-KO", "#CC79A7", "RPMI-R-KO_DMSO", "RPMI-R-KO_BTZ", None, None),
    ("RPMI-RY",   "#E69F00", "RPMI-RY_DMSO",   "RPMI_RY_BTZ",   None, None),
]
TREATMENTS = ["DMSO", "BTZ"]   # DMSO first = control
DSER = {c[0]: c[4] for c in CELL_LINES}
BSER = {c[0]: c[5] for c in CELL_LINES}

METRICS = [
    ("sgpos",   "SG+ cells (%)",   "Fig_percelline_SGpositive_percent", (0, 101),
     [0, 20, 40, 60, 80, 100]),
    ("percell", "SG per cell",     "Fig_percelline_SG_per_cell",        None, None),
    ("persg",   "SG per SG+ cell", "Fig_percelline_SG_per_SGpos_cell",  None, None),
]

# ================================= CONFIG ===================================
BAR_WIDTH   = 0.38
GROUP_GAP   = 1.0            # x-distance between cell-line groups
DMSO_FACE   = "white"       # control bar drawn open (white face, coloured edge)
EDGE_WIDTH  = 2.2
FIGSIZE     = (9.6, 6.8)
FONT_FAMILY = "DejaVu Sans"
FONT_BASE   = 20
FONT_AXLABEL = 24
FONT_TICKS  = 20
ERRORBAR    = "SD"          # "SD" or "SEM"
BW_MODE     = False         # True = grayscale bars (print / B&W): black edges,
                            # white-hatched DMSO vs solid grey BTZ, cell-line
                            # colours ignored. Set on the module before calling
                            # figure() to render a B&W version.
BW_EDGE     = "black"
BW_BTZ_FACE = "0.45"        # grey fill for the BTZ bar in B&W mode
SHOW_POINTS = True
POINT_SIZE  = 34
POINT_ALPHA = 0.5
POINT_JITTER = 0.09         # in x-units
CAPSIZE     = 4
SPINE_WIDTH = 1.6
FORMATS     = ("png", "svg", "pdf")
DPI         = 200
OUT_XLSX    = "SG_percelline_treatment_summary.xlsx"
PPT_SVG     = True
PPT_FONT    = "Arial"
PPT_SUFFIX  = "_PPT"
# ---- SG+ (%) as a LINE / slope figure --------------------------------------
SGPOS_AS_LINE = True        # draw the SG+ (%) metric as connected DMSO->BTZ lines
                            # (one per cell line) instead of grouped bars. The
                            # other two metrics stay as bars.
LINE_LW       = 3.0
LINE_MS       = 11
LINE_XOFFSET  = 0.05        # horizontal spread between cell-line lines (x-units)
LINE_JITTER   = 0.018       # per-image point jitter around each treatment
LEGEND_TITLE  = "Cell line (BTZ vs DMSO)"
# ============================================================================

# indices inside the per-image record tuple
# rec = (line, trt, image, SG+%, SG/cell, SG/SG+cell, cell_rows, slices,
#        n_slices, csv, acq, ana, series)
MIDX = {"sgpos": 3, "percell": 4, "persg": 5}
SERIES_RE = re.compile(r"Series(\d+)", re.IGNORECASE)


def is_data(row, ci):
    try:
        int(float(row[ci])); return True
    except (ValueError, IndexError):
        return False


def image_stats(path):
    """-> (image, SG+%, SG/cell, SG/SG+cell, cell_rows, slices, n_slices) or None."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return None
    h = rows[0]
    try:
        ci = h.index("Cell_ID"); pi = h.index("Punctae_Count"); si = h.index("Slice")
    except ValueError:
        return None
    data = [r for r in rows[1:] if is_data(r, ci)]
    if not data:
        return None
    cellp, cell_sl = {}, {}
    for r in data:
        cid = r[ci].strip()
        if cid == "0":
            continue
        cellp[cid] = cellp.get(cid, 0) + int(float(r[pi]))
        cell_sl.setdefault(cid, set()).add(r[si].strip())
    if not cellp:
        return None
    ncell = len(cellp)
    withp = sum(1 for v in cellp.values() if v > 0)
    tot = sum(cellp.values())
    sl_str, n_sl = slices_used(cell_sl)
    return (data[0][0], withp / ncell * 100, tot / ncell,
            tot / withp if withp else 0.0, cell_rows(cellp, cell_sl), sl_str, n_sl)


def classify(fn):
    """-> (line_label, treatment) or (None, None). First matching (line,trt) wins."""
    for label, _c, dsub, bsub, _ds, _bs in CELL_LINES:
        if dsub in fn:
            return label, "DMSO"
        if bsub in fn:
            return label, "BTZ"
    return None, None


def series_no(fn):
    m = SERIES_RE.search(fn)
    return int(m.group(1)) if m else 10 ** 9


def collect():
    recs = []
    for path in sorted(glob.glob(os.path.join(FOLDER, "*.csv"))):
        fn = os.path.basename(path)
        line, trt = classify(fn)
        if line is None:
            print("  unmatched (skipped):", fn)
            continue
        st = image_stats(path)
        if st is None:
            print("  no data (skipped):", fn)
            continue
        ana, acq = parse_dates(fn)
        recs.append((line, trt) + st + (fn, acq, ana, series_no(fn)))
    return recs


def paired(recs, line, key):
    """Aligned (dmso_vals, btz_vals) for one line+metric, paired by sorted Series.

    DMSO is filtered by DSER[line] first; both arms are sorted by Series number
    and truncated to the shorter length so the pairing is 1:1.
    """
    i = MIDX[key]
    dser, bser = DSER[line], BSER[line]
    d = sorted([r for r in recs if r[0] == line and r[1] == "DMSO"
                and (dser is None or r[12] in dser)], key=lambda r: r[12])
    b = sorted([r for r in recs if r[0] == line and r[1] == "BTZ"
                and (bser is None or r[12] in bser)], key=lambda r: r[12])
    n = min(len(d), len(b))
    d, b = d[:n], b[:n]
    return np.array([x[i] for x in d]), np.array([x[i] for x in b])


def used_in_pairing(recs):
    """Set of Source_CSV names that actually enter the paired analysis."""
    used = set()
    for line, *_ in CELL_LINES:
        dser, bser = DSER[line], BSER[line]
        d = sorted([r for r in recs if r[0] == line and r[1] == "DMSO"
                    and (dser is None or r[12] in dser)], key=lambda r: r[12])
        b = sorted([r for r in recs if r[0] == line and r[1] == "BTZ"
                    and (bser is None or r[12] in bser)], key=lambda r: r[12])
        n = min(len(d), len(b))
        for r in d[:n] + b[:n]:
            used.add(r[9])
    return used


def sig(p):
    if not np.isfinite(p):
        return "n/a"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def ptest(recs, line, key):
    d, b = paired(recs, line, key)
    if len(d) >= 2 and len(d) == len(b):
        t, p = stats.ttest_rel(b, d)
        return float(t), float(p), len(d)
    return float("nan"), float("nan"), len(d)


def figure(recs, key, ylabel, stem, ylim, yticks, font_family=None):
    fam = font_family or FONT_FAMILY
    plt.rcParams.update({"font.size": FONT_BASE, "font.family": fam,
                         "mathtext.fontset": "custom", "mathtext.rm": fam,
                         "mathtext.it": f"{fam}:italic", "mathtext.bf": f"{fam}:bold",
                         "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
                         "savefig.dpi": DPI, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(figsize=FIGSIZE)
    rng = np.random.default_rng(0)
    centers = np.arange(len(CELL_LINES)) * GROUP_GAP
    off = {"DMSO": -BAR_WIDTH / 2 - 0.02, "BTZ": BAR_WIDTH / 2 + 0.02}
    top = 0.0

    for gi, (label, color, _d, _b, _ds, _bs) in enumerate(CELL_LINES):
        d_vals, b_vals = paired(recs, label, key)
        vals = {"DMSO": d_vals, "BTZ": b_vals}
        pair_max = 0.0
        for trt in TREATMENTS:
            v = vals[trt]
            if not len(v):
                continue
            x = centers[gi] + off[trt]
            m = v.mean()
            e = v.std(ddof=1) if len(v) > 1 else 0.0
            if ERRORBAR.upper() == "SEM" and len(v) > 1:
                e = e / np.sqrt(len(v))
            edge = BW_EDGE if BW_MODE else color
            if trt == "DMSO":
                face = DMSO_FACE
            else:
                face = BW_BTZ_FACE if BW_MODE else color
            ax.bar(x, m, BAR_WIDTH, color=face, edgecolor=edge,
                   linewidth=EDGE_WIDTH, zorder=2,
                   hatch=("//" if trt == "DMSO" else None))
            ax.errorbar(x, m, yerr=e, fmt="none", ecolor="0.2",
                        capsize=CAPSIZE, capthick=2.0, lw=2.0, zorder=4)
            if SHOW_POINTS:
                jit = rng.uniform(-POINT_JITTER, POINT_JITTER, size=len(v))
                ax.scatter(np.full(len(v), x) + jit, v, s=POINT_SIZE,
                           color="0.15", alpha=POINT_ALPHA, edgecolors="none",
                           zorder=5)
            pair_max = max(pair_max, m + e, v.max())
        # paired t-test bracket
        _t, p, n = ptest(recs, label, key)
        if np.isfinite(p):
            bump = 2.0 if key == "sgpos" else pair_max * 0.02
            y = pair_max * 1.06 + bump
            x0, x1 = centers[gi] + off["DMSO"], centers[gi] + off["BTZ"]
            ax.plot([x0, x0, x1, x1], [y, y + bump, y + bump, y], lw=1.4,
                    color="0.2", zorder=6)
            ax.text((x0 + x1) / 2, y + bump * 1.2, sig(p), ha="center",
                    va="bottom", fontsize=FONT_BASE)
            top = max(top, y + bump * (5 if key == "sgpos" else 6))
        top = max(top, pair_max)

    ax.set_xticks(centers)
    ax.set_xticklabels([l[0] for l in CELL_LINES])
    ax.set_ylabel(ylabel, fontsize=FONT_AXLABEL)
    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        ax.set_ylim(0, top * 1.12 if top else 1)
    if yticks is not None:
        ax.set_yticks(yticks)
    ax.tick_params(labelsize=FONT_TICKS, width=SPINE_WIDTH, length=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(SPINE_WIDTH)

    handles = [
        mpatches.Patch(facecolor="white", edgecolor="0.2", hatch="//",
                       linewidth=EDGE_WIDTH, label="DMSO (control)"),
        mpatches.Patch(facecolor="0.55", edgecolor="0.2",
                       linewidth=EDGE_WIDTH, label="BTZ 1 $\\mu$M, 3 h"),
    ]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.01),
              fontsize=FONT_BASE - 2, frameon=False, ncol=2, handlelength=1.4,
              borderaxespad=0.0)
    fig.tight_layout()
    return fig, ax


def line_figure(recs, key, ylabel, stem, ylim, yticks, font_family=None):
    """SG+ (%) as connected DMSO->BTZ lines, one per cell line (paired design)."""
    fam = font_family or FONT_FAMILY
    plt.rcParams.update({"font.size": FONT_BASE, "font.family": fam,
                         "mathtext.fontset": "custom", "mathtext.rm": fam,
                         "mathtext.it": f"{fam}:italic", "mathtext.bf": f"{fam}:bold",
                         "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
                         "savefig.dpi": DPI, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(figsize=FIGSIZE)
    rng = np.random.default_rng(0)
    xpos = {"DMSO": 0.0, "BTZ": 1.0}
    n_lines = len(CELL_LINES)
    offs = np.linspace(-LINE_XOFFSET * (n_lines - 1) / 2,
                       LINE_XOFFSET * (n_lines - 1) / 2, n_lines)

    for gi, (label, color, _d, _b, _ds, _bs) in enumerate(CELL_LINES):
        d_vals, b_vals = paired(recs, label, key)
        dx = offs[gi]
        means, errs, xs = [], [], []
        for trt, v in (("DMSO", d_vals), ("BTZ", b_vals)):
            if not len(v):
                continue
            x = xpos[trt] + dx
            m = v.mean()
            e = v.std(ddof=1) if len(v) > 1 else 0.0
            if ERRORBAR.upper() == "SEM" and len(v) > 1:
                e = e / np.sqrt(len(v))
            xs.append(x); means.append(m); errs.append(e)
            if SHOW_POINTS:
                jit = rng.uniform(-LINE_JITTER, LINE_JITTER, size=len(v))
                ax.scatter(np.full(len(v), x) + jit, v, s=POINT_SIZE,
                           color=color, alpha=POINT_ALPHA, edgecolors="none",
                           zorder=3)
        # paired-test significance goes into the legend label (endpoints are
        # too crowded to annotate individually)
        _t, p, _n = ptest(recs, label, key)
        lab = f"{label} ({sig(p)})"
        ax.errorbar(xs, means, yerr=errs, marker="o", ms=LINE_MS, lw=LINE_LW,
                    color=color, capsize=CAPSIZE, capthick=2.0, label=lab,
                    zorder=4)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["DMSO", "BTZ 1 $\\mu$M, 3 h"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel(ylabel, fontsize=FONT_AXLABEL)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if yticks is not None:
        ax.set_yticks(yticks)
    ax.tick_params(labelsize=FONT_TICKS, width=SPINE_WIDTH, length=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(SPINE_WIDTH)

    leg = ax.legend(title=LEGEND_TITLE, loc="lower left",
                    bbox_to_anchor=(0.0, 1.01), fontsize=FONT_BASE - 3,
                    title_fontsize=FONT_BASE - 1, frameon=False, ncol=2,
                    handlelength=1.4, borderaxespad=0.0)
    leg._legend_box.align = "left"
    fig.tight_layout()
    return fig, ax


def save_all(recs):
    for key, ylabel, stem, ylim, yticks in METRICS:
        builder = line_figure if (key == "sgpos" and SGPOS_AS_LINE) else figure
        fig, _ = builder(recs, key, ylabel, stem, ylim, yticks)
        for ext in FORMATS:
            fig.savefig(os.path.join(OUT_DIR, f"{stem}.{ext}"), dpi=DPI,
                        bbox_inches="tight")
        plt.close(fig)
        print("Saved:", os.path.join(OUT_DIR, f"{stem}.png"), "(+svg/pdf)")
        if PPT_SVG:
            pf, _ = builder(recs, key, ylabel, stem, ylim, yticks, font_family=PPT_FONT)
            pf.savefig(os.path.join(OUT_DIR, f"{stem}{PPT_SUFFIX}.svg"), dpi=DPI,
                       bbox_inches="tight")
            plt.close(pf)


def workbook(recs):
    wb = openpyxl.Workbook()
    order = {l[0]: i for i, l in enumerate(CELL_LINES)}
    ordered = sorted(recs, key=lambda r: (order[r[0]], r[1] != "DMSO", r[12]))
    used = used_in_pairing(recs)

    ws = wb.active; ws.title = "Per_image"
    ws.append(["Cell_line", "Treatment", "Series", "Used_in_pairing", "Image",
               "Source_CSV", "Acquisition_date", "Analysis_date", "Slices_used",
               "N_slices", "SG_positive_(%)", "SG_per_cell", "SG_per_SG+cell"])
    for r in ordered:
        ws.append([r[0], r[1], r[12], "Yes" if r[9] in used else "No", r[2], r[9],
                   r[10], r[11], r[7], r[8], round(r[3], 3), round(r[4], 4),
                   round(r[5], 4)])
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    for i, w in enumerate([12, 11, 8, 15, 46, 46, 17, 15, 16, 10, 15, 14, 16], start=1):
        ws.column_dimensions[get_column_letter(i)].width = w

    write_sheet(
        wb.create_sheet("Per_cell_counts"),
        [("Cell_line", "line", 12), ("Treatment", "trt", 11),
         ("Image", "image", 46), ("Source_CSV", "csv", 46)] + PER_CELL_TAIL,
        [dict(c, line=r[0], trt=r[1], image=r[2], csv=r[9], acq_date=r[10],
              ana_date=r[11]) for r in ordered if r[9] in used for c in r[6]])

    ws2 = wb.create_sheet("Group_means")
    ws2.append(["Cell_line", "Treatment", "n_paired", "Mean_SG+_(%)", "SD_SG+_(%)",
                "Mean_SG_per_cell", "SD_SG_per_cell", "Mean_SG_per_SG+cell",
                "SD_SG_per_SG+cell"])
    for label, *_ in CELL_LINES:
        dsg, bsg = paired(recs, label, "sgpos")
        dpc, bpc = paired(recs, label, "percell")
        dps, bps = paired(recs, label, "persg")
        for trt, a1, a2, a3 in (("DMSO", dsg, dpc, dps), ("BTZ", bsg, bpc, bps)):
            if not len(a1):
                continue
            def ms(a):
                return (round(a.mean(), 4),
                        round(a.std(ddof=1), 4) if len(a) > 1 else 0.0)
            m1, s1 = ms(a1); m2, s2 = ms(a2); m3, s3 = ms(a3)
            ws2.append([label, trt, len(a1), m1, s1, m2, s2, m3, s3])
    ws2.freeze_panes = "A2"
    ws2.auto_filter.ref = ws2.dimensions
    for i, w in enumerate([12, 11, 10, 14, 14, 16, 16, 18, 18], start=1):
        ws2.column_dimensions[get_column_letter(i)].width = w

    ws3 = wb.create_sheet("Stats")
    ws3.append(["Cell_line", "Metric", "n_pairs", "Mean_DMSO", "Mean_BTZ",
                "t", "p_value", "Significance"])
    hdr_font = Font(bold=True, color="FFFFFF"); hdr_fill = PatternFill("solid", fgColor="404040")
    for c in range(1, 9):
        ws3.cell(1, c).font = hdr_font; ws3.cell(1, c).fill = hdr_fill
        ws3.cell(1, c).alignment = Alignment(horizontal="center")
    for label, *_ in CELL_LINES:
        for key, ylabel, *_ in METRICS:
            d, b = paired(recs, label, key)
            t, p, n = ptest(recs, label, key)
            ws3.append([label, ylabel, n,
                        round(d.mean(), 4) if len(d) else None,
                        round(b.mean(), 4) if len(b) else None,
                        round(t, 4) if np.isfinite(t) else None,
                        round(p, 6) if np.isfinite(p) else None, sig(p)])
    ws3.append([])
    ws3.append(["Test: PAIRED Student's t-test (ttest_rel), BTZ vs DMSO control, within each cell line."])
    ws3.append(["Pairing: by Series number after sorting; unit = image."])
    ws3.append(["PARENTAL RPMI 8226 SLIDE SWAP: files named 'RPMI_DMSO_3h' are the BTZ sample (high SG+), "
                "'RPMI_BTZ_1_uM_3hr' are the DMSO control (low SG+). Treatment column reflects the corrected label, "
                "not the filename. BTZ arm thinned to Series 6-10 to pair 5v5 with DMSO."])
    ws3.append(["*** p<0.001, ** p<0.01, * p<0.05, ns = not significant."])
    ws3.freeze_panes = "A2"
    for i, w in enumerate([12, 18, 9, 12, 12, 10, 12, 13], start=1):
        ws3.column_dimensions[get_column_letter(i)].width = w

    out = os.path.join(OUT_DIR, OUT_XLSX)
    wb.save(out)
    return out


def main():
    recs = collect()
    print(f"\n{'cell line':12} {'DMSO(paired)':>16} {'BTZ(paired)':>16}")
    for label, *_ in CELL_LINES:
        d, b = paired(recs, label, "sgpos")
        print(f"{label:12} {f'{d.mean():5.1f}% (n={len(d)})':>16} "
              f"{f'{b.mean():5.1f}% (n={len(b)})':>16}")
    print("\n=== BTZ vs DMSO, PAIRED t-test, within each cell line ===")
    for label, *_ in CELL_LINES:
        for key, ylabel, *_ in METRICS:
            t, p, n = ptest(recs, label, key)
            print(f"  {label:11} {ylabel:16} n={n}  p={p:.5f}  {sig(p)}")
    save_all(recs)
    print("\nSaved:", workbook(recs))


if __name__ == "__main__":
    main()
