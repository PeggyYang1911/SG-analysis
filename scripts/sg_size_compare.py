"""
Cell size and nucleus size in SG+ vs SG- cells, across the NaAsO2 time course.

Written for the RPMI 8226 panel (parental / HDAC6 KO / BTZ-resistant) in
    260624_PY_260611_RPMI_KO_R_AS\\v15 analysis\\5slices_subset
but it takes the folder as argv[1] and derives everything else from it.

------------------------------------------------------------------------------
WHERE THE TWO SIZES COME FROM
------------------------------------------------------------------------------
Cell size    `Cyto_Area_px2` in the per-image CSVs. Despite the column name the
             values are in MICRONS SQUARED, not pixels: a cell's Cyto_Area is a
             fraction of its own nucleus area measured in pixels, which is
             impossible, and lands on ~1x the cell footprint measured on the SG
             channel once the ROI pixel counts are converted with the TIFF
             calibration. Cyto_Area is one value per z-slice, so the per-cell
             cell size used here is the MEAN over that cell's subset slices.

Nucleus size NOT in the CSVs. It is recovered from the Fiji macro's own
             `<Series>_v14_ROIs.zip`, whose `Cell_<N>.roi` entries are the
             traced DAPI outlines that defined Cell_ID N (verified by overlaying
             them on channel 0 of the extracted hyperstacks). Area is the
             shoelace area of the traced polygon -- ImageJ traced ROIs store
             pixel-CORNER coordinates, so this is the exact enclosed pixel count
             -- times the pixel area read from the hyperstack TIFF calibration.

             Caveat, and it is a real one: the macro saves ONE outline per cell,
             not one per z-slice, so nucleus area is a single-plane measurement
             while cell area is a multi-slice mean. Nucleus area is therefore
             sensitive to which plane the macro segmented on, and should be read
             as a comparison BETWEEN SG classes within an image (which is what
             the statistics below do) rather than as an absolute size.

------------------------------------------------------------------------------
COUNTING AND STATISTICS
------------------------------------------------------------------------------
Counting rule (unchanged): a cell is a unique Cell_ID, Cell_ID 0 is background
and is dropped, per-cell punctae = SUM of Punctae_Count over that cell's slices,
and the cell is SG+ if that sum > 0.

Statistical unit is the IMAGE, not the cell. For every image the SG+ cells and
the SG- cells are each averaged, giving one pair of numbers per image, and the
two are compared with a PAIRED t-test (plus Wilcoxon signed-rank) at each
cell line x timepoint. Pairing is by image, so illumination, focus and density
differences between images cancel out. An image contributes a pair only if it
has at least MIN_CELLS_PER_CLASS cells on both sides; timepoints where one class
is nearly absent (0 min has few SG+ cells, 45 min few SG-) will therefore show a
small n_pairs, which is reported rather than hidden.

The per-cell Mann-Whitney U on the pooled cells is also reported, clearly marked
DESCRIPTIVE: it treats cells as independent and so understates its own p-values.

------------------------------------------------------------------------------
WHY THE ABSOLUTE SIZES MUST NOT BE COMPARED ACROSS TIMEPOINTS
------------------------------------------------------------------------------
In this dataset the apparent cell AND nucleus area track the depth of the z
stack, which was not held constant between the .lif files: RPMI parental was
imaged with 5 z-slices at 0 min and 12 at 30 min, and its median nucleus area
falls from ~80 to ~30 um2 over the same span, with all six series of a timepoint
agreeing closely. Nuclei do not shrink threefold in half an hour. Deeper stacks
place the macro's segmentation plane further from the cell equator, so the 2D
cross-section it outlines is smaller. The effect is per-.lif, i.e. per timepoint.

Consequence: the ABSOLUTE cell/nucleus area curves are a property of the
acquisition as much as of the biology and should not be read as a time course.
Everything that compares SG+ with SG- WITHIN an image is unaffected, because
both classes come from the same stack -- that is the paired test above, and the
SG+/SG- RATIO per image, which is the quantity that can honestly be compared
across timepoints. The ratio figure is the one to use for the time course.

Outputs (into the CSV folder unless OUT_DIR is changed):
  SG_size_analysis.xlsx   Per_cell / Per_image / Group_means / Statistics /
                          Ratio_stats / Provenance / Figures
  Figure_CellArea_SGpos_vs_SGneg.(png|svg|pdf|_PPT.svg)
  Figure_NucleusArea_SGpos_vs_SGneg.(png|svg|pdf|_PPT.svg)
  Figure_NCratio_SGpos_vs_SGneg.(png|svg|pdf|_PPT.svg)
  Figure_SizeRatio_SGpos_over_SGneg.(png|svg|pdf|_PPT.svg)   <- z-independent

Usage:
    python sg_size_compare.py "PATH\\TO\\5slices_subset"
    python sg_size_compare.py "PATH\\TO\\5slices_subset" --show

Requires: numpy, scipy, matplotlib, openpyxl, pillow.
"""
import csv, glob, math, os, re, struct, sys, textwrap, zipfile
import numpy as np
from scipy import stats
import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpyxl
from openpyxl.drawing.image import Image as XLImage
from PIL import Image

from sg_provenance import parse_dates, slice_sort, write_sheet

CURRENT_FOLDER = (r"J:\FATHO_LAB\YANG Pei-Tzu\Biological evaluation\Immunostaining"
                  r"\260624_PY_260611_RPMI_KO_R_AS\v15 analysis\5slices_subset")

# ================================= CONFIG ===================================
# ---- where the non-CSV inputs live, relative to the CSV folder --------------
ROI_DIR_REL   = ".."                       # holds <Series>_v14_ROIs.zip
TIFF_DIR_REL  = os.path.join("..", "..", "Extracted_Hyperstacks")
OUT_DIR       = None                       # None = write next to the CSVs
PIXEL_UM_FALLBACK = 0.3611                 # used only if a TIFF is unreadable

# ---- grouping --------------------------------------------------------------
# (legend label, colour, include-substring). The most specific pattern wins, so
# "_RPMI_AS_" cannot swallow "_RPMI_KO_AS_".
CELL_LINES = [
    ("Parental",      "#0072B2", "_RPMI_AS_"),
    ("HDAC6 KO",      "#E69F00", "_RPMI_KO_AS_"),
    ("BTZ-resistant", "#009E73", "_RPMI_R_AS_"),
]
TIMEPOINTS = [0, 15, 30, 45]
TIME_RE   = re.compile(r"(\d+)\s*_?\s*min", re.IGNORECASE)
SERIES_RE = re.compile(r"Series(\d+)")

# ---- how a cell's size is reduced to ONE number ----------------------------
# Cyto_Area is measured once per z-slice, so each cell has several values.
# "max"  : the largest cross-section the cell shows across its slices. That is
#          the slice nearest its equator, so it estimates the cell's true width
#          and is the least sensitive to WHICH slices happened to be sampled.
# "mean" : average over the cell's slices. Smoother, but it drags the value down
#          whenever the subset includes slices near the top or bottom of a cell,
#          and how often that happens depends on the stack depth.
# Nucleus size is unaffected either way: the macro saves one outline per cell,
# so there is only ever one value to take.
CELL_AREA_STAT = "max"

# ---- statistics ------------------------------------------------------------
MIN_CELLS_PER_CLASS = 3     # an image needs >= this many SG+ AND SG- cells to
                            # contribute a pair at its timepoint

# ---- axes ------------------------------------------------------------------
XLABEL    = "1 mM NaAsO$_2$ (min)"
XLIM      = (-5, 50)
XTICKS    = [0, 15, 30, 45]
YLIM_CELL = None            # None = auto
YLIM_NUC  = None
YLIM_NC   = None

# ---- SG classes ------------------------------------------------------------
SG_NEG_LABEL, SG_NEG_COLOR = "SG\u2212", "#0072B2"    # Okabe-Ito blue
SG_POS_LABEL, SG_POS_COLOR = "SG+",      "#D55E00"    # Okabe-Ito vermillion

# ---- sizes / fonts ---------------------------------------------------------
FIGSIZE      = (15.0, 5.4)   # 3 panels side by side
FONT_FAMILY  = "DejaVu Sans"
FONT_BASE    = 13
FONT_AXLABEL = 15
FONT_TICKS   = 13
FONT_LEGEND  = 13
FONT_TITLE   = 15

# ---- marks -----------------------------------------------------------------
ERRORBAR     = "SD"          # "SD" or "SEM", across images
LINEWIDTH    = 2.4
MARKERSIZE   = 9
CAPSIZE      = 4
SHOW_POINTS  = True          # jittered per-image means
POINT_SIZE   = 34
POINT_ALPHA  = 0.35
POINT_JITTER = 1.2           # x-units (minutes)
SPINE_WIDTH  = 1.4
HIDE_TOP_RIGHT_SPINES = True
STAR_ON_PLOT = True          # significance mark above each timepoint

# Every figure has a caption (see figure_captions()). The captions are always
# written to CAPTIONS_TXT beside the figures; SHOW_FOOTNOTE additionally prints
# each one under its own axes. Off by default, so the figures drop straight into
# a slide and the wording lives in one editable text file instead.
SHOW_FOOTNOTE = False
CAPTIONS_TXT  = "Figure_captions.txt"
FONT_FOOTNOTE = 10.5

# ---- output ----------------------------------------------------------------
FORMATS    = ("png", "svg", "pdf")
DPI        = 200
OUT_XLSX   = "SG_size_analysis.xlsx"
PPT_SVG    = True
PPT_FONT   = "Arial"
PPT_SUFFIX = "_PPT"
# ============================================================================

# (record key, y-axis label, short name, y-limits, output stem)
METRICS = [
    ("cell_area", "Cell area (\u00b5m$^2$)",    "Cell size",
     YLIM_CELL, "Figure_CellArea_SGpos_vs_SGneg"),
    ("nuc_area",  "Nucleus area (\u00b5m$^2$)", "Nucleus size",
     YLIM_NUC,  "Figure_NucleusArea_SGpos_vs_SGneg"),
    ("nc_ratio",  "Nucleus / cell area",        "N/C ratio",
     YLIM_NC,   "Figure_NCratio_SGpos_vs_SGneg"),
]

# ---- the z-independent time course: SG+ / SG- within each image ------------
# (record key, legend label, colour)
RATIO_SERIES = [
    ("cell_area", "Cell area",    "#0072B2"),
    ("nuc_area",  "Nucleus area", "#D55E00"),
    ("nc_ratio",  "N/C ratio",    "#009E73"),
]
RATIO_STEM  = "Figure_SizeRatio_SGpos_over_SGneg"
RATIO_YLIM  = None
RATIO_DX    = {"cell_area": -1.1, "nuc_area": 0.0, "nc_ratio": 1.1}  # x offsets

# ---- the two-group bar figure: SG- and SG+ side by side at each timepoint ---
# Two columns per timepoint, all four timepoints shown. A bar is drawn from
# whatever cells that class actually has, so a group with no cells at a
# timepoint (typically SG+ at 0 min) simply has no bar there rather than
# suppressing the whole timepoint.
# (record key, axis label, y-limits)
BAR_METRICS = [
    ("cell_area", "Cell area (µm$^2$)",    None),
    ("nuc_area",  "Nucleus area (µm$^2$)", None),
]
BAR_STEM      = "Figure_TwoGroup_bars_SGpos_vs_SGneg"
BAR_WIDTH     = 0.38
BAR_ALPHA     = 0.85
BAR_FIGSIZE   = (15.0, 8.8)

# "cell"  : every scored cell in the group contributes. Bars are mean +/- SD
#           over CELLS, n is cells, and the p-value is an unpaired Welch t-test
#           over all cells. Nothing is discarded for want of a partner.
# "image" : bars are mean +/- SD over per-image means and the p-value is the
#           paired test, which needs MIN_CELLS_PER_CLASS cells of both classes
#           in the same image and therefore drops images.
#
# Cell-level is the default here by request: maximise the cells behind each bar.
# It is the more optimistic of the two -- cells inside one image share focus,
# staining and z depth, so they are not independent and the per-cell p-value
# understates its own uncertainty. Both tests are written to the workbook
# (Group_stats_allcells vs Statistics) so the two can be compared directly.
# Within a timepoint this does NOT reintroduce the z-depth confound: all six
# series of a timepoint come from one .lif at one stack depth.
BAR_UNIT = "cell"


# ------------------------------- ImageJ ROI ---------------------------------
def read_roi(buf):
    """(x, y) vertex arrays of an ImageJ .roi, in image pixel coordinates.

    Only the polygon-like types carry an explicit vertex list; rectangles and
    ovals are rebuilt from their bounding box so no ROI is silently skipped.
    """
    if buf[:4] != b"Iout":
        raise ValueError("not an ImageJ ROI")
    roi_type = buf[6]
    top, left, bottom, right, n = struct.unpack(">hhhhh", buf[8:18])
    if roi_type in (0, 4, 5, 7, 8) and n > 0:          # polygon/freehand/traced
        xs = np.array(struct.unpack(">%dh" % n, buf[64:64 + 2 * n]), float) + left
        ys = np.array(struct.unpack(">%dh" % n, buf[64 + 2 * n:64 + 4 * n]), float) + top
        return xs, ys
    if roi_type == 1:                                   # rectangle
        return (np.array([left, right, right, left], float),
                np.array([top, top, bottom, bottom], float))
    if roi_type == 2:                                   # oval
        t = np.linspace(0, 2 * np.pi, 72, endpoint=False)
        cx, cy = (left + right) / 2.0, (top + bottom) / 2.0
        return (cx + (right - left) / 2.0 * np.cos(t),
                cy + (bottom - top) / 2.0 * np.sin(t))
    raise ValueError("unsupported ROI type %d" % roi_type)


def polygon_area_px(xs, ys):
    """Shoelace area. ImageJ traced/freehand ROIs store pixel-CORNER
    coordinates, so for those this is exactly the enclosed pixel count."""
    return 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(ys, np.roll(xs, -1)))


def nucleus_areas_px(zip_path):
    """{Cell_ID -> nucleus area in pixels} from one <Series>_v14_ROIs.zip."""
    out = {}
    with zipfile.ZipFile(zip_path) as z:
        for nm in z.namelist():
            m = re.search(r"Cell[_ ]?(\d+)", os.path.basename(nm))
            if not m:
                continue
            try:
                xs, ys = read_roi(z.read(nm))
            except (ValueError, struct.error) as e:
                print("    unreadable ROI %s in %s (%s)"
                      % (nm, os.path.basename(zip_path), e))
                continue
            out[int(m.group(1))] = polygon_area_px(xs, ys)
    return out


def tiff_info(tif_path):
    """(pixel size in microns, z-slices in the stack); either may be None.

    The z count is not used in any calculation. It is recorded because it is the
    confound described in the header: apparent 2D areas shrink as the stack gets
    deeper, so a reader can see at a glance which timepoints are comparable.
    """
    try:
        with Image.open(tif_path) as im:
            tags = im.tag_v2
            desc = tags.get(270, "") or ""
            meta = dict(l.split("=", 1) for l in desc.split("\n") if "=" in l)
            try:
                n_z = int(meta.get("slices", ""))
            except ValueError:
                n_z = None
            if "unit=micron" not in desc and "unit=um" not in desc:
                return None, n_z
            xres = tags.get(282)
            if xres is None:
                return None, n_z
            val = (float(xres[0]) / float(xres[1])
                   if isinstance(xres, tuple) else float(xres))
            return (1.0 / val if val else None), n_z
    except (OSError, ValueError, ZeroDivisionError):
        return None, None


# ------------------------------- parsing ------------------------------------
def cell_line_of(fname):
    for label, _, pat in sorted(CELL_LINES, key=lambda c: -len(c[2])):
        if pat.upper() in fname.upper():
            return label
    return None


def series_of(fname):
    m = SERIES_RE.search(fname)
    return int(m.group(1)) if m else -1


def time_of(fname):
    m = TIME_RE.search(fname)
    return int(m.group(1)) if m else None


def is_data(row, ci):
    try:
        int(float(row[ci]))
        return True
    except (ValueError, IndexError):
        return False


def read_image(path):
    """Per-cell raw values for one image CSV, or None if it has no usable data.

    Rows are recognised by an integer Cell_ID -- never by the File column ending
    in .tif, which would drop .czi datasets entirely.
    """
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return None
    h = rows[0]
    try:
        ci, pi, si, ai = (h.index("Cell_ID"), h.index("Punctae_Count"),
                          h.index("Slice"), h.index("Cyto_Area_px2"))
    except ValueError:
        return None
    data = [r for r in rows[1:] if is_data(r, ci)]
    if not data:
        return None

    punct, slices, areas = {}, {}, {}
    for r in data:
        cid = int(float(r[ci]))
        if cid == 0:                       # background
            continue
        punct[cid] = punct.get(cid, 0) + int(float(r[pi]))
        slices.setdefault(cid, set()).add(r[si].strip())
        try:
            areas.setdefault(cid, []).append(float(r[ai]))
        except ValueError:
            pass
    if not punct:
        return None
    return {"image": data[0][0], "punct": punct, "slices": slices, "areas": areas}


def sig_mark(p):
    if p != p:
        return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ------------------------------- analysis -----------------------------------
def collect(folder):
    roi_dir = os.path.normpath(os.path.join(folder, ROI_DIR_REL))
    tif_dir = os.path.normpath(os.path.join(folder, TIFF_DIR_REL))

    per_cell, per_image, provenance = [], [], []
    for path in sorted(glob.glob(os.path.join(folder, "*.csv"))):
        fname = os.path.basename(path)
        line, t, series = cell_line_of(fname), time_of(fname), series_of(fname)
        if line is None or t is None:
            print("  skipping (no cell line / timepoint in name): %s" % fname)
            continue
        img = read_image(path)
        if img is None:
            print("  skipping (no parseable data rows): %s" % fname)
            continue

        stem = fname[:-4]
        base = stem[:-4] if stem.endswith("_v14") else stem
        zip_path = os.path.join(roi_dir, base + "_v14_ROIs.zip")
        tif_path = os.path.join(tif_dir, base + ".tif")

        nuc_px = nucleus_areas_px(zip_path) if os.path.exists(zip_path) else {}
        if not nuc_px:
            print("  no ROI zip for %s -- nucleus size will be blank" % fname)
        px, n_z = tiff_info(tif_path) if os.path.exists(tif_path) else (None, None)
        px_src = "TIFF"
        if px is None:
            px, px_src = PIXEL_UM_FALLBACK, "fallback"
        px2 = px * px

        ana_date, acq_date = parse_dates(fname)
        img_slices = sorted({s for ss in img["slices"].values() for s in ss},
                            key=slice_sort)

        recs = []
        for cid in sorted(img["punct"]):
            sl = sorted(img["slices"][cid], key=slice_sort)
            areas = img["areas"].get(cid, [])
            n_um2 = nuc_px[cid] * px2 if cid in nuc_px else None
            # one number per cell, per CELL_AREA_STAT
            c_max = float(np.max(areas)) if areas else None
            c_avg = float(np.mean(areas)) if areas else None
            c_um2 = c_max if CELL_AREA_STAT == "max" else c_avg
            both = n_um2 is not None and c_um2 is not None
            recs.append({
                "line": line, "time": t, "series": series,
                "image": img["image"], "csv": fname,
                "acq_date": acq_date, "ana_date": ana_date,
                "cell_id": cid, "slices": ",".join(sl), "n_slices": len(sl),
                "punctae": img["punct"][cid],
                "sg_pos": "Yes" if img["punct"][cid] > 0 else "No",
                "cell_area": c_um2,
                "cell_area_max": c_max, "cell_area_mean": c_avg,
                "n_area_slices": len(areas),
                "nuc_area": n_um2,
                "nc_ratio": (n_um2 / c_um2) if (both and c_um2) else None,
                "cyto_only": (c_um2 - n_um2) if both else None,
                "pixel_um": round(px, 5),
            })
        per_cell.extend(recs)

        matched = sum(1 for r in recs if r["nuc_area"] is not None)
        provenance.append({
            "csv": fname, "image": img["image"], "line": line, "time": t,
            "series": series, "acq_date": acq_date, "ana_date": ana_date,
            "slices": ",".join(img_slices), "n_slices": len(img_slices),
            "roi_zip": (os.path.basename(zip_path)
                        if os.path.exists(zip_path) else "MISSING"),
            "n_rois": len(nuc_px), "n_cells": len(recs), "n_matched": matched,
            "pixel_um": round(px, 5), "pixel_src": px_src, "n_z": n_z,
        })

        # one pair of means per image: SG+ cells and SG- cells averaged separately
        row = {"csv": fname, "image": img["image"], "line": line, "time": t,
               "series": series, "acq_date": acq_date, "ana_date": ana_date,
               "slices": ",".join(img_slices), "n_slices": len(img_slices),
               "n_cells": len(recs), "n_z": n_z}
        for cls, want in (("pos", "Yes"), ("neg", "No")):
            sub = [r for r in recs if r["sg_pos"] == want]
            row["n_" + cls] = len(sub)
            for key in ("cell_area", "nuc_area", "nc_ratio"):
                v = [r[key] for r in sub if r[key] is not None]
                row["%s_%s" % (key, cls)] = float(np.mean(v)) if v else None
        # SG+ / SG- within the same stack: immune to the z-depth confound, so
        # this is the quantity that can be compared across timepoints.
        for key in ("cell_area", "nuc_area", "nc_ratio"):
            p, q = row[key + "_pos"], row[key + "_neg"]
            row["ratio_" + key] = (p / q) if (p and q) else None
        per_image.append(row)

    return per_cell, per_image, provenance


def usable(d, cls, key):
    """Does this image contribute a usable mean for one SG class?"""
    return (d["n_" + cls] >= MIN_CELLS_PER_CLASS
            and d["%s_%s" % (key, cls)] is not None)


def paired_tests(per_image, key):
    """Per (line, time): paired t-test and Wilcoxon of SG+ vs SG- image means."""
    out = {}
    for line, _, _ in CELL_LINES:
        for t in TIMEPOINTS:
            pos, neg = [], []
            for d in per_image:
                if d["line"] != line or d["time"] != t:
                    continue
                if usable(d, "pos", key) and usable(d, "neg", key):
                    pos.append(d[key + "_pos"])
                    neg.append(d[key + "_neg"])
            n = len(pos)
            res = {"n_pairs": n,
                   "mean_pos": float(np.mean(pos)) if n else float("nan"),
                   "mean_neg": float(np.mean(neg)) if n else float("nan"),
                   "t": float("nan"), "p": float("nan"),
                   "w_p": float("nan"), "d": float("nan")}
            if n >= 2:
                diff = np.array(pos) - np.array(neg)
                res["t"], res["p"] = stats.ttest_rel(pos, neg)
                sd = diff.std(ddof=1)
                res["d"] = diff.mean() / sd if sd else float("nan")
                if n >= 5 and np.any(diff != 0):
                    try:
                        res["w_p"] = stats.wilcoxon(pos, neg).pvalue
                    except ValueError:
                        pass
            out[(line, t)] = res
    return out


def image_ratios(per_image, line, t, key):
    """Per-image SG+/SG- ratios for one group, using only usable image pairs."""
    return [d["ratio_" + key] for d in per_image
            if d["line"] == line and d["time"] == t
            and usable(d, "pos", key) and usable(d, "neg", key)
            and d["ratio_" + key]]


def ratio_tests(per_image, key):
    """Per (line, time): is the SG+/SG- size ratio different from 1?

    Tested on log2 ratios with a one-sample t-test, so that a 2x increase and a
    2x decrease weigh the same; the reported ratio is the geometric mean, which
    is what that test is actually centring.
    """
    out = {}
    for line, _, _ in CELL_LINES:
        for t in TIMEPOINTS:
            r = np.array(image_ratios(per_image, line, t, key), float)
            res = {"n": len(r), "ratio": float("nan"), "sd": float("nan"),
                   "t": float("nan"), "p": float("nan")}
            if len(r):
                lg = np.log2(r)
                res["ratio"] = float(2 ** lg.mean())
                res["sd"] = float(r.std(ddof=1)) if len(r) > 1 else 0.0
                if len(r) >= 2 and lg.std(ddof=1) > 0:
                    res["t"], res["p"] = stats.ttest_1samp(lg, 0.0)
            out[(line, t)] = res
    return out


def cell_values(per_cell, line, t, key, cls):
    """Every scored cell of one class at one cell line x timepoint.

    No image-level gate at all: a cell counts wherever it was measured, which is
    the point of the cell-level bars.
    """
    want = "Yes" if cls == "pos" else "No"
    return [r[key] for r in per_cell
            if r["line"] == line and r["time"] == t
            and r["sg_pos"] == want and r[key] is not None]


def image_means(per_image, line, t, key, cls):
    """Per-image means of one class: every image with any cells of that class."""
    return [d["%s_%s" % (key, cls)] for d in per_image
            if d["line"] == line and d["time"] == t
            and d["n_" + cls] >= 1
            and d["%s_%s" % (key, cls)] is not None]


def allcell_tests(per_cell, per_image, key):
    """Unpaired SG+ vs SG- using every cell, per cell line x timepoint.

    Welch's t-test (unequal variances, unequal n) is the headline because it
    matches how the bars are built. The Mann-Whitney U and the unpaired
    image-level test are carried alongside as sanity checks: if the cell-level
    p-value is small only because cells within an image repeat each other, the
    image-level one will not agree.
    """
    out = {}
    for line, _, _ in CELL_LINES:
        for t in TIMEPOINTS:
            pos = np.array(cell_values(per_cell, line, t, key, "pos"), float)
            neg = np.array(cell_values(per_cell, line, t, key, "neg"), float)
            res = {"n_pos": len(pos), "n_neg": len(neg),
                   "mean_pos": float(pos.mean()) if len(pos) else float("nan"),
                   "mean_neg": float(neg.mean()) if len(neg) else float("nan"),
                   "sd_pos": float(pos.std(ddof=1)) if len(pos) > 1 else float("nan"),
                   "sd_neg": float(neg.std(ddof=1)) if len(neg) > 1 else float("nan"),
                   "t": float("nan"), "p": float("nan"), "mw_p": float("nan"),
                   "d": float("nan"),
                   "img_n_pos": 0, "img_n_neg": 0, "img_p": float("nan")}
            if len(pos) >= 2 and len(neg) >= 2:
                res["t"], res["p"] = stats.ttest_ind(pos, neg, equal_var=False)
                sp = math.sqrt((pos.var(ddof=1) + neg.var(ddof=1)) / 2.0)
                res["d"] = float((pos.mean() - neg.mean()) / sp) if sp else float("nan")
            if len(pos) >= 3 and len(neg) >= 3:
                res["mw_p"] = float(stats.mannwhitneyu(
                    pos, neg, alternative="two-sided").pvalue)
            ip = np.array(image_means(per_image, line, t, key, "pos"), float)
            ineg = np.array(image_means(per_image, line, t, key, "neg"), float)
            res["img_n_pos"], res["img_n_neg"] = len(ip), len(ineg)
            if len(ip) >= 2 and len(ineg) >= 2:
                res["img_p"] = float(stats.ttest_ind(ip, ineg, equal_var=False).pvalue)
            out[(line, t)] = res
    return out


def descriptive_tests(per_cell, key):
    """Per (line, time): Mann-Whitney U on pooled cells. Pseudo-replicated --
    reported for description only, never as the headline test."""
    out = {}
    for line, _, _ in CELL_LINES:
        for t in TIMEPOINTS:
            pos = [r[key] for r in per_cell
                   if r["line"] == line and r["time"] == t
                   and r["sg_pos"] == "Yes" and r[key] is not None]
            neg = [r[key] for r in per_cell
                   if r["line"] == line and r["time"] == t
                   and r["sg_pos"] == "No" and r[key] is not None]
            u, p = float("nan"), float("nan")
            if len(pos) >= 3 and len(neg) >= 3:
                u, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
            out[(line, t)] = {
                "n_pos": len(pos), "n_neg": len(neg),
                "median_pos": float(np.median(pos)) if pos else float("nan"),
                "median_neg": float(np.median(neg)) if neg else float("nan"),
                "U": u, "p": p}
    return out


# -------------------------------- figures -----------------------------------
ABS_SUBJECT = {
    "cell_area": "Cell area, taken as the largest cross-section across each "
                 "cell's slices",
    "nuc_area":  "Nucleus area, from each cell's single traced DAPI outline",
    "nc_ratio":  "Nucleus-to-cell area ratio per cell",
}


def abs_caption(key):
    return (
        "%s, in SG− and SG+ cells over the NaAsO2 time course, one panel per "
        "cell line. Mean ± SD across images, with the individual image means "
        "jittered behind. Stars are the paired t-test of SG+ against SG− image "
        "means (n = images with at least %d cells of both classes). Compare SG+ "
        "with SG− within a timepoint only: z-stack depth differs between "
        "timepoints and scales the apparent area, so the left-to-right trend is "
        "a property of the acquisition as much as of the biology. For the time "
        "course use the SG+/SG− ratio figure."
        % (ABS_SUBJECT[key], MIN_CELLS_PER_CLASS))

RATIO_CAPTION = (
    "SG+ / SG− size ratio within each image, one panel per cell line. Each image "
    "contributes one ratio, so both numbers come from the same z stack and the "
    "z-depth differences between timepoints cancel — this is the figure that can "
    "honestly be read across time. Points are the geometric mean across images "
    "with SD; the dashed line at 1 is no difference. Ratio > 1 means SG+ cells "
    "are larger. Tested as a one-sample t-test of log2(ratio) against 0.")


def bar_caption():
    if BAR_UNIT == "cell":
        return (
            "Cell area and nucleus area in SG− and SG+ cells at each timepoint, "
            "one panel per cell line. Every scored cell contributes: bars are "
            "mean ± SD over ALL cells, grey numbers are n cells per bar "
            "(SG−/SG+), and the two classes are deliberately unequal. Dots are "
            "the per-image means. Cell size is the LARGEST Cyto_Area across that "
            "cell's slices, one number per cell; nucleus size is its single "
            "traced DAPI outline. Stars are an unpaired Welch t-test over all "
            "cells — cells within an image share focus, staining and z depth and "
            "are not independent, so this p-value is optimistic; the paired "
            "by-image test is in the Statistics sheet and the by-image unpaired "
            "test in Group_stats_allcells. Compare the two columns within a "
            "timepoint, not across timepoints.")
    return (
        "Cell area and nucleus area in SG− and SG+ cells at each timepoint, one "
        "panel per cell line. Bars are mean ± SD over images of each class's "
        "per-image mean; dots are the individual images; grey numbers are n "
        "images per bar (SG−/SG+). Stars are the paired t-test on image means, "
        "so they appear only where at least %d cells of BOTH classes share an "
        "image. Compare the two columns within a timepoint, not across "
        "timepoints — z-stack depth differs between timepoints and scales the "
        "apparent area." % MIN_CELLS_PER_CLASS)


def figure_captions():
    """[(figure stem, caption)] in the order the figures are written."""
    out = [(stem, abs_caption(k)) for k, _y, _s, _l, stem in METRICS]
    out.append((RATIO_STEM, RATIO_CAPTION))
    out.append((BAR_STEM, bar_caption()))
    return out


def write_captions(out_dir, folder):
    """Captions for every figure, as a plain text file beside them."""
    path = os.path.join(out_dir, CAPTIONS_TXT)
    with open(path, "w", encoding="utf-8") as f:
        f.write("Figure captions - SG+ vs SG- cell and nucleus size\n")
        f.write("Generated by sg_size_compare.py\n")
        f.write("Source: %s\n" % folder)
        f.write("Cell size statistic: %s of each cell's per-slice Cyto_Area\n"
                % CELL_AREA_STAT)
        f.write("=" * 78 + "\n")
        for stem, text in figure_captions():
            f.write("\n%s\n%s\n" % (stem, "-" * len(stem)))
            f.write("\n".join(textwrap.wrap(text, 78)) + "\n")
    return path


def footnote(fig, text, wrap=115):
    """Caption under the axes, hard-wrapped.

    Wrapping is not cosmetic: the figures are saved with bbox_inches="tight",
    which grows the canvas to fit the widest artist, so one long unwrapped line
    silently stretches the whole figure and squashes the panels.
    """
    if not text or not SHOW_FOOTNOTE:
        return
    fig.text(0.5, 0.005, "\n".join(textwrap.wrap(text, wrap)),
             ha="center", va="bottom", fontsize=FONT_FOOTNOTE, color="0.30")


def apply_style():
    plt.rcParams.update({
        "font.family": FONT_FAMILY, "font.size": FONT_BASE,
        "axes.labelsize": FONT_AXLABEL, "axes.titlesize": FONT_TITLE,
        "xtick.labelsize": FONT_TICKS, "ytick.labelsize": FONT_TICKS,
        "legend.fontsize": FONT_LEGEND, "svg.fonttype": "none",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def make_figure(per_image, tests, key, ylabel, ylim, out_dir, stem, show):
    apply_style()
    fig, axes = plt.subplots(1, len(CELL_LINES), figsize=FIGSIZE, sharey=True)
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    handles = None

    for ax, (line, _, _) in zip(axes, CELL_LINES):
        for cls, lbl, col, dx in (("neg", SG_NEG_LABEL, SG_NEG_COLOR, -0.7),
                                  ("pos", SG_POS_LABEL, SG_POS_COLOR, +0.7)):
            xs, ms, es = [], [], []
            for t in TIMEPOINTS:
                v = [d["%s_%s" % (key, cls)] for d in per_image
                     if d["line"] == line and d["time"] == t and usable(d, cls, key)]
                if not v:
                    continue
                v = np.array(v, float)
                xs.append(t + dx)
                ms.append(v.mean())
                sd = v.std(ddof=1) if len(v) > 1 else 0.0
                es.append(sd / math.sqrt(len(v)) if ERRORBAR == "SEM" else sd)
                if SHOW_POINTS:
                    ax.scatter(t + dx + rng.uniform(-POINT_JITTER, POINT_JITTER, len(v)),
                               v, s=POINT_SIZE, color=col, alpha=POINT_ALPHA,
                               edgecolors="none", zorder=2)
            ax.errorbar(xs, ms, yerr=es, label=lbl, color=col, lw=LINEWIDTH,
                        marker="o", ms=MARKERSIZE, capsize=CAPSIZE, zorder=3)

        ax.set_title(line)
        ax.set_xlabel(XLABEL)
        if XLIM:
            ax.set_xlim(*XLIM)
        if XTICKS:
            ax.set_xticks(XTICKS)
        if ylim:
            ax.set_ylim(*ylim)
        if HIDE_TOP_RIGHT_SPINES:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_linewidth(SPINE_WIDTH)
        if handles is None:
            handles = ax.get_legend_handles_labels()

    axes[0].set_ylabel(ylabel)

    if STAR_ON_PLOT:
        # shared y-axis: grow the room once, then annotate every panel
        lo, hi = axes[0].get_ylim()
        axes[0].set_ylim(lo, hi + (hi - lo) * 0.12)
        lo, hi = axes[0].get_ylim()
        for ax, (line, _, _) in zip(axes, CELL_LINES):
            for t in TIMEPOINTS:
                r = tests.get((line, t))
                if not r or r["n_pairs"] < 2:
                    continue
                ax.text(t, hi - (hi - lo) * 0.05, sig_mark(r["p"]),
                        ha="center", va="center", fontsize=FONT_BASE)

    fig.legend(*handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=(0, 0.06 if SHOW_FOOTNOTE else 0, 1, 0.93))
    footnote(fig, abs_caption(key))

    paths = []
    for ext in FORMATS:
        p = os.path.join(out_dir, "%s.%s" % (stem, ext))
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        paths.append(p)
    if PPT_SVG:
        # Same figure in a font Windows actually has: DejaVu Sans ships with
        # matplotlib and is not installed, so PowerPoint would substitute it and
        # shift every label.
        plt.rcParams["font.family"] = PPT_FONT
        for ax in axes:
            for txt in ([ax.title, ax.xaxis.label, ax.yaxis.label]
                        + list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                        + list(ax.texts)):
                txt.set_fontfamily(PPT_FONT)
        for lg in fig.legends:
            for txt in lg.get_texts():
                txt.set_fontfamily(PPT_FONT)
        for txt in fig.texts:
            txt.set_fontfamily(PPT_FONT)
        p = os.path.join(out_dir, "%s%s.svg" % (stem, PPT_SUFFIX))
        fig.savefig(p, bbox_inches="tight")
        paths.append(p)
    if show:
        plt.show()
    plt.close(fig)
    return paths


def make_ratio_figure(per_image, rtests, out_dir, show):
    """SG+ / SG- size ratio vs time -- the z-independent version of the figures
    above. A ratio of 1 (dashed line) means SG+ and SG- cells are the same size."""
    apply_style()
    fig, axes = plt.subplots(1, len(CELL_LINES), figsize=FIGSIZE, sharey=True)
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    handles = None

    for ax, (line, _, _) in zip(axes, CELL_LINES):
        ax.axhline(1.0, color="0.35", lw=1.2, ls="--", zorder=1)
        for key, lbl, col in RATIO_SERIES:
            dx = RATIO_DX[key]
            xs, ms, es = [], [], []
            for t in TIMEPOINTS:
                v = np.array(image_ratios(per_image, line, t, key), float)
                if not len(v):
                    continue
                xs.append(t + dx)
                ms.append(rtests[key][(line, t)]["ratio"])
                sd = v.std(ddof=1) if len(v) > 1 else 0.0
                es.append(sd / math.sqrt(len(v)) if ERRORBAR == "SEM" else sd)
                if SHOW_POINTS:
                    ax.scatter(t + dx + rng.uniform(-POINT_JITTER, POINT_JITTER, len(v)),
                               v, s=POINT_SIZE, color=col, alpha=POINT_ALPHA,
                               edgecolors="none", zorder=2)
            ax.errorbar(xs, ms, yerr=es, label=lbl, color=col, lw=LINEWIDTH,
                        marker="o", ms=MARKERSIZE, capsize=CAPSIZE, zorder=3)

        ax.set_title(line)
        ax.set_xlabel(XLABEL)
        if XLIM:
            ax.set_xlim(*XLIM)
        if XTICKS:
            ax.set_xticks(XTICKS)
        if RATIO_YLIM:
            ax.set_ylim(*RATIO_YLIM)
        if HIDE_TOP_RIGHT_SPINES:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_linewidth(SPINE_WIDTH)
        if handles is None:
            handles = ax.get_legend_handles_labels()

    axes[0].set_ylabel("SG+ / SG− (per image)")
    fig.legend(*handles, loc="upper center", ncol=len(RATIO_SERIES), frameon=False,
               bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=(0, 0.06 if SHOW_FOOTNOTE else 0, 1, 0.93))
    footnote(fig, RATIO_CAPTION)

    paths = []
    for ext in FORMATS:
        p = os.path.join(out_dir, "%s.%s" % (RATIO_STEM, ext))
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        paths.append(p)
    if PPT_SVG:
        plt.rcParams["font.family"] = PPT_FONT
        for ax in axes:
            for txt in ([ax.title, ax.xaxis.label, ax.yaxis.label]
                        + list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                        + list(ax.texts)):
                txt.set_fontfamily(PPT_FONT)
        for lg in fig.legends:
            for txt in lg.get_texts():
                txt.set_fontfamily(PPT_FONT)
        for txt in fig.texts:
            txt.set_fontfamily(PPT_FONT)
        p = os.path.join(out_dir, "%s%s.svg" % (RATIO_STEM, PPT_SUFFIX))
        fig.savefig(p, bbox_inches="tight")
        paths.append(p)
    if show:
        plt.show()
    plt.close(fig)
    return paths


def make_bars_figure(per_cell, per_image, tests, atests, out_dir, show):
    """SG- and SG+ as two columns at each of the four timepoints.

    Rows are the two sizes asked for (cell, nucleus), columns the three cell
    lines. With BAR_UNIT == "cell" every scored cell contributes: the bar is the
    mean +/- SD over cells, the grey numbers are cell counts, and the stars are
    the unpaired Welch t-test over all cells. Per-image means are still drawn as
    dots, so an effect carried by a single image stays visible.
    """
    apply_style()
    by_cell = BAR_UNIT == "cell"
    fig, axes = plt.subplots(len(BAR_METRICS), len(CELL_LINES),
                             figsize=BAR_FIGSIZE, sharey="row", squeeze=False)
    rng = np.random.default_rng(0)
    xpos = np.arange(len(TIMEPOINTS), dtype=float)
    handles = None

    for r, (key, ylabel, ylim) in enumerate(BAR_METRICS):
        for c, (line, _, _) in enumerate(CELL_LINES):
            ax = axes[r][c]
            counts = []
            for cls, lbl, col, dx in (("neg", SG_NEG_LABEL, SG_NEG_COLOR, -BAR_WIDTH / 2),
                                      ("pos", SG_POS_LABEL, SG_POS_COLOR, +BAR_WIDTH / 2)):
                bx, bm, be, ns = [], [], [], []
                for i, t in enumerate(TIMEPOINTS):
                    cells = cell_values(per_cell, line, t, key, cls)
                    imgs = image_means(per_image, line, t, key, cls)
                    v = np.array(cells if by_cell else imgs, float)
                    ns.append(len(v))
                    if not len(v):                 # no cells of this class here
                        continue
                    sd = v.std(ddof=1) if len(v) > 1 else 0.0
                    if ERRORBAR == "SEM" and len(v) > 1:
                        sd /= math.sqrt(len(v))
                    bx.append(xpos[i] + dx); bm.append(v.mean()); be.append(sd)
                    if SHOW_POINTS:
                        # dots are always the per-image means: with thousands of
                        # cells a per-cell swarm would just fill the panel
                        d = np.array(imgs, float)
                        ax.scatter(xpos[i] + dx + rng.uniform(-BAR_WIDTH / 3,
                                                              BAR_WIDTH / 3, len(d)),
                                   d, s=POINT_SIZE, color="black", alpha=0.45,
                                   edgecolors="none", zorder=4)
                ax.bar(bx, bm, BAR_WIDTH, yerr=be, capsize=CAPSIZE, label=lbl,
                       color=col, alpha=BAR_ALPHA, edgecolor="black",
                       linewidth=0.8, zorder=2)
                counts.append(ns)

            ax.set_xticks(xpos)
            ax.set_xticklabels(["%d" % t for t in TIMEPOINTS])
            ax.set_xlim(-0.6, len(TIMEPOINTS) - 0.4)
            if ylim:
                ax.set_ylim(*ylim)
            if HIDE_TOP_RIGHT_SPINES:
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
            for s in ("left", "bottom"):
                ax.spines[s].set_linewidth(SPINE_WIDTH)
            if r == 0:
                ax.set_title(line)
            if r == len(BAR_METRICS) - 1:
                ax.set_xlabel(XLABEL, labelpad=26)   # clear of the n labels
            if c == 0:
                ax.set_ylabel(ylabel)
            if handles is None:
                handles = ax.get_legend_handles_labels()

            # n images per bar, under the axis
            for i in range(len(TIMEPOINTS)):
                ax.annotate("%d/%d" % (counts[0][i], counts[1][i]),
                            xy=(xpos[i], 0), xycoords=("data", "axes fraction"),
                            xytext=(0, -26), textcoords="offset points",
                            ha="center", va="top", fontsize=FONT_BASE - 3,
                            color="0.35", annotation_clip=False)

        # Headroom and stars are done once the whole row is drawn: the panels
        # share a y-axis, so expanding inside the panel loop would stretch the
        # row again for every column and leave the stars at stepped heights.
        row_ax = axes[r][0]
        lo, hi = row_ax.get_ylim()
        row_ax.set_ylim(lo, hi + (hi - lo) * 0.13)
        lo, hi = row_ax.get_ylim()
        for c, (line, _, _) in enumerate(CELL_LINES):
            for i, t in enumerate(TIMEPOINTS):
                if by_cell:
                    res = atests[key].get((line, t))
                    ok = res and res["n_pos"] >= 2 and res["n_neg"] >= 2
                else:
                    res = tests[key].get((line, t))
                    ok = res and res["n_pairs"] >= 2
                if not ok:
                    continue
                axes[r][c].text(xpos[i], hi - (hi - lo) * 0.05, sig_mark(res["p"]),
                                ha="center", va="center", fontsize=FONT_BASE)

    fig.legend(*handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.005))
    fig.tight_layout(rect=(0, 0.085 if SHOW_FOOTNOTE else 0.02, 1, 0.945))
    footnote(fig, bar_caption())

    paths = []
    for ext in FORMATS:
        p = os.path.join(out_dir, "%s.%s" % (BAR_STEM, ext))
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        paths.append(p)
    if PPT_SVG:
        plt.rcParams["font.family"] = PPT_FONT
        for row in axes:
            for ax in row:
                for txt in ([ax.title, ax.xaxis.label, ax.yaxis.label]
                            + list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                            + list(ax.texts)):
                    txt.set_fontfamily(PPT_FONT)
        for lg in fig.legends:
            for txt in lg.get_texts():
                txt.set_fontfamily(PPT_FONT)
        for txt in fig.texts:
            txt.set_fontfamily(PPT_FONT)
        p = os.path.join(out_dir, "%s%s.svg" % (BAR_STEM, PPT_SUFFIX))
        fig.savefig(p, bbox_inches="tight")
        paths.append(p)
    if show:
        plt.show()
    plt.close(fig)
    return paths


# -------------------------------- workbook ----------------------------------
def r4(v, nd=4):
    return round(v, nd) if isinstance(v, float) and v == v else v


def rounded(records):
    return [{k: r4(v) for k, v in d.items()} for d in records]


def write_workbook(out_path, per_cell, per_image, provenance, tests, desc,
                   rtests, atests, figures):
    wb = openpyxl.Workbook()

    ws = wb.active
    ws.title = "Per_cell"
    write_sheet(ws, [
        ("Cell_line", "line", 15), ("Time_min", "time", 10),
        ("Series", "series", 8), ("Image", "image", 46),
        ("Source_CSV", "csv", 46), ("Acquisition_date", "acq_date", 17),
        ("Analysis_date", "ana_date", 15), ("Cell_ID", "cell_id", 9),
        ("Slices", "slices", 16), ("N_slices", "n_slices", 10),
        ("Punctae_per_cell", "punctae", 18), ("SG_positive", "sg_pos", 12),
        ("Cell_area_um2", "cell_area", 15),
        ("Cell_area_max_um2", "cell_area_max", 19),
        ("Cell_area_mean_um2", "cell_area_mean", 20),
        ("N_slices_with_area", "n_area_slices", 19),
        ("Nucleus_area_um2", "nuc_area", 18),
        ("NC_ratio", "nc_ratio", 11),
        ("Cytoplasm_only_um2", "cyto_only", 20),
        ("Pixel_um", "pixel_um", 10),
    ], rounded(per_cell))

    ws2 = wb.create_sheet("Per_image")
    write_sheet(ws2, [
        ("Cell_line", "line", 15), ("Time_min", "time", 10),
        ("Series", "series", 8), ("Source_CSV", "csv", 46),
        ("Acquisition_date", "acq_date", 17), ("Analysis_date", "ana_date", 15),
        ("Slices_used", "slices", 16), ("N_slices", "n_slices", 10),
        ("N_z_in_stack", "n_z", 14),
        ("Total_cells", "n_cells", 12),
        ("N_SGpos", "n_pos", 10), ("N_SGneg", "n_neg", 10),
        ("Cell_area_SGpos", "cell_area_pos", 17),
        ("Cell_area_SGneg", "cell_area_neg", 17),
        ("Nucleus_area_SGpos", "nuc_area_pos", 19),
        ("Nucleus_area_SGneg", "nuc_area_neg", 19),
        ("NC_ratio_SGpos", "nc_ratio_pos", 16),
        ("NC_ratio_SGneg", "nc_ratio_neg", 16),
        ("Ratio_cell_area", "ratio_cell_area", 17),
        ("Ratio_nucleus_area", "ratio_nuc_area", 20),
        ("Ratio_NC", "ratio_nc_ratio", 11),
    ], rounded(per_image))

    ws3 = wb.create_sheet("Group_means")
    recs = []
    for line, _, _ in CELL_LINES:
        for t in TIMEPOINTS:
            for cls, lbl in (("neg", "SG-"), ("pos", "SG+")):
                sub = [d for d in per_image
                       if d["line"] == line and d["time"] == t
                       and d["n_" + cls] >= MIN_CELLS_PER_CLASS]
                rec = {"line": line, "time": t, "sg": lbl, "n_images": len(sub),
                       "n_cells": sum(d["n_" + cls] for d in sub)}
                for key, nm in (("cell_area", "cell"), ("nuc_area", "nuc"),
                                ("nc_ratio", "nc")):
                    v = [d["%s_%s" % (key, cls)] for d in sub
                         if d["%s_%s" % (key, cls)] is not None]
                    rec[nm + "_mean"] = r4(float(np.mean(v))) if v else None
                    rec[nm + "_sd"] = r4(float(np.std(v, ddof=1))) if len(v) > 1 else None
                recs.append(rec)
    write_sheet(ws3, [
        ("Cell_line", "line", 15), ("Time_min", "time", 10),
        ("SG_class", "sg", 10), ("n_images", "n_images", 10),
        ("n_cells", "n_cells", 10),
        ("Cell_area_mean", "cell_mean", 16), ("Cell_area_SD", "cell_sd", 14),
        ("Nucleus_area_mean", "nuc_mean", 18), ("Nucleus_area_SD", "nuc_sd", 16),
        ("NC_ratio_mean", "nc_mean", 15), ("NC_ratio_SD", "nc_sd", 13),
    ], recs)

    ws4 = wb.create_sheet("Statistics")
    srecs = []
    for mkey, mlabel in (("cell_area", "Cell area (um2)"),
                         ("nuc_area", "Nucleus area (um2)"),
                         ("nc_ratio", "N/C ratio")):
        for line, _, _ in CELL_LINES:
            for t in TIMEPOINTS:
                r, d = tests[mkey][(line, t)], desc[mkey][(line, t)]
                srecs.append({
                    "metric": mlabel, "line": line, "time": t,
                    "n_pairs": r["n_pairs"],
                    "mean_pos": r4(r["mean_pos"]), "mean_neg": r4(r["mean_neg"]),
                    "delta": r4(r["mean_pos"] - r["mean_neg"]),
                    "t": r4(r["t"]), "p": r4(r["p"], 6), "sig": sig_mark(r["p"]),
                    "d": r4(r["d"], 3), "w_p": r4(r["w_p"], 6),
                    "cells_pos": d["n_pos"], "cells_neg": d["n_neg"],
                    "mw_p": r4(d["p"], 6),
                })
    write_sheet(ws4, [
        ("Metric", "metric", 20), ("Cell_line", "line", 15),
        ("Time_min", "time", 10), ("n_pairs_images", "n_pairs", 15),
        ("Mean_SGpos", "mean_pos", 13), ("Mean_SGneg", "mean_neg", 13),
        ("Difference", "delta", 12), ("t", "t", 10),
        ("p_paired_t", "p", 12), ("Significance", "sig", 13),
        ("Cohens_dz", "d", 12), ("p_Wilcoxon", "w_p", 13),
        ("n_cells_SGpos", "cells_pos", 15), ("n_cells_SGneg", "cells_neg", 15),
        ("p_MannWhitney_DESCRIPTIVE", "mw_p", 26),
    ], srecs)
    for note in (
        "",
        "PRIMARY TEST: paired t-test of SG+ vs SG- image MEANS. Statistical unit "
        "= image; an image contributes a pair only if it has >= %d cells in BOTH "
        "classes." % MIN_CELLS_PER_CLASS,
        "Significance is read from p_paired_t: *** p<0.001, ** p<0.01, * p<0.05, "
        "ns = not significant.",
        "p_MannWhitney_DESCRIPTIVE pools individual cells and treats them as "
        "independent. It is pseudo-replicated and is NOT the headline test.",
        "Blank t / p = fewer than 2 usable image pairs at that timepoint.",
        "Nucleus area is a single-plane ROI per cell (the Fiji macro saves one "
        "outline per cell); cell area is the mean over that cell's subset slices.",
        "DO NOT compare Mean_SGpos / Mean_SGneg ACROSS timepoints: apparent area "
        "tracks the z-stack depth, which differs per .lif (see N_z_in_stack on "
        "Per_image). Use the Ratio_stats sheet for the time course.",
    ):
        ws4.append([note])

    ws7 = wb.create_sheet("Ratio_stats")
    rrecs = []
    for mkey, mlabel in (("cell_area", "Cell area"),
                         ("nuc_area", "Nucleus area"),
                         ("nc_ratio", "N/C ratio")):
        for line, _, _ in CELL_LINES:
            for t in TIMEPOINTS:
                r = rtests[mkey][(line, t)]
                rrecs.append({
                    "metric": mlabel, "line": line, "time": t, "n": r["n"],
                    "ratio": r4(r["ratio"]), "sd": r4(r["sd"]),
                    "pct": r4((r["ratio"] - 1) * 100, 2) if r["ratio"] == r["ratio"] else None,
                    "t": r4(r["t"]), "p": r4(r["p"], 6), "sig": sig_mark(r["p"]),
                })
    write_sheet(ws7, [
        ("Metric", "metric", 18), ("Cell_line", "line", 15),
        ("Time_min", "time", 10), ("n_images", "n", 11),
        ("Ratio_SGpos_over_SGneg", "ratio", 24), ("SD_of_ratio", "sd", 13),
        ("Percent_difference", "pct", 19), ("t", "t", 10),
        ("p_vs_ratio_1", "p", 14), ("Significance", "sig", 13),
    ], rrecs)
    for note in (
        "",
        "Each image gives ONE ratio: mean of its SG+ cells / mean of its SG- cells. "
        "Both come from the same z stack, so this is unaffected by the z-depth "
        "confound that makes the absolute areas incomparable across timepoints.",
        "Ratio > 1 means SG+ cells are LARGER than SG- cells in the same image.",
        "Test: one-sample t-test of log2(ratio) against 0, so a 2x increase and a "
        "2x decrease carry equal weight. Ratio_SGpos_over_SGneg is the geometric "
        "mean across images, which is the value that test centres.",
        "Same image-pair inclusion rule as the Statistics sheet: >= %d cells in "
        "BOTH classes." % MIN_CELLS_PER_CLASS,
    ):
        ws7.append([note])

    ws8 = wb.create_sheet("Group_stats_allcells")
    arecs = []
    for mkey, mlabel in (("cell_area", "Cell area (um2)"),
                         ("nuc_area", "Nucleus area (um2)"),
                         ("nc_ratio", "N/C ratio")):
        for line, _, _ in CELL_LINES:
            for t in TIMEPOINTS:
                a = atests[mkey][(line, t)]
                arecs.append({
                    "metric": mlabel, "line": line, "time": t,
                    "n_neg": a["n_neg"], "n_pos": a["n_pos"],
                    "mean_neg": r4(a["mean_neg"]), "sd_neg": r4(a["sd_neg"]),
                    "mean_pos": r4(a["mean_pos"]), "sd_pos": r4(a["sd_pos"]),
                    "delta": r4(a["mean_pos"] - a["mean_neg"]),
                    "pct": r4((a["mean_pos"] / a["mean_neg"] - 1) * 100, 2)
                           if a["mean_neg"] else None,
                    "t": r4(a["t"]), "p": r4(a["p"], 6), "sig": sig_mark(a["p"]),
                    "d": r4(a["d"], 3), "mw_p": r4(a["mw_p"], 6),
                    "img_n": "%d/%d" % (a["img_n_neg"], a["img_n_pos"]),
                    "img_p": r4(a["img_p"], 6),
                })
    write_sheet(ws8, [
        ("Metric", "metric", 20), ("Cell_line", "line", 15),
        ("Time_min", "time", 10),
        ("n_cells_SGneg", "n_neg", 15), ("n_cells_SGpos", "n_pos", 15),
        ("Mean_SGneg", "mean_neg", 13), ("SD_SGneg", "sd_neg", 11),
        ("Mean_SGpos", "mean_pos", 13), ("SD_SGpos", "sd_pos", 11),
        ("Difference", "delta", 12), ("Percent_difference", "pct", 19),
        ("t_Welch", "t", 11), ("p_Welch", "p", 12), ("Significance", "sig", 13),
        ("Cohens_d", "d", 11), ("p_MannWhitney", "mw_p", 15),
        ("n_images_neg/pos", "img_n", 18), ("p_unpaired_by_image", "img_p", 21),
    ], arecs)
    for note in (
        "",
        "EVERY scored cell is used: no image-level gate, no pairing, unequal n "
        "between the two classes is expected and fine. This sheet is what the "
        "bar figure shows.",
        "Cell area per cell = the LARGEST Cyto_Area across that cell's slices "
        "(CELL_AREA_STAT = '%s'), i.e. the slice nearest the cell's equator. "
        "Nucleus area per cell = its single traced DAPI ROI." % CELL_AREA_STAT,
        "p_Welch is an unpaired t-test with unequal variances over all cells. "
        "Cells inside one image share focus, staining and z depth, so they are "
        "not independent and p_Welch is optimistic.",
        "p_unpaired_by_image repeats the test on per-image means (n = images "
        "with any cells of that class). Where it disagrees with p_Welch, the "
        "cell-level result is being carried by within-image correlation.",
        "The Statistics sheet holds the stricter paired-by-image test, which "
        "discards images lacking %d cells of both classes." % MIN_CELLS_PER_CLASS,
        "Still compare only WITHIN a timepoint: z-stack depth differs between "
        "timepoints and scales the apparent area.",
    ):
        ws8.append([note])

    ws5 = wb.create_sheet("Provenance")
    write_sheet(ws5, [
        ("Source_CSV", "csv", 46), ("Image", "image", 46),
        ("Cell_line", "line", 15), ("Time_min", "time", 10),
        ("Series", "series", 8), ("Acquisition_date", "acq_date", 17),
        ("Analysis_date", "ana_date", 15),
        ("Slices_used", "slices", 16), ("N_slices", "n_slices", 10),
        ("ROI_zip", "roi_zip", 46), ("N_ROIs_in_zip", "n_rois", 15),
        ("N_cells_counted", "n_cells", 16),
        ("N_cells_with_nucleus", "n_matched", 21),
        ("Pixel_um", "pixel_um", 10), ("Pixel_source", "pixel_src", 13),
        ("N_z_in_stack", "n_z", 14),
    ], provenance)

    ws6 = wb.create_sheet("Figures")
    row = 1
    for png in figures:
        ws6.add_image(XLImage(png), "A%d" % row)
        row += 32

    wb.save(out_path)


# ---------------------------------- main ------------------------------------
def main(folder, show=False):
    out_dir = OUT_DIR or folder
    print("CSV folder :", folder)
    print("ROI zips   :", os.path.normpath(os.path.join(folder, ROI_DIR_REL)))
    print("Hyperstacks:", os.path.normpath(os.path.join(folder, TIFF_DIR_REL)))

    per_cell, per_image, provenance = collect(folder)
    if not per_cell:
        print("No usable data found.")
        return
    print("\n%d images, %d cells" % (len(per_image), len(per_cell)))
    n_nuc = sum(1 for r in per_cell if r["nuc_area"] is not None)
    print("%d/%d cells have a nucleus ROI (%.1f%%)"
          % (n_nuc, len(per_cell), 100.0 * n_nuc / len(per_cell)))
    print("pixel sizes in use:", sorted({p["pixel_um"] for p in provenance}))

    tests = {m[0]: paired_tests(per_image, m[0]) for m in METRICS}
    desc = {m[0]: descriptive_tests(per_cell, m[0]) for m in METRICS}
    rtests = {m[0]: ratio_tests(per_image, m[0]) for m in METRICS}
    atests = {m[0]: allcell_tests(per_cell, per_image, m[0]) for m in METRICS}

    figures = []
    for key, ylabel, _short, ylim, stem in METRICS:
        paths = make_figure(per_image, tests[key], key, ylabel, ylim,
                            out_dir, stem, show)
        figures.append(paths[0])
        print("wrote", paths[0])
    paths = make_ratio_figure(per_image, rtests, out_dir, show)
    figures.append(paths[0])
    print("wrote", paths[0])
    paths = make_bars_figure(per_cell, per_image, tests, atests, out_dir, show)
    figures.append(paths[0])
    print("wrote", paths[0])

    caps = write_captions(out_dir, folder)
    print("wrote", caps)

    xlsx = os.path.join(out_dir, OUT_XLSX)
    write_workbook(xlsx, per_cell, per_image, provenance, tests, desc,
                   rtests, atests, figures)
    print("wrote", xlsx)

    print("\n=== paired t-test, SG+ vs SG- image means ===")
    for key, _ylabel, short, _ylim, _stem in METRICS:
        print("\n-- %s --" % short)
        for line, _, _ in CELL_LINES:
            for t in TIMEPOINTS:
                r = tests[key][(line, t)]
                if r["n_pairs"] < 2:
                    print("  %-14s %2d min  n=%d  (too few paired images)"
                          % (line, t, r["n_pairs"]))
                    continue
                print("  %-14s %2d min  n=%d  SG-=%8.3f  SG+=%8.3f  diff=%+8.3f  "
                      "p=%.4g %s"
                      % (line, t, r["n_pairs"], r["mean_neg"], r["mean_pos"],
                         r["mean_pos"] - r["mean_neg"], r["p"], sig_mark(r["p"])))

    print("\n=== SG+ / SG- ratio per image (z-independent) ===")
    for key, _ylabel, short, _ylim, _stem in METRICS:
        print("\n-- %s --" % short)
        for line, _, _ in CELL_LINES:
            for t in TIMEPOINTS:
                r = rtests[key][(line, t)]
                if r["n"] < 2:
                    print("  %-14s %2d min  n=%d  (too few paired images)"
                          % (line, t, r["n"]))
                    continue
                print("  %-14s %2d min  n=%d  ratio=%.3f (%+.1f%%)  p=%.4g %s"
                      % (line, t, r["n"], r["ratio"], (r["ratio"] - 1) * 100,
                         r["p"], sig_mark(r["p"])))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args[0] if args else CURRENT_FOLDER, show="--show" in sys.argv)
