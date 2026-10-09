"""
Styled SG+ % line figure (Okabe-Ito, RPMI8226 labels) + two-way ANOVA for the
5-SLICE-SUBSET three-cell-line comparison.

Reads the Summary_per_image sheet of ThreeCellLine_comparison_5slices.xlsx (the
full 5-slice image set), so it matches the 5-slice figure exactly. Design is
unbalanced (5-7 images/group) -> statsmodels type-II two-way ANOVA.

Outputs into the 5slices folder:
  Fig_line_5slices_SGpositive_percent_styled.png
  TwoWay_ANOVA_5slices.xlsx  (one sheet per metric)
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

FOLDER = r"J:\FATHO_LAB\YANG Pei-Tzu\Biological evaluation\Immunostaining\260601_PY_260522 Sterallis KMS WT and KO AS\v15 analysis\3celline_comparison_WT_KO_R_5slices"
WORKBOOK = os.path.join(FOLDER, "ThreeCellLine_comparison_5slices.xlsx")
SHEET = "Summary_per_image"

LINES = [
    ("KMS WT", "Parental",             "#0072B2"),
    ("KMS KO", "HDAC6 Knock-out",      "#E69F00"),
    ("KMS-R",  "bortezomib-resistant", "#009E73"),
]
TIMEPOINTS = [0, 15, 30, 45]
METRICS = [
    ("SG_positive_(%)", "SG+ cells (%)"),
    ("Mean_punctae_per_cell", "SG per cell"),
    ("Mean_punctae_per_SG+_cell", "SG per SG+ cell"),
]


def load_df():
    wb = openpyxl.load_workbook(WORKBOOK, data_only=True)
    ws = wb[SHEET]
    rows = list(ws.iter_rows(values_only=True))
    header = [str(h).strip() for h in rows[0]]
    idx = {h: i for i, h in enumerate(header)}
    recs = []
    for r in rows[1:]:
        if r is None or r[idx["Cell_line"]] is None:
            continue
        recs.append({
            "line": str(r[idx["Cell_line"]]).strip(),
            "tp": int(r[idx["Timepoint_min"]]),
            "SG_positive_(%)": float(r[idx["SG_positive_(%)"]]),
            "Mean_punctae_per_cell": float(r[idx["Mean_punctae_per_cell"]]),
            "Mean_punctae_per_SG+_cell": float(r[idx["Mean_punctae_per_SG+_cell"]]),
        })
    return pd.DataFrame(recs)


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


def two_way(df, key):
    d = df.rename(columns={key: "y"}).copy()
    d["CellLine"] = d["line"].astype("category")
    d["Time"] = d["tp"].astype("category")
    model = ols("y ~ C(CellLine) * C(Time)", data=d).fit()
    tbl = sm.stats.anova_lm(model, typ=2)
    labels = {"C(CellLine)": "Cell line", "C(Time)": "NaAsO2",
              "C(CellLine):C(Time)": "Interaction"}
    return {name: (float(tbl.loc[term, "F"]), float(tbl.loc[term, "PR(>F)"]))
            for term, name in labels.items()}


def styled_figure(df, tw):
    plt.rcParams.update({"font.size": 22, "font.family": "DejaVu Sans"})
    fig, ax = plt.subplots(figsize=(8.4, 6.6))
    rng = np.random.default_rng(0)
    for orig, label, color in LINES:
        xs, means, sems = [], [], []
        for tp in TIMEPOINTS:
            v = df[(df["line"] == orig) & (df["tp"] == tp)]["SG_positive_(%)"].values
            if len(v):
                xs.append(tp)
                means.append(v.mean())
                sems.append(v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)
                jit = rng.uniform(-1.2, 1.2, size=len(v))
                ax.scatter(np.array([tp] * len(v)) + jit, v, color=color,
                           s=40, alpha=0.35, edgecolors="none", zorder=2)
        ax.errorbar(xs, means, yerr=sems, marker="o", ms=11, lw=3.0,
                    color=color, capsize=4, capthick=2.2, label=label, zorder=3)
    ax.set_xlabel("1 mM NaAsO$_2$ (minute)", fontsize=26)
    ax.set_ylabel("SG+ cells (%)", fontsize=26)
    ax.set_xticks(TIMEPOINTS)
    # axis terminates at 100; legend + ANOVA box sit just ABOVE the plot area
    ax.set_ylim(-2, 101)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.tick_params(labelsize=24, width=1.6, length=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.6)

    leg = ax.legend(title="KMS-12-BM", loc="lower left", bbox_to_anchor=(0.0, 1.01),
                    fontsize=21, title_fontsize=22, frameon=False,
                    handlelength=1.4, borderaxespad=0.0)
    leg._legend_box.align = "left"

    box = ("Two-way ANOVA\n"
           f"Cell line: {sig(tw['Cell line'][1])}\n"
           f"NaAsO$_2$: {sig(tw['NaAsO2'][1])}\n"
           f"Interaction: {sig(tw['Interaction'][1])}")
    ax.text(1.0, 1.01, box, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=20, bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                                   edgecolor="0.5"))
    fig.tight_layout()
    out = os.path.join(FOLDER, "Fig_line_5slices_SGpositive_percent_styled.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def write_stats(results):
    wb = openpyxl.Workbook()
    hdr_font = Font(bold=True, color="FFFFFF")
    hdr_fill = PatternFill("solid", fgColor="404040")
    first = True
    for key, ylabel in METRICS:
        ws = wb.active if first else wb.create_sheet()
        ws.title = ylabel.replace("+", "pos").replace("%", "pct")[:31]
        first = False
        ws.cell(1, 1, f"Two-way ANOVA (type-II, unbalanced)  -  {ylabel}  [5-slice subset]").font = Font(bold=True)
        ws.append(["Factor", "F", "p_value", "Significance"])
        for c in range(1, 5):
            ws.cell(2, c).font = hdr_font
            ws.cell(2, c).fill = hdr_fill
            ws.cell(2, c).alignment = Alignment(horizontal="center")
        for name in ("Cell line", "NaAsO2", "Interaction"):
            F, p = results[key][name]
            ws.append([name, round(F, 4), round(p, 8), sig(p)])
        ws.append([])
        ws.append(["Model: y ~ C(CellLine) * C(Time);  type-II SS (unbalanced n=5-7/cell)"])
        ws.append(["Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant"])
        for i, w in enumerate([16, 12, 14, 14], start=1):
            ws.column_dimensions[get_column_letter(i)].width = w
    out = os.path.join(FOLDER, "TwoWay_ANOVA_5slices.xlsx")
    wb.save(out)
    return out


def main():
    df = load_df()
    for orig, label, _ in LINES:
        c = {tp: int(((df["line"] == orig) & (df["tp"] == tp)).sum()) for tp in TIMEPOINTS}
        print(f"  {label}: " + ", ".join(f"{tp}min n={c[tp]}" for tp in TIMEPOINTS))
    results = {key: two_way(df, key) for key, _ in METRICS}
    fig = styled_figure(df, results["SG_positive_(%)"])
    xlsx = write_stats(results)

    print("\n=== Two-way ANOVA (5-slice subset) ===")
    for key, ylabel in METRICS:
        print(f"\n{ylabel}:")
        for name in ("Cell line", "NaAsO2", "Interaction"):
            F, p = results[key][name]
            print(f"  {name:12s} F={F:8.3f}  p={p:.3e}  {sig(p)}")
    print("\nSaved:")
    print("  ", fig)
    print("  ", xlsx)


if __name__ == "__main__":
    main()
