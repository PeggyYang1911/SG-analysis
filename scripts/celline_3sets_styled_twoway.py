"""
Restyle the SG+ % line figure (selected "3 sets per group") to the Okabe-Ito
colour theme with RPMI8226 biological labels, and run a two-way ANOVA
(Cell line x NaAsO2 exposure, with interaction) for each metric.

Reads the "3 sets per group" sheet of ThreeCellLine_comparison.xlsx (metrics
already computed), so it matches the user's selected rows exactly.

Outputs into the comparison folder:
  Fig_line_3sets_SGpositive_percent.png   (RESTYLED, overwrites)
  TwoWay_ANOVA_3sets.xlsx                  (one sheet per metric)
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

FOLDER = r"J:\FATHO_LAB\YANG Pei-Tzu\Biological evaluation\Immunostaining\260601_PY_260522 Sterallis KMS WT and KO AS\v15 analysis\3celline_comparison_WT_KO_R"
WORKBOOK = os.path.join(FOLDER, "ThreeCellLine_comparison.xlsx")
SHEET = "3 sets per group"

# original label -> (display label, Okabe-Ito colour)
LINES = [
    ("KMS WT", "Parental",             "#0072B2"),  # blue
    ("KMS KO", "HDAC6 Knock-out",      "#E69F00"),  # orange
    ("KMS-R",  "bortezomib-resistant", "#009E73"),  # bluish green
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
    out = {}
    labels = {"C(CellLine)": "Cell line", "C(Time)": "NaAsO2",
              "C(CellLine):C(Time)": "Interaction"}
    for term, name in labels.items():
        out[name] = (float(tbl.loc[term, "F"]), float(tbl.loc[term, "PR(>F)"]))
    return out


def styled_figure(df, twoway_pct):
    plt.rcParams.update({"font.size": 18})
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
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
                           s=34, alpha=0.35, edgecolors="none", zorder=2)
        ax.errorbar(xs, means, yerr=sems, marker="o", ms=9, lw=2.6,
                    color=color, capsize=4, capthick=2, label=label, zorder=3)
    ax.set_xlabel("1 mM NaAsO$_2$ (minute)", fontsize=22)
    ax.set_ylabel("SG+ cells (%)", fontsize=24)
    ax.set_xticks(TIMEPOINTS)
    ax.set_ylim(-3, 108)
    ax.tick_params(labelsize=20, width=1.4, length=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.4)

    leg = ax.legend(title="KMS-12-BM", loc="upper left", fontsize=17,
                    title_fontsize=19, frameon=False,
                    handlelength=1.4, borderaxespad=0.3)
    leg._legend_box.align = "left"

    box = ("Two-way ANOVA\n"
           f"Cell line: {sig(twoway_pct['Cell line'][1])}\n"
           f"NaAsO$_2$: {sig(twoway_pct['NaAsO2'][1])}\n"
           f"Interaction: {sig(twoway_pct['Interaction'][1])}")
    # our data plateaus high, so the top-right is occupied -> put the box in the
    # empty lower-right region instead of the template's top-right.
    ax.text(0.97, 0.03, box, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=16, bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                                   edgecolor="0.5"))
    fig.tight_layout()
    out = os.path.join(FOLDER, "Fig_line_3sets_SGpositive_percent.png")
    fig.savefig(out, dpi=200)
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
        ws.cell(1, 1, f"Two-way ANOVA (balanced, n=3/cell)  -  {ylabel}").font = Font(bold=True)
        ws.append(["Factor", "F", "p_value", "Significance"])
        for c in range(1, 5):
            ws.cell(2, c).font = hdr_font
            ws.cell(2, c).fill = hdr_fill
            ws.cell(2, c).alignment = Alignment(horizontal="center")
        for name in ("Cell line", "NaAsO2", "Interaction"):
            F, p = results[key][name]
            ws.append([name, round(F, 4), round(p, 8), sig(p)])
        ws.append([])
        ws.append(["Model: y ~ C(CellLine) * C(Time);  type-II SS"])
        ws.append(["Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant"])
        for i, w in enumerate([16, 12, 14, 14], start=1):
            ws.column_dimensions[get_column_letter(i)].width = w
    out = os.path.join(FOLDER, "TwoWay_ANOVA_3sets.xlsx")
    wb.save(out)
    return out


def main():
    df = load_df()
    results = {key: two_way(df, key) for key, _ in METRICS}
    fig = styled_figure(df, results["SG_positive_(%)"])
    xlsx = write_stats(results)

    print("=== Two-way ANOVA (Cell line x NaAsO2, with interaction) ===")
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
