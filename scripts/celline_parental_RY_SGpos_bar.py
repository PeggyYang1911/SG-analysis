"""
SG+ (%) BAR figure for PARENTAL RPMI 8226 + RPMI-RY only (DMSO vs BTZ).

A trimmed companion to celline_treatment_compare.py: it reuses that script's
data collection, pairing, paired t-test and bar drawing, but restricts the panel
to the two requested cell lines and forces the SG+ metric to render as bars
(the main script now draws SG+ as a line). Distinct output stem, so the 4-line
figures are left untouched.

The parental RPMI 8226 slide swap is preserved (see the memory note / the main
script): files named 'RPMI_DMSO_3h' are BTZ, 'RPMI_BTZ_1_uM_3hr' are DMSO.

Outputs into the same 5slices_subset folder:
  Fig_parental_RY_SGpositive_percent.png / .svg / .pdf  (+ _PPT.svg, Arial)
"""
import os
import matplotlib.pyplot as plt

import celline_treatment_compare as C

# two cell lines only; parental keeps its crossed (swapped) substrings
C.CELL_LINES = [
    ("RPMI 8226", "#0072B2", "RPMI_BTZ",     "RPMI_DMSO",   None, set(range(6, 11))),
    ("RPMI-RY",   "#E69F00", "RPMI-RY_DMSO", "RPMI_RY_BTZ", None, None),
]
C.DSER = {c[0]: c[4] for c in C.CELL_LINES}
C.BSER = {c[0]: c[5] for c in C.CELL_LINES}
C.FIGSIZE = (6.4, 6.8)          # narrower canvas for just two groups

STEM = "Fig_parental_RY_SGpositive_percent"
KEY, YLABEL, _S, YLIM, YTICKS = C.METRICS[0]   # sgpos


def render(font_family=None):
    # C.figure is the BAR builder (C.line_figure is the slope one)
    return C.figure(C.RECS, KEY, YLABEL, STEM, YLIM, YTICKS, font_family=font_family)


def main():
    C.RECS = C.collect()
    print(f"{'cell line':12} {'DMSO':>16} {'BTZ':>16}   paired p")
    for label, *_ in C.CELL_LINES:
        d, b = C.paired(C.RECS, label, KEY)
        _t, p, _n = C.ptest(C.RECS, label, KEY)
        print(f"{label:12} {f'{d.mean():5.1f}% (n={len(d)})':>16} "
              f"{f'{b.mean():5.1f}% (n={len(b)})':>16}   p={p:.5f} {C.sig(p)}")

    fig, _ = render()
    for ext in C.FORMATS:
        fig.savefig(os.path.join(C.OUT_DIR, f"{STEM}.{ext}"), dpi=C.DPI,
                    bbox_inches="tight")
    plt.close(fig)
    print("\nSaved:", os.path.join(C.OUT_DIR, f"{STEM}.png"), "(+svg/pdf)")

    if C.PPT_SVG:
        pf, _ = render(font_family=C.PPT_FONT)
        pf.savefig(os.path.join(C.OUT_DIR, f"{STEM}{C.PPT_SUFFIX}.svg"), dpi=C.DPI,
                   bbox_inches="tight")
        plt.close(pf)
        print("Saved:", os.path.join(C.OUT_DIR, f"{STEM}{C.PPT_SUFFIX}.svg"),
              f"({C.PPT_FONT})")


if __name__ == "__main__":
    main()
