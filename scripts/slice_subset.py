"""
Slice-subset helper for the SG pipeline.

Copies every *.csv in a folder into a `<N>slices_subset` subfolder, keeping only
the N slices with the HIGHEST cell coverage per file. Priority per slice = number
of distinct Cell_IDs it contains (the slice covering all cells first, then one
cell less, two cells less, ...); ties are broken by the lowest slice number.
Files with fewer than N slices keep all their slices.

Run this AFTER clean_csv.py. Then run sg_analysis.py on the created subfolder.

Usage:
    python slice_subset.py "PATH\\TO\\vNN analysis" [N]

N defaults to 4. Data rows are detected by an integer Cell_ID (works for .tif,
.czi, ... and skips the appended SUMMARY block); files with no data are skipped.
Originals are never modified -- output goes only to the new subfolder.
"""
import csv
import glob
import os
import sys

DEFAULT_N = 4


def is_data(row, header, cell_idx):
    if len(row) < len(header):
        return False
    try:
        int(float(row[cell_idx]))
        return True
    except (ValueError, IndexError):
        return False


def make_subset(folder, n):
    out = os.path.join(folder, f"{n}slices_subset")
    os.makedirs(out, exist_ok=True)
    files = sorted(glob.glob(os.path.join(folder, "*.csv")))
    print(f"{'File':30} {'#cellIDs':>8}  selected_slices (cells)")
    for path in files:
        with open(path, newline="") as f:
            rows = list(csv.reader(f))
        if not rows:
            continue
        h = rows[0]
        si, ci = h.index("Slice"), h.index("Cell_ID")
        data = [r for r in rows[1:] if is_data(r, h, ci)]
        name = os.path.basename(path)
        if not data:
            print(f"{name:30} SKIP (no data rows)")
            continue

        all_ids = {r[ci] for r in data}
        perslice = {}
        for r in data:
            perslice.setdefault(r[si], set()).add(r[ci])
        # highest coverage first, ties by lowest slice number
        ranked = sorted(perslice.items(), key=lambda kv: (-len(kv[1]), int(kv[0])))
        chosen = {s for s, _ in ranked[:n]}

        out_rows = [h] + [r for r in data if r[si] in chosen]
        with open(os.path.join(out, name), "w", newline="") as f:
            csv.writer(f).writerows(out_rows)

        desc = ", ".join(f"{s}({len(perslice[s])})" for s, _ in ranked[:n])
        print(f"{name:30} {len(all_ids):>8}  {desc}")

    print(f"\nNew subfolder: {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit('Usage: python slice_subset.py "PATH\\TO\\vNN analysis" [N]')
    folder = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_N
    make_subset(folder, n)
