"""
Step 1 - Clean v14 cytoplasm CSV files in a folder.

For every *.csv in the target folder:
  1. Remove all rows with Cell_ID == 0.
  2. Remove an entire Slice when it has fewer than 5 cells remaining
     (counted AFTER the Cell_ID==0 rows are removed).

Files are edited in place. Any trailing "SUMMARY" block written by
sg_analysis.py is ignored (only genuine data rows -- File column ending in
.tif -- are processed), so this can be run before or after the analysis step.

Usage:
    python clean_csv.py "PATH\\TO\\v14 analysis"

If no path is given it uses the CURRENT_FOLDER constant below.
"""
import csv
import glob
import os
import sys

# Edit this if you prefer to hard-code a default folder.
CURRENT_FOLDER = r"J:\FATHO_LAB\YANG Pei-Tzu\Biological evaluation\Immunostaining\260518_PY_ Sterallis autofluorescent 260508 RPMI 1-45 min As\Extracted_Hyperstacks\RPMI\v14 analysis"

MIN_CELLS_PER_SLICE = 5


def clean_folder(folder):
    files = sorted(glob.glob(os.path.join(folder, "*.csv")))
    print(f"Found {len(files)} CSV files in:\n  {folder}\n")

    for path in files:
        with open(path, "r", newline="") as f:
            rows = list(csv.reader(f))
        if not rows:
            continue

        header = rows[0]
        slice_idx = header.index("Slice")
        cell_idx = header.index("Cell_ID")

        # A genuine data row has the full column count and an integer Cell_ID.
        # (Do NOT key off the File column ending in .tif -- some datasets use a
        # different File value, and that would wrongly drop every row.)
        def is_data(r):
            if len(r) < len(header):
                return False
            try:
                int(float(r[cell_idx]))
                return True
            except (ValueError, IndexError):
                return False

        data = [r for r in rows[1:] if is_data(r)]
        has_content = any(any(c.strip() for c in r) for r in rows[1:])

        # SAFETY GUARD: never overwrite a file we couldn't parse. If the file has
        # non-empty rows but none look like data, leave it completely untouched.
        if not data:
            name = os.path.basename(path)
            if has_content:
                print(f"!! ABORT (unparseable rows; file left UNTOUCHED): {name}")
            else:
                print(f"(already header-only, skipped): {name}")
            continue

        # 1. remove rows with Cell_ID == 0
        kept = [r for r in data if r[cell_idx].strip() != "0"]
        removed_cell0 = len(data) - len(kept)

        # 2. count cells per slice, drop slices with < MIN_CELLS_PER_SLICE
        counts = {}
        for r in kept:
            counts[r[slice_idx]] = counts.get(r[slice_idx], 0) + 1
        valid = {s for s, c in counts.items() if c >= MIN_CELLS_PER_SLICE}
        dropped = sorted(set(counts) - valid, key=lambda x: int(x))

        final = [r for r in kept if r[slice_idx] in valid]

        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(final)

        print(f"{os.path.basename(path)}")
        print(f"  rows: {len(data)} -> {len(final)} | "
              f"Cell_ID=0 removed: {removed_cell0} | "
              f"slices dropped (<{MIN_CELLS_PER_SLICE} cells): {dropped}")


if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else CURRENT_FOLDER
    clean_folder(folder)
