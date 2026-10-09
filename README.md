# Stress-granule (SG) quantification — analysis code

Code for quantifying stress granules (G3BP1 puncta) in confocal immunostaining
of multiple myeloma cell lines (KMS-12-BM and RPMI 8226, including HDAC6-KO and
bortezomib-resistant variants), treated with sodium arsenite (NaAsO₂ time course)
or bortezomib / other drugs.

This repository holds **only the code** — the Fiji/ImageJ macros that detect and
count puncta, and the Python scripts that clean, subset, score and plot the
per-image CSVs. Raw images and analysis outputs are **not** tracked (see
`.gitignore`).

## Layout

```
macros/     Fiji/ImageJ macros (.ijm) — run inside Fiji
scripts/    Python analysis + figure scripts
```

## Pipeline

1. **Extract** per-Series hyperstacks from the raw `.lif` files
   (`macros/Batch_LIF_Extractor.ijm`).
2. **Detect & count** puncta per cell per Z-slice with the Fiji macro
   → one CSV per image (`macros/Pipeline v16.ijm`, batch).
   Use `Pipeline v16-single imag for tuning.ijm` to tune parameters on one open
   image (it can also export per-section QC overlays and a per-punctum
   diagnostic CSV).
3. **Clean** → `scripts/clean_csv.py <folder>` (drops `Cell_ID==0`, drops slices
   with < 5 cells; edits CSVs in place).
4. **Subset** → `scripts/slice_subset.py <folder> 5` (keeps the 5 z-slices with
   the most distinct cells per image; copies to `5slices_subset/`).
5. **Score** → `scripts/sg_analysis.py <folder>/5slices_subset` (writes
   `SG_analysis_summary.xlsx` + SG+% / SG-per-cell figures; paired t-test vs
   control).
6. **Compare** → the `celline_*` / `rpmi_*` figure scripts overlay several
   conditions on shared axes.

### Counting rules (invariants)

- Total cells = unique `Cell_ID`; `Cell_ID == 0` is background, always excluded.
- Per-cell puncta = **sum** of `Punctae_Count` across that cell's kept slices.
- A cell is **SG+** if that sum > 0.
- The statistical unit is the **image (Series)**, not the cell.
- Detect data rows by an integer `Cell_ID`, never by a `.tif` filename filter.

### Puncta detection (Fiji macro v16)

Candidates are all objects above an intensity threshold (`PUNCTA_NOISE`), then
each filter is applied only if its `USE_*` toggle is on:

| filter | keeps |
|---|---|
| size      | `PUNCTA_MIN_SIZE` ≤ area ≤ `PUNCTA_MAX_SIZE` (µm²) |
| circularity | Circ ≥ `PUNCTA_MIN_CIRC` |
| roundness | Round ≥ `PUNCTA_MIN_ROUND` |
| solidity  | Solidity ≥ `PUNCTA_MIN_SOLIDITY` |

## Requirements

- **Fiji / ImageJ** with Bio-Formats (for the macros).
- **Python 3** with `numpy`, `scipy`, `matplotlib`, `openpyxl`, `statsmodels`.
  ```
  pip install numpy scipy matplotlib openpyxl statsmodels
  ```

## Notes when moving between computers (Windows ↔ macOS)

- The figure scripts and the macros contain **hard-coded absolute paths** (e.g.
  `J:\...` on Windows) in their CONFIG / header blocks. Update these to the local
  paths (e.g. `/Volumes/...` or `~/...` on macOS) before running.
- `scripts/sg_provenance.py` is a shared helper imported by the counting scripts;
  keep it in the same folder as the scripts that import it.
- Paths with spaces must be quoted.
