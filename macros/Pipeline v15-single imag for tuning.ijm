// ============================================================
// Pipeline v15 — SINGLE IMAGE (for tuning)
//
// Runs on the image already OPEN and active in Fiji (open the .tif
// first, then run this macro). Same detection engine as the batch
// v15, for iterating parameters on one field.
//
// NEW in v15: a SHAPE FILTER on punctae. After thresholding and the
// circularity test, each candidate punctum must also have roundness
// >= PUNCTA_MIN_ROUND (axis-ratio based: 1 = circle, ~0.5 = 2:1 rod),
// which rejects elongated / rod-shaped bright patches that slip past
// circularity. Round is robust to the pixelation that makes
// circularity unreliable on small objects.
//
// Ring-artifact fix (from v14): rolling-ball runs on the FULL G3BP1
// slice (no zeroing) so no bright ridge forms at the nucleus boundary.
//
// CSV schema is unchanged (same columns); only Punctae_Count changes.
//
// Tuning aids (toggle in USER PARAMETERS):
//   DIAGNOSTIC     -> <name>_v15_punctae.csv: one row per candidate punctum
//                     (Area, Circ, AR, Round, Solidity) for shape analysis.
//   SAVE_SG_IMAGES -> <name>_SG_sections/zN_SG.png: per-section raw G3BP1 with
//                     the puncta COUNTED as SG filled in yellow (visual check).
// ============================================================

macro "Analyze G3BP1 RPMI Single v15" {

// ── USER PARAMETERS ─────────────────────────────────────────
    DAPI_CH          = 1;
    G3BP1_CH         = 2;

    DAPI_BLUR        = 1;
    NUC_THRESHOLD    = "Triangle";
    NUC_ERODE_PX     = 4.5;
    NUC_MIN_SIZE     = 20;
    NUC_MAX_SIZE     = 350;
    NUC_MIN_CIRC     = 0.20;

    MAX_MATCH_DIST   = 20;
    EXCLUDE_EDGES    = true;

    CELL_EXPAND_PX   = 3;

    BG_RADIUS        = 10;
    PUNCTA_NOISE     = 550;  // minimum intensity above background to count as a punctum
    PUNCTA_MIN_SIZE  = 0.5;   // min area (px²)
    PUNCTA_MAX_SIZE  = 1.5;     // max area (px²)
    PUNCTA_MIN_CIRC  = 0.90;  // min circularity (0=any shape, 1=perfect circle)
    PUNCTA_MIN_ROUND = 0.00;  // min roundness ("Round": 1=circle, ~0.5=2:1 rod).
                              // Rejects elongated/rod-shaped patches that pass
                              // circularity. Lower = more permissive; 0 = off.

    DIAGNOSTIC       = true;  // also write <name>_v15_punctae.csv: one row per
                              // detected punctum (Area, Circ, AR, Round, Solidity),
                              // logging ALL candidates regardless of the Round gate,
                              // for tuning the shape filter. Does NOT change the
                              // Punctae_Count in the main CSV.

    SAVE_SG_IMAGES   = true;  // save one PNG per Z-section into <name>_SG_sections/:
                              // the raw G3BP1 slice with the puncta COUNTED as SG
                              // (those passing size, circularity AND roundness)
                              // filled in yellow, for visual verification.
// ─────────────────────────────────────────────────────────────

    if (EXCLUDE_EDGES) { excludeStr = " exclude"; } else { excludeStr = ""; }

    run("Set Measurements...", "area mean standard min max shape redirect=None decimal=3");
    setOption("BlackBackground", true);
    run("Options...", "iterations=1 count=1 black");

    origID = getImageID();
    title  = getTitle();
    getDimensions(imgW, imgH, nCh, nZ, nT);

    if (nCh < 2)
        exit("ERROR: Need at least 2 channels. Found: " + nCh);

    // ── 1. Sum Z projection ──────────────────────────────────
    selectImage(origID);
    run("Z Project...", "projection=[Sum Slices]");
    projID = getImageID();
    rename("SumZ_proj");

    selectImage(projID);
    run("Duplicate...", "title=DAPI_proj duplicate channels=" + DAPI_CH);
    dapiProjID = getImageID();

    // ── 2. Segment DAPI projection → reference centroids ────
    selectImage(dapiProjID);
    getStatistics(area, mean, rawMin, rawMax, std);
    print("DAPI proj: min=" + rawMin + " max=" + rawMax + " mean=" + mean);
    run("Gaussian Blur...", "sigma=" + DAPI_BLUR);
    getStatistics(area, mean, blurMin, blurMax, std);
    run("32-bit");
    run("Subtract...", "value=" + blurMin);
    if (blurMax > blurMin)
        run("Multiply...", "value=" + d2s(65535.0 / (blurMax - blurMin), 6));
    run("16-bit");
    setAutoThreshold(NUC_THRESHOLD + " dark");
    run("Convert to Mask");
    run("Fill Holes");
    if (NUC_ERODE_PX > 0) {
        run("Options...", "iterations=" + NUC_ERODE_PX + " count=1 black");
        run("Erode");
        run("Options...", "iterations=1 count=1 black");
    }
    run("Watershed");

    roiManager("reset");
    run("Analyze Particles...",
        "size=" + NUC_MIN_SIZE + "-" + NUC_MAX_SIZE +
        " circularity=" + NUC_MIN_CIRC + "-1.00" + excludeStr + " add");
    nRef = roiManager("count");
    print("Reference cells (projection): " + nRef);

    if (nRef == 0)
        exit("No nuclei found in projection. Adjust segmentation parameters.");

    refX = newArray(nRef);
    refY = newArray(nRef);
    for (i = 0; i < nRef; i++) {
        roiManager("select", i);
        roiManager("rename", "Cell_" + (i + 1));
        getSelectionBounds(bx, by, bw, bh);
        refX[i] = bx + bw / 2.0;
        refY[i] = by + bh / 2.0;
    }

    // Ask user to choose output folder
    dir      = getDirectory("Choose folder to save results");
    // Derive base name from image title (reliable even without a saved file path)
    baseName = title;
    dotIdx   = lastIndexOf(baseName, ".");
    if (dotIdx > 0) baseName = substring(baseName, 0, dotIdx);

    // Save reference ROI zip now — before the Z-loop resets the manager
    roiPath  = dir + baseName + "_v15_ROIs.zip";
    roiManager("deselect");
    roiManager("save", roiPath);
    print("ROI set saved: " + roiPath);

    // ── 3. Open CSV ──────────────────────────────────────────
    outPath  = dir + baseName + "_v15.csv";
    if (File.exists(outPath)) File.delete(outPath);
    File.append("File,Slice,Cell_ID,Cyto_Area_px2,Max,Mean,StdDev,CV_%,Punctae_Count", outPath);
    diagPath = dir + baseName + "_v15_punctae.csv";
    if (DIAGNOSTIC) {
        if (File.exists(diagPath)) File.delete(diagPath);
        File.append("File,Slice,Cell_ID,Punctum_Area,Circ,AR,Round,Solidity", diagPath);
    }
    if (SAVE_SG_IMAGES) {
        sgDir = dir + baseName + "_SG_sections" + File.separator;
        File.makeDirectory(sgDir);
    }
    print("\\Clear");
    print("Output: " + outPath);
    print("Reference cells: " + nRef + "  Slices: " + nZ);

    // ── 4. Per-slice loop ────────────────────────────────────
    for (z = 1; z <= nZ; z++) {

        // Extract and segment DAPI slice
        selectImage(origID);
        run("Duplicate...", "title=DAPI_slice duplicate channels=" + DAPI_CH +
            " slices=" + z + " frames=1");
        dapiSliceID = getImageID();

        selectImage(dapiSliceID);
        getStatistics(area, mean, rawMin, rawMax, std);
        run("Gaussian Blur...", "sigma=" + DAPI_BLUR);
        getStatistics(area, mean, blurMin, blurMax, std);
        run("32-bit");
        run("Subtract...", "value=" + blurMin);
        if (blurMax > blurMin)
            run("Multiply...", "value=" + d2s(65535.0 / (blurMax - blurMin), 6));
        run("16-bit");
        setAutoThreshold(NUC_THRESHOLD + " dark");
        run("Convert to Mask");
        run("Fill Holes");
        if (NUC_ERODE_PX > 0) {
            run("Options...", "iterations=" + NUC_ERODE_PX + " count=1 black");
            run("Erode");
            run("Options...", "iterations=1 count=1 black");
        }
        run("Watershed");

        roiManager("reset");
        run("Analyze Particles...",
            "size=" + NUC_MIN_SIZE + "-" + NUC_MAX_SIZE +
            " circularity=" + NUC_MIN_CIRC + "-1.00" + excludeStr + " add");
        nLocal = roiManager("count");
        print("Z=" + z + " local nuclei=" + nLocal);

        // Extract G3BP1 slice
        selectImage(origID);
        run("Duplicate...", "title=G3BP1_slice duplicate channels=" + G3BP1_CH +
            " slices=" + z + " frames=1");
        g3SliceID = getImageID();

        // Pre-compute background-subtracted image once per slice (full image,
        // no zeroing) so the rolling ball has no hard boundary artifacts.
        selectImage(g3SliceID);
        run("Duplicate...", "title=G3BP1_bgsub");
        bgsubID = getImageID();
        run("Subtract Background...", "rolling=" + BG_RADIUS);

        // Per-slice canvas that accumulates the counted-SG puncta (filled white).
        if (SAVE_SG_IMAGES) {
            newImage("SG_mask", "8-bit black", imgW, imgH, 1);
            sgMaskID = getImageID();
        }

        // Build cytoplasm ROIs and track indices
        cytoRoiIdx   = newArray(nLocal);
        cellIDForLoc = newArray(nLocal);
        tempCellIdx  = newArray(nLocal);
        for (loc = 0; loc < nLocal; loc++) {
            cytoRoiIdx[loc]   = -1;
            cellIDForLoc[loc] = 0;
            tempCellIdx[loc]  = -1;
        }

        for (loc = 0; loc < nLocal; loc++) {
            roiManager("select", loc);
            getSelectionBounds(bx, by, bw, bh);
            cx = bx + bw / 2.0;
            cy = by + bh / 2.0;

            minDist = 999999;
            cellID  = 0;
            for (r = 0; r < nRef; r++) {
                dx = cx - refX[r];
                dy = cy - refY[r];
                dist = sqrt(dx*dx + dy*dy);
                if (dist < minDist) {
                    minDist = dist;
                    cellID  = r + 1;
                }
            }
            if (minDist > MAX_MATCH_DIST) {
                print("Z=" + z + " nucleus " + (loc+1) +
                      " unmatched (dist=" + d2s(minDist,1) + ")");
                cellID = 0;
            }
            cellIDForLoc[loc] = cellID;

            selectImage(g3SliceID);
            roiManager("select", loc);
            run("Enlarge...", "enlarge=" + CELL_EXPAND_PX);
            roiManager("add");
            cIdx = roiManager("count") - 1;
            tempCellIdx[loc] = cIdx;

            roiManager("select", newArray(loc, cIdx));
            roiManager("XOR");

            if (selectionType() == -1) {
                roiManager("select", cIdx); roiManager("delete");
                tempCellIdx[loc] = -1;
                print("Z=" + z + " Cell=" + cellID + " SKIPPED: empty cytoplasm");
            } else {
                roiManager("add");
                cytoRoiIdx[loc] = roiManager("count") - 1;
            }
        }

        // Measure and write — iterate in reverse to keep indices stable when deleting
        for (loc = nLocal - 1; loc >= 0; loc--) {
            cellID  = cellIDForLoc[loc];
            cytoIdx = cytoRoiIdx[loc];
            cIdx    = tempCellIdx[loc];

            if (cytoIdx >= 0) {
                // ── Measure raw signal ───────────────────────────
                selectImage(g3SliceID);
                roiManager("select", cytoIdx);
                run("Clear Results");
                run("Measure");
                cytoArea = getResult("Area",   0);
                cytoMean = getResult("Mean",   0);
                cytoMax  = getResult("Max",    0);
                cytoStd  = getResult("StdDev", 0);
                if (cytoMean > 0) { cytoCV = (cytoStd / cytoMean) * 100; } else { cytoCV = 0; }

                // ── Punctae count — no zeroing before rolling ball ─
                // Duplicate the FULL bgsubID with no selection active so
                // Fiji does not crop to a bounding box, which would shift
                // coordinates when we later re-select the cytoplasm ROI.
                selectImage(bgsubID);
                run("Select None");
                run("Duplicate...", "title=G3BP1_cyto");
                g3cytoID = getImageID();

                // Measure bgMax within cytoplasm on the full-image copy.
                roiManager("select", cytoIdx);
                getStatistics(area, mean, bgMin, bgMax, bgStd);

                if (bgMax < PUNCTA_NOISE) {
                    punctaeCount = 0;
                } else {
                    // Zero outside cytoplasm (rolling ball already ran on the
                    // full image, so this zero does not create a ring artifact).
                    run("Make Inverse");
                    run("Set...", "value=0");
                    run("Select None");
                    setThreshold(PUNCTA_NOISE, 65535);
                    run("Convert to Mask");
                    run("Clear Results");
                    apOpt = "size=" + PUNCTA_MIN_SIZE + "-" + PUNCTA_MAX_SIZE +
                            " circularity=" + PUNCTA_MIN_CIRC + "-1.00 display";
                    if (SAVE_SG_IMAGES) {
                        startIdx = roiManager("count");   // puncta ROIs appended after this
                        apOpt = apOpt + " add";
                    }
                    run("Analyze Particles...", apOpt);
                    // shape filter: keep only round puncta, reject rods
                    punctaeCount = 0;
                    for (rr = 0; rr < nResults; rr++) {
                        rnd = getResult("Round", rr);
                        kept = (rnd >= PUNCTA_MIN_ROUND);
                        if (kept) punctaeCount++;
                        if (DIAGNOSTIC) {
                            File.append(title + "," + z + "," + cellID + "," +
                                getResult("Area", rr) + "," +
                                getResult("Circ.", rr) + "," +
                                getResult("AR", rr) + "," + rnd + "," +
                                getResult("Solidity", rr), diagPath);
                        }
                        // paint kept SG puncta onto the per-slice mask
                        if (SAVE_SG_IMAGES && kept) {
                            selectImage(sgMaskID);
                            roiManager("select", startIdx + rr);
                            setForegroundColor(255, 255, 255);
                            run("Fill", "slice");
                            run("Select None");
                        }
                    }
                    if (SAVE_SG_IMAGES) {
                        // remove the puncta ROIs we appended, restoring the manager
                        while (roiManager("count") > startIdx) {
                            roiManager("select", roiManager("count") - 1);
                            roiManager("delete");
                        }
                    }
                }
                selectImage(g3cytoID); close();

                csvLine = title + "," + z + "," + cellID + "," +
                          cytoArea + "," + cytoMax + "," + cytoMean + "," +
                          cytoStd + "," + cytoCV + "," + punctaeCount;
                File.append(csvLine, outPath);
                print("Z=" + z + " Cell=" + cellID +
                      " Mean=" + cytoMean + " Punctae=" + punctaeCount);

                roiManager("select", cytoIdx); roiManager("delete");
            }

            if (cIdx >= 0) {
                roiManager("select", cIdx); roiManager("delete");
            }
        }

        // Save this section's counted-SG puncta overlaid (yellow) on raw G3BP1.
        if (SAVE_SG_IMAGES) {
            selectImage(g3SliceID);
            run("Select None");
            run("Duplicate...", "title=SG_overlay");
            sgOverlayID = getImageID();
            run("Enhance Contrast", "saturated=0.35");
            run("8-bit");
            run("RGB Color");
            selectImage(sgMaskID);
            getStatistics(mArea, mMean, mMin, mMax);
            if (mMax > 0) {
                setThreshold(128, 255);
                run("Create Selection");
                resetThreshold();
                selectImage(sgOverlayID);
                run("Restore Selection");
                setForegroundColor(255, 255, 0);
                run("Fill", "slice");
                run("Select None");
            }
            selectImage(sgOverlayID);
            saveAs("PNG", sgDir + "z" + z + "_SG.png");
            close();
            selectImage(sgMaskID); close();
        }

        selectImage(bgsubID);     close();
        selectImage(dapiSliceID); close();
        selectImage(g3SliceID);   close();

    } // end Z-loop

    // ── 5. Clean up projection images ───────────────────────
    selectImage(dapiProjID); close();
    selectImage(projID);     close();

    // Reload reference ROIs so Cell_1…Cell_N are visible in ROI Manager
    roiManager("reset");
    roiManager("open", roiPath);

    print("Done.");
    print("Results : " + outPath);
    print("ROI set : " + roiPath);
    showMessage("Analysis complete",
        "Results : " + outPath + "\nROI set : " + roiPath);
}
