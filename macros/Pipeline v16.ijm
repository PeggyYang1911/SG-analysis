// ============================================================
// Pipeline v16.ijm  —  BATCH VERSION
//
// Processes all .tif files in a chosen input folder.
// For each image: segments nuclei from the DAPI Sum Z projection,
// measures G3BP1 per cell per Z-slice, counts punctae using
// full-image rolling-ball subtraction (no ring artifact).
//
// INCREMENTAL FILTERS (same as v16 single-image): candidates are ALL
// objects above the intensity threshold; each size/shape filter is
// applied only if its USE_* toggle is true. A punctum is counted only
// if it passes EVERY enabled filter. Set the parameters below to match
// whatever you tuned on the single-image macro.
//
// Pipeline per candidate:
//   1. intensity threshold  (PUNCTA_NOISE, always on)
//   2. size        if USE_SIZE      (PUNCTA_MIN_SIZE..PUNCTA_MAX_SIZE, um^2)
//   3. circularity if USE_CIRC      (>= PUNCTA_MIN_CIRC)
//   4. roundness   if USE_ROUND     (>= PUNCTA_MIN_ROUND)
//   5. solidity    if USE_SOLIDITY  (>= PUNCTA_MIN_SOLIDITY)
//
// CSV schema is unchanged (same columns) so the Python pipeline is
// unaffected; only the Punctae_Count value changes.
//
// Output per image (saved to chosen output folder):
//   <imagename>_v16.csv
//   <imagename>_v16_ROIs.zip
//   <imagename>_v16_punctae.csv   (only if DIAGNOSTIC = true)
// ============================================================

macro "Analyze G3BP1 RPMI Batch v16" {

// ── USER PARAMETERS ─────────────────────────────────────────
    DAPI_CH          = 1;
    G3BP1_CH         = 2;

    DAPI_BLUR        = 1;
    NUC_THRESHOLD    = "Triangle";
    NUC_ERODE_PX     = 4.5;
    NUC_MIN_SIZE     = 10;
    NUC_MAX_SIZE     = 350;
    NUC_MIN_CIRC     = 0.20;

    MAX_MATCH_DIST   = 20;
    EXCLUDE_EDGES    = true;

    CELL_EXPAND_PX   = 3;

    BG_RADIUS        = 10;
    PUNCTA_NOISE     = 500;   // intensity threshold above background (always on)

    // ---- puncta filters: each applied only if its USE_* is true ----
    USE_SIZE         = true;  // keep PUNCTA_MIN_SIZE <= area <= PUNCTA_MAX_SIZE
    PUNCTA_MIN_SIZE  = 0.5;    //   min area (um^2)
    PUNCTA_MAX_SIZE  = 1.5;    //   max area (um^2)

    USE_CIRC         = true;  // keep circularity >= PUNCTA_MIN_CIRC
    PUNCTA_MIN_CIRC  = 0.90;

    USE_ROUND        = true;  // keep roundness ("Round") >= PUNCTA_MIN_ROUND
    PUNCTA_MIN_ROUND = 0.50;  //   1=circle, ~0.5=2:1 rod

    USE_SOLIDITY     = false; // keep solidity >= PUNCTA_MIN_SOLIDITY
    PUNCTA_MIN_SOLIDITY = 0.90; // 1=convex, lower=bent/branched

    // ---- optional output ----
    DIAGNOSTIC       = true;  // write <name>_v16_punctae.csv (all candidates + Kept).
                              // Set false to skip it. Does NOT change Punctae_Count.
// ─────────────────────────────────────────────────────────────

    if (EXCLUDE_EDGES) { excludeStr = " exclude"; } else { excludeStr = ""; }

    run("Set Measurements...", "area mean standard min max shape redirect=None decimal=3");
    setOption("BlackBackground", true);
    run("Options...", "iterations=1 count=1 black");

    // ── Choose folders ───────────────────────────────────────
    inputDir  = getDirectory("Choose INPUT folder (images to analyse)");
    outputDir = getDirectory("Choose OUTPUT folder (CSV files)");

    fileList  = getFileList(inputDir);
    nFiles    = fileList.length;
    nDone     = 0;

    print("\\Clear");
    print("Input  : " + inputDir);
    print("Output : " + outputDir);
    print("Files found: " + nFiles);

    // ── File loop ────────────────────────────────────────────
    for (f = 0; f < nFiles; f++) {
        fileName = fileList[f];

        // Skip non-tif files and sub-folders
        lc = toLowerCase(fileName);
        if (!endsWith(lc, ".tif") && !endsWith(lc, ".tiff")) continue;

        print("\n── Processing: " + fileName + " ──");
        open(inputDir + fileName);
        origID = getImageID();
        title  = getTitle();
        getDimensions(imgW, imgH, nCh, nZ, nT);

        if (nCh < 2) {
            print("SKIP " + fileName + ": need 2 channels, found " + nCh);
            close();
            continue;
        }

        // Strip extension for output file names
        baseName = title;
        dotIdx   = lastIndexOf(baseName, ".");
        if (dotIdx > 0) baseName = substring(baseName, 0, dotIdx);

        // ── 1. Sum Z projection ──────────────────────────────
        selectImage(origID);
        run("Z Project...", "projection=[Sum Slices]");
        projID = getImageID();
        rename("SumZ_proj");

        selectImage(projID);
        run("Duplicate...", "title=DAPI_proj duplicate channels=" + DAPI_CH);
        dapiProjID = getImageID();

        // ── 2. Segment DAPI projection → reference centroids ─
        selectImage(dapiProjID);
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
        nRef = roiManager("count");
        print("Reference cells (projection): " + nRef);

        if (nRef == 0) {
            print("SKIP " + fileName + ": no nuclei found in projection.");
            selectImage(dapiProjID); close();
            selectImage(projID);     close();
            selectImage(origID);     close();
            continue;
        }

        refX = newArray(nRef);
        refY = newArray(nRef);
        for (i = 0; i < nRef; i++) {
            roiManager("select", i);
            roiManager("rename", "Cell_" + (i + 1));
            getSelectionBounds(bx, by, bw, bh);
            refX[i] = bx + bw / 2.0;
            refY[i] = by + bh / 2.0;
        }

        // Save reference ROI zip before Z-loop resets the manager
        roiPath = outputDir + baseName + "_v16_ROIs.zip";
        roiManager("deselect");
        roiManager("save", roiPath);

        // ── 3. Open CSV ──────────────────────────────────────
        outPath = outputDir + baseName + "_v16.csv";
        if (File.exists(outPath)) File.delete(outPath);
        File.append("File,Slice,Cell_ID,Cyto_Area_px2,Max,Mean,StdDev,CV_%,Punctae_Count", outPath);
        diagPath = outputDir + baseName + "_v16_punctae.csv";
        if (DIAGNOSTIC) {
            if (File.exists(diagPath)) File.delete(diagPath);
            File.append("File,Slice,Cell_ID,Punctum_Area,Circ,AR,Round,Solidity,Kept", diagPath);
        }
        print("Output: " + outPath);
        print("Reference cells: " + nRef + "  Slices: " + nZ);

        // ── 4. Per-slice loop ────────────────────────────────
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

            // Full-image background subtraction (no zeroing → no ring artifact)
            selectImage(g3SliceID);
            run("Duplicate...", "title=G3BP1_bgsub");
            bgsubID = getImageID();
            run("Subtract Background...", "rolling=" + BG_RADIUS);

            // Build cytoplasm ROIs
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

            // Measure and write (reverse order keeps indices stable when deleting)
            for (loc = nLocal - 1; loc >= 0; loc--) {
                cellID  = cellIDForLoc[loc];
                cytoIdx = cytoRoiIdx[loc];
                cIdx    = tempCellIdx[loc];

                if (cytoIdx >= 0) {
                    // Raw signal measurements
                    selectImage(g3SliceID);
                    roiManager("select", cytoIdx);
                    run("Clear Results");
                    run("Measure");
                    cytoArea = getResult("Area",   0);
                    cytoMean = getResult("Mean",   0);
                    cytoMax  = getResult("Max",    0);
                    cytoStd  = getResult("StdDev", 0);
                    if (cytoMean > 0) { cytoCV = (cytoStd / cytoMean) * 100; } else { cytoCV = 0; }

                    // Punctae count on full-image bgsub copy
                    selectImage(bgsubID);
                    run("Select None");
                    run("Duplicate...", "title=G3BP1_cyto");
                    g3cytoID = getImageID();
                    roiManager("select", cytoIdx);
                    getStatistics(area, mean, bgMin, bgMax, bgStd);

                    if (bgMax < PUNCTA_NOISE) {
                        punctaeCount = 0;
                    } else {
                        run("Make Inverse");
                        run("Set...", "value=0");
                        run("Select None");
                        setThreshold(PUNCTA_NOISE, 65535);
                        run("Convert to Mask");
                        run("Clear Results");
                        // candidates = ALL objects above threshold; filters applied below
                        run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 display");
                        // apply each enabled filter one by one
                        punctaeCount = 0;
                        for (rr = 0; rr < nResults; rr++) {
                            pArea = getResult("Area", rr);
                            pCirc = getResult("Circ.", rr);
                            pAR   = getResult("AR", rr);
                            pRnd  = getResult("Round", rr);
                            pSol  = getResult("Solidity", rr);
                            kept = true;
                            if (USE_SIZE && (pArea < PUNCTA_MIN_SIZE || pArea > PUNCTA_MAX_SIZE)) kept = false;
                            if (USE_CIRC && pCirc < PUNCTA_MIN_CIRC) kept = false;
                            if (USE_ROUND && pRnd < PUNCTA_MIN_ROUND) kept = false;
                            if (USE_SOLIDITY && pSol < PUNCTA_MIN_SOLIDITY) kept = false;
                            if (kept) punctaeCount++;
                            if (DIAGNOSTIC) {
                                keptFlag = 0; if (kept) keptFlag = 1;
                                File.append(title + "," + z + "," + cellID + "," +
                                    pArea + "," + pCirc + "," + pAR + "," + pRnd + "," +
                                    pSol + "," + keptFlag, diagPath);
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

            selectImage(bgsubID);     close();
            selectImage(dapiSliceID); close();
            selectImage(g3SliceID);   close();

        } // end Z-loop

        // ── 5. Clean up this image ───────────────────────────
        selectImage(dapiProjID); close();
        selectImage(projID);     close();
        selectImage(origID);     close();

        nDone++;
        print("Done: " + baseName + "  (" + nDone + " of " + nFiles + " tif files)");

    } // end file loop

    // Reload last ROI set so manager shows named reference cells
    if (nDone > 0) {
        roiManager("reset");
        roiManager("open", roiPath);
    }

    showMessage("Batch complete",
        "Processed " + nDone + " image(s).\nOutput folder: " + outputDir);
}
