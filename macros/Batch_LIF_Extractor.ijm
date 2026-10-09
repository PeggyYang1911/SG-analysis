// Batch LIF Extractor for Fiji/ImageJ
// Opens all .lif files in a directory, saves each series as an individual
// hyperstack TIFF in a new subfolder, then closes all images.
//
// Requirements: Bio-Formats plugin (bundled with Fiji)

macro "Batch LIF Extractor" {

    // --- Select input folder ---
    inputDir = getDirectory("Select folder containing .lif files");

    // --- Create output subfolder ---
    outputDir = inputDir + "Extracted_Hyperstacks" + File.separator;
    if (!File.exists(outputDir)) {
        File.makeDirectory(outputDir);
    }

    // --- Initialize Bio-Formats macro extensions ---
    run("Bio-Formats Macro Extensions");

    fileList = getFileList(inputDir);
    totalSaved = 0;

    print("\\Clear");
    print("=== Batch LIF Extractor ===");
    print("Input : " + inputDir);
    print("Output: " + outputDir);
    print("---------------------------");

    for (i = 0; i < fileList.length; i++) {

        // Process only .lif files (case-insensitive)
        if (!endsWith(toLowerCase(fileList[i]), ".lif")) continue;

        lifPath  = inputDir + fileList[i];
        lifBase  = substring(fileList[i], 0, lastIndexOf(fileList[i], "."));

        print("File " + (i+1) + ": " + fileList[i]);

        // --- Get number of series without opening an image window ---
        Ext.setId(lifPath);
        Ext.getSeriesCount(nSeries);
        Ext.close();

        print("  -> " + nSeries + " series found");

        // --- Open, save and close each series ---
        for (s = 1; s <= nSeries; s++) {

            run("Bio-Formats Importer",
                "open=[" + lifPath + "] " +
                "autoscale color_mode=Default " +
                "view=Hyperstack stack_order=XYCZT " +
                "series_" + s);

            // Build output filename from the window title
            imgTitle = getTitle();
            outName  = imgTitle;
            outName  = replace(outName, ".lif",  "");
            outName  = replace(outName, ".LIF",  "");
            outName  = replace(outName, " - ",   "_");
            outName  = replace(outName, " ",     "_");
            outName  = replace(outName, "/",     "_");
            outName  = replace(outName, ":",     "-");
            outName  = replace(outName, "|",     "_");

            // Safety fallback if title ends up empty
            if (lengthOf(outName) == 0) {
                outName = lifBase + "_series" + IJ.pad(s, 3);
            }

            outPath = outputDir + outName + ".tif";
            saveAs("Tiff", outPath);
            totalSaved++;

            print("  [" + s + "/" + nSeries + "] Saved: " + outName + ".tif");

            close();
        }
    }

    print("---------------------------");
    print("Done!  " + totalSaved + " hyperstacks saved.");

    showMessage("Batch Complete",
        "Saved " + totalSaved + " hyperstacks to:\n\n" + outputDir);
}
