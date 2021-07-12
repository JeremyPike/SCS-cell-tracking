/* ImageJ macro script to run ilastik pixel classifier on movies for Lok at al. 
 *  also saves some useful info on the files as a csv file
 *  
 *  Author: Jeremy Pike, Image analyst for COMPARE, j.a.pike@bham.ac.uk
*/


#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ File (label = "Ilastik project file", style = "file") ilp
#@ String (label = "File suffix", value = ".tif") suffix


setBatchMode(false);
close("*");

// clear results table
count = 0;
run("Clear Results");

// process all data
processFolder(input);

// save summary table of data as csv file
saveAs("Results", input + File.separator + "meta.csv");

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {

	
	path = input + File.separator + file;
	
	// point bio-formats at raw data file
	run("Bio-Formats Macro Extensions");
	Ext.setId(path);
	Ext.getSeriesCount(seriesCount);
	Ext.getDimensionOrder(dimOrder);

	// loop thrugh series
	for (i = 1; i < = seriesCount; i++) { 
		
		Ext.setSeries(i - 1);
		Ext.getSizeC(sizeC);
		Ext.getSizeZ(sizeZ);
		Ext.getSizeT(sizeT);
		Ext.getSeriesName(seriesName);
		
		
		// only process movies if more than 10 slices and 10 time-points
		if (sizeZ >= 10 && sizeT >= 10) { 
	
			// dont open first channel load everything else
			options = "open=[" + path + "] view=[Hyperstack] stack_order=" + dimOrder + " specify_range series_" + i ;
			cOpts = "c_begin_" + i + "=" + 2 + " c_end_" + i + "=" + sizeC + " c_step_" + i + "=1";
			zOpts = "z_begin_" + i + "=" + 1 + " z_end_" + i + "=" + sizeZ + " z_step_" + i + "=1";
			tOpts = "t_begin_" + i + "=" + 1 + " t_end_" + i + "=" + sizeT + " t_step_" + i + "=1";
			options = options + " " + cOpts + " " + zOpts + " " + tOpts;
			print(path + ":", i, sizeC, sizeZ, sizeT);
		    run("Bio-Formats Importer", options);
 			rename("raw");
		   
		    // output path for ilastik probability file
		   	exportPath = output + File.separator + replace(file, suffix, "") + "_" + seriesName + "_S" + i + "_C" + sizeC + "_Z" + sizeZ + "_T" + sizeT + "_CSIlastikProb.tif";
			print(exportPath);
			
		   	// put some useful info in the file in the results table
			getVoxelSize(width, height, depth, unit);
			setResult("filePath", count, input + File.separator + file);
			setResult("ilastikPath", count, exportPath);
			setResult("fileName", count, replace(file, suffix, ""));
		    setResult("seriesName", count, seriesName);
		    setResult("seriesIndex", count, i);
		   	setResult("sizeC", count, sizeC);
		    setResult("sizeZ", count, sizeZ);
		    setResult("sizeT", count, sizeT);
		    setResult("sizeC", count, sizeC);
		    setResult("voxelX", count, width);
		    setResult("voxelY", count, height);
		    setResult("voxelZ", count, depth);

		    // downsample data by factor of 4 and close raw data
		  	run("Scale...", "x=0.25 y=0.25 z=1.0 interpolation=Bilinear average create");
		    rename("scaled");
		    close("raw");

		    // put scaled voxel size in results table
		    getVoxelSize(width, height, depth, unit);
		    setResult("voxelX_scaled", count, width);
		    setResult("voxelY_scaled", count, height);
		    setResult("voxelZ_scaled", count, depth);
		    updateResults();
		    
			// run pixel classifer (need ilastik plugin installed and configured)
		    pixelClassificationArgs = "projectfilename=" + ilp + " saveonly=false inputimage=scaled pixelclassificationtype=Probabilities";
		    run("Run Pixel Classification Prediction", pixelClassificationArgs);
			
			close("scaled");
		  	// save probability file
		  	saveAs("Tiff", exportPath);
		
			close("*");
		    // add one to count
			count ++;

		}
	}
	Ext.close();
}
