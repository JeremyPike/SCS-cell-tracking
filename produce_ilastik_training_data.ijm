/* ImageJ macro script to generate training exampled for ilastik classifier for Lok at al. 
 *  
 *  Author: Jeremy Pike, Image analyst for COMPARE, j.a.pike@bham.ac.uk
*/


#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix


setBatchMode(true);
close("*");
count = 0;
processFolder(input);

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

	// point bio-formats at raw data file
	path = input + File.separator + file;
	run("Bio-Formats Macro Extensions");
	Ext.setId(path);
	Ext.getSeriesCount(seriesCount);
	Ext.getDimensionOrder(dimOrder);

	// loop thrugh series
	for (i = 1; i < =seriesCount; i++) { 
		
		Ext.setSeries(i - 1);
		Ext.getSizeC(sizeC);
		Ext.getSizeZ(sizeZ);
		Ext.getSizeT(sizeT);

		// only process movies if more than 10 slices and 10 time-points
		if (sizeZ >= 10 && sizeT >= 10) { 
			
			
			options = "open=[" + path + "] view=[Hyperstack] stack_order=" + dimOrder + " specify_range series_" + i ;
			// only load second channel
			cOpts = "c_begin_" + i + "=" + 2 + " c_end_" + i + "=" + sizeC + " c_step_" + i + "=1";
			// all slices
			zOpts = "z_begin_" + i + "=" + 1 + " z_end_" + i + "="+ sizeZ +" z_step_" + i + "=1";

			// alternate loading cropped movie at start, middle and end of full movie
			count ++;
			if (count%3 == 0) {
				tOpts = "t_begin_" + i + "=" + (sizeT - 1) + " t_end_" + i + "=" + sizeT + " t_step_" + i + "=1";
				options = options + " " + cOpts + " " + zOpts + " " + tOpts;
				print("test3");
			} else if (count%3 == 1) {
				
				tOpts = "t_begin_" + i + "=" + floor(sizeT/2) + " t_end_" + i + "=" + (floor(sizeT/2) + 1) + " t_step_" + i + "=1";
				options = options + " " + cOpts + " " + zOpts + " " + tOpts;
				print("test2");
			} else {
				tOpts = "t_begin_" + i + "=" + 1 + " t_end_" + i + "=" + 2 + " t_step_" + i + "=1";
				options = options + " " + cOpts + " " + zOpts + " " + tOpts;
				print("test1");
			}
			// load subset of data
			print(path + ":", i, sizeC, sizeZ, sizeT);
			print(options);
			print(tOpts);
		    run("Bio-Formats Importer", options);
		    rename("raw");
		    // downsample by factor of four
		    run("Scale...", "x=0.25 y=0.25 z=1.0 width=355 height=356 depth=17 interpolation=Bilinear average create");
		    rename(count);
		    close("raw");
		    // export as ilastik compatible hdf5 file
		    exportPath = output + File.separator + count + ".h5";
		    run("Export HDF5", "select=" + exportPath + " exportpath=" + exportPath + " datasetname=" + count + " compressionlevel=0");
			
			close(count);

		    
		

		}
	}
	Ext.close();
}
