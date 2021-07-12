

/* ImageJ groovy script to perform automated cell tracking for Lok at al. 
 *  
 *  Makes extenstive use of the TrackMate plugin: 
 *  Tinevez, J.-Y., Perry, N., Schindelin, J., Hoopes, G. M., Reynolds, G. D., Laplantine, E., … Eliceiri, K. W. (2017). 
 *  TrackMate: An open and extensible platform for single-particle tracking. 
 *  Methods, 115, 80–90. doi:10.1016/j.ymeth.2016.09.016
 *  
 *  Author: Jeremy Pike, Image analyst for COMPARE, j.a.pike@bham.ac.uk
*/



#@File (label="Select a directory containing the files to process", style="directory") directory
#@File (label="Select an output directory for the excel files", style="directory") outDirectory
#@String(label="Specify a file extension", value="tif") ext
#@Double(label="Specify spot radius (microns)", value=10) spotRadius
#@Double(label="Specify quality threshold for spots", value=15) spotThreshQual
#@Double(label="Specify mean intensity threshold for spots", value=0) spotThreshMeanInt
#@Double(label="Specify maximum track linking distance (microns)", value=20) maxLinkDist
#@Double(label="Specify the max linking distance (gap closing)", value= 20) maxLinkingDistGC
#@Integer(label="Specify the max gap in frames", value=2) maxFrameGap

#@Integer(label="Specify the number of time-points", value=36) t_crop
#@Integer(label="Specify the  number of frames", value=13) z_crop
#@Boolean(label="Display tracks", value=false) display
#@Boolean(label="Save spot statistics", value=true) save


import java.io.File
import java.util.Collection

import org.apache.commons.io.FileUtils

import ij.IJ
import ij.ImagePlus
import ij.plugin.HyperStackConverter
import ij.measure.Calibration
import ij.ImageStack
import ij.process.ImageProcessor
import ij.gui.WaitForUserDialog
import ij.measure.ResultsTable
import ij.plugin.Duplicator

import fiji.plugin.trackmate.Model
import fiji.plugin.trackmate.Settings
import fiji.plugin.trackmate.detection.LogDetectorFactory
import fiji.plugin.trackmate.features.FeatureFilter
import fiji.plugin.trackmate.TrackMate
import fiji.plugin.trackmate.SelectionModel
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory
import fiji.plugin.trackmate.SpotCollection
import fiji.plugin.trackmate.Spot
import fiji.plugin.trackmate.tracking.sparselap.SimpleSparseLAPTrackerFactory
import fiji.plugin.trackmate.SelectionModel
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer
import fiji.plugin.trackmate.FeatureModel
import fiji.plugin.trackmate.Spot

import loci.formats.ChannelSeparator
import loci.formats.meta.MetadataRetrieve
import loci.formats.services.OMEXMLServiceImpl
import loci.plugins.util.ImageProcessorReader
import loci.plugins.util.LociPrefs

import ome.units.UNITS

// create empty TrackMate Settings object
Settings settings = new Settings();
// use a Laplacian of Gaussian detector
settings.detectorFactory = new LogDetectorFactory();
// specify that inensity based features should be calculated for spots
settings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<>());
// specify the quality threshold for spots
settings.addSpotFilter(new FeatureFilter('QUALITY', spotThreshQual, true))
// specify the mean intensity threshold for spots
settings.addSpotFilter(new FeatureFilter('MEAN_INTENSITY', spotThreshMeanInt, true))
// create and fill map with detector settings
HashMap detectMap = settings.detectorFactory.getDefaultSettings();
detectMap.put('RADIUS', spotRadius);
detectMap.put('DO_MEDIAN_FILTERING', true);
detectMap.put('DO_SUBPIXEL_LOCALIZATION', false)
detectMap.put('THRESHOLD', 0.0d);
detectMap.put('TARGET_CHANNEL', 1);
settings.detectorSettings = detectMap; 
// start with the simple LAP tracker
settings.trackerFactory = new SimpleSparseLAPTrackerFactory();
// moddifcy with specified tracking settings
HashMap trackMap = settings.trackerFactory.getDefaultSettings();
trackMap.put('LINKING_MAX_DISTANCE', maxLinkDist);
trackMap.put('GAP_CLOSING_MAX_DISTANCE', maxLinkingDistGC);
trackMap.put('MAX_FRAME_GAP', maxFrameGap);
settings.trackerSettings = trackMap;

// find all files in specified directory with specified extension
String[] exts = new String[1];
exts[0] = ext;
List<File> files = FileUtils.listFiles(directory, exts, true);

// loop through all files
for (int f = 0; f < files.size(); f++) {

	// reader to load specified file
	ImageProcessorReader  r = new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
	// set up metadata reading
	OMEXMLServiceImpl OMEXMLService = new OMEXMLServiceImpl();
	r.setMetadataStore(OMEXMLService.createOMEXMLMetadata());
	// set reader to specified file
	r.setId(files.get(f).getAbsolutePath());
	// get metadata
	MetadataRetrieve meta = (MetadataRetrieve) r.getMetadataStore();
	// get number of series in specified file
	int numSeries = meta.getImageCount();
	
	// loop through all series in file
	for (int i = 0; i < numSeries; i++) {
		println ("Processing file " + (f + 1) + " of " + files.size() + ", series " + (i + 1) + " of " + numSeries)
		// set reader to current series
		r.setSeries(i);

		// only load series if more than 10 slices and more than 10 frames (avoids any snapshots which have been aquired and saved with the movies)
		if (r.getSizeZ() >= 10 && r.getSizeT() >= 10) { 

			// read data
			ImagePlus imp_raw = readSeries(r, meta, i);
			// crop to a specified size in Z and time so all movies are same size
			ImagePlus imp = new Duplicator().run(imp_raw, 1, r.getSizeC(), 1, z_crop, 1, t_crop);

			// track cells using TrackMate
			ResultsTable spotTable = trackCells(imp, settings, display); 
		
			// if csv file save is requested
			if (save) {
				// save csv results table
				spotTable.save(outDirectory.getAbsolutePath() + File.separator + files.get(f).getName().substring(0, files.get(f).getName().length() - ext.length() - 8) + 'spotStats.csv');
			}
				}
				
			}
	// close the reader
	r.close();
	

}

public ResultsTable trackCells(ImagePlus imp, Settings settings, boolean display) {

	settings.setFrom(imp);
	// create empty TrackMate Model
	Model model = new fiji.plugin.trackmate.Model();
	// create TrackMate object to peform spot detection and tracking
    TrackMate trackmate = new TrackMate(model, settings);
    // find spots
    ok = trackmate.execDetection();
    if (ok == false) {
        println(trackmate.getErrorMessage());
    }
   // compute spot featues (using first channel)
   ok = trackmate.computeSpotFeatures(true);
    if (ok == false) {
        println(trackmate.getErrorMessage());
    }
    // filter spots
   ok = trackmate.execSpotFiltering(true);
    if (ok == false) {
        println(trackmate.getErrorMessage());
    }
	// tracking
   ok = trackmate.execTracking();
    if (ok == false) {
        println(trackmate.getErrorMessage());
    }

	// export all spot features to a ResultsTable
	ResultsTable spotTable = calcAllSpotFeatures(trackmate);

	// if display is reuqested
		if (display) {
			// display spots and tracks on dataset
		    SelectionModel selectionModel = new SelectionModel(model);
			HyperStackDisplayer displayer =  new HyperStackDisplayer(model, selectionModel, imp);
			displayer.render();
			displayer.refresh();
			// show spot statistics
			spotTable.show("All Spots statistics");
			// wait for user to move on
			new WaitForUserDialog("ok?").show();
			//close window
			imp.close();	
		}
	
	return spotTable
}
public ImagePlus readSeries(ImageProcessorReader r, MetadataRetrieve meta, int series) {
	r.setSeries(series);
	// load data for current series as a stack
	ImageStack stack = new ImageStack(r.getSizeX(), r.getSizeY());
	for (int n = 0; n < r.getImageCount(); n++) {
		ImageProcessor ip = r.openProcessors(n)[0];
		stack.addSlice("" + (n + 1), ip);
	}
	// create ImagePlus using stack
	ImagePlus imp = new ImagePlus("", stack);
	// convert to HyperStack with correct dimensions
	imp = HyperStackConverter.toHyperStack(imp, r.getSizeC(), r.getSizeZ(), r.getSizeT());
	// set calibration for data using metadata
	Calibration cali = new Calibration();
	cali.pixelWidth = meta.getPixelsPhysicalSizeX(series).value(UNITS.MICROMETER).doubleValue();
	cali.pixelHeight = meta.getPixelsPhysicalSizeY(series).value(UNITS.MICROMETER).doubleValue();
	if (r.getSizeZ() > 1) {
		cali.pixelDepth = meta.getPixelsPhysicalSizeZ(series).value(UNITS.MICROMETER).doubleValue();
	}
	cali.setUnit("micron");
	imp.setGlobalCalibration(cali);
	// set imp title to series name
	imp.setTitle(meta.getImageName(series));
	
	return imp
}



public ResultsTable calcAllSpotFeatures(TrackMate trackmate) {
	
	// create table
	ResultsTable spotTable = new ResultsTable();
	FeatureModel fm = trackmate.getModel().getFeatureModel();
	// get list of all spot features
	Collection< String > spotFeatures = trackmate.getModel().getFeatureModel().getSpotFeatures();

	// iterate though all visible spots
	Iterable< Spot > iterable = trackmate.getModel().getSpots().iterable( true );
	for (Spot spot : iterable) {

			// add spot name and ID to table
			spotTable.incrementCounter();
			spotTable.addLabel( spot.getName() );
			spotTable.addValue( "ID", "" + spot.ID() );

			// c// loop through all series in fileheck if current spot is in a track and add Track ID to table
			Integer trackID = trackmate.getModel().getTrackModel().trackIDOf( spot );
			if ( null != trackID ) {
				spotTable.addValue( "TRACK_ID", "" + trackID.intValue() );
			}
			else {
				spotTable.addValue( "TRACK_ID", -1);
			}
			// add remaining spot features to table
			for (String feature : spotFeatures ) {
				Double val = spot.getFeature( feature );
				if ( null == val ) {
					spotTable.addValue( feature, "None" );
				}
				else {
					if ( fm.getSpotFeatureIsInt().get( feature ).booleanValue() ) {
						spotTable.addValue( feature, "" + val.intValue() );
					}
					else {
						spotTable.addValue( feature, val.doubleValue() );
					}
				}
			}
		}
	return spotTable;
}