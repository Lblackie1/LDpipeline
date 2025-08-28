//Lipid droplet counting pipeline
//Written by Laura Blackie, Miguel-Aliaga lab, The Francis Crick Institute, 2025

//Currently it is not possible to use two sets of macro extensions (bio-formats and 3D manager) one after the other in a loop 
//so save the images from the lif into tif and then reopen for 3D manager
args = getArgument();
tokens = split(args, ",");
userChosenDirectory = tokens[0] + File.separator;
file = tokens[1];
lowthresh = tokens[2];
saveFolder = "LipidDropletCount/";

run("Bio-Formats Macro Extensions");
Ext.setId(userChosenDirectory + file);
Ext.getSeriesCount(seriesCount);

for (series = 0; series <=seriesCount-1 ; series++) {
	//Following commands open only lightning images inside the lif file by name:
	Ext.setSeries(series);
	Ext.getSeriesName(seriesName);
	if(indexOf(seriesName, "Lng")<=0){
	}
	else {
		run("Bio-Formats Importer", "open=[" + userChosenDirectory + file + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+(series+1));
		name = getTitle();
		savename = split(seriesName, "(\/)");
		savename = file+savename[2];
		//savename = savename[2];
		//savename = seriesName;
		//savename = file+savename;
		saveAs("Tiff", userChosenDirectory+saveFolder+savename+"_original.tif");
		rename(name);

		//get mask of gut 2d
		getDimensions(width, height, channels, slices, frames);
		run("Z Project...", "projection=[Max Intensity]");
		projection = getTitle();
		Stack.setChannel(1);
		run("Duplicate...", " ");
		rename("projection_ch2");
		run("8-bit");
				
		selectImage(projection);
		Stack.setChannel(2);
		run("Duplicate...", " ");
		rename("projection_lpd");
		imageCalculator("Min create", "projection_ch2","projection_lpd");
		rename("projection_min");
		resetMinAndMax();
		setAutoThreshold("Li dark no-reset");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Dilate");
		run("Dilate");
		run("Erode");
		run("Erode");
		run("Analyze Particles...", "size=5000-Infinity pixel show=Masks clear");
		run("Convert to Mask");
			
		rename("GutMask");
		saveAs("Tiff", userChosenDirectory+saveFolder+savename+"_GutMask.tif");
		rename("GutMask.tif");
		close(projection);
		close("* projection_ch2"); 
		close("projection_ch2_*");
		close("projection_lpd");
		close("*projection_min*");
				
		//optional crop step for large images
		//3D object counter doesn't work for images above a certain size so we can crop to gut size
		selectImage(name);
		getDimensions(width, height, channels, slices, frames);
		if (width*height*slices >= 2147483647){
			selectWindow("GutMask.tif");
			run("Create Selection");
			getSelectionBounds(x, y, width, height);
			makeRectangle(x, y, width, height);
			run("Crop");
			selectWindow(name);
			makeRectangle(x, y, width, height);
			run("Crop");
			run("Stack to Hyperstack...", "order=xyczt(default) channels=2 slices="+slices+" frames=1 display=Color");
			saveAs("Tiff", userChosenDirectory+saveFolder+savename+"_original_crop.tif");
			rename(name);
			selectWindow("projection_ch2");
			makeRectangle(x, y, width, height);
			run("Crop");
		}
	
		//process original image and watershed particles by distance transfrom watershed
		selectWindow(name);
		run("Duplicate...", "duplicate channels=2");
		rename("OGchannel2");
		close(name);
		selectWindow("OGchannel2");
		run("Gaussian Blur 3D...", "x=1 y=1 z=1");
		setThreshold(lowthresh, 65535);
		run("Convert to Mask", "background=Dark black");
		run("Distance Transform Watershed 3D", "distances=[Borgefors (3,4,5)] output=[16 bits] normalize dynamic=0.5 connectivity=6");
		saveAs("Tiff", userChosenDirectory+saveFolder+savename+"_processedDTlabelledSpots.tif");
		rename("labelledSpots");
		close("labelledSpots");
		close("OGchannel2");
			
		//open user drawn centreline and convert to points
		selectWindow("projection_ch2");
		roiManager("open", userChosenDirectory + saveFolder + savename + "_centreline.roi");
		roiManager("Select", 0);
		getSelectionCoordinates(xp, yp);
		List.setMeasurements();
		total = List.getValue("Length");
		steps = 20;
		xm=newArray(steps);
		ym=newArray(steps);
		for (i=0;i<xm.length;i++) {
			xm[i]=xp[i*xp.length/steps];
			ym[i]=yp[i*xp.length/steps];
		}
		makeSelection("Point medium dot white", xm, ym);
		roiManager("deselect");
		roiManager("delete");
		roiManager("Add");
			
		getDimensions(width, height, channels, slices, frames);
		newImage("Untitled", "8-bit black", width, height, 1);
		roiManager("Select", 0);
		run("Flatten");
		run("Select None");
		run("Convert to Mask");
		run("Convert to Mask");
		run("Analyze Particles...", "size=150-Infinity show=Masks clear");
		run("Convert to Mask");
		rename("points");
		saveAs("Tiff", userChosenDirectory+saveFolder+savename+"_points.tif");
		rename("points.tif");
		close("Untitled*");
		close("projection_ch2");
			
		//for watershed to work, use gut mask as input image, as mask and use binary seeds
		run("Marker-controlled Watershed", "input=GutMask.tif marker=[points.tif] mask=GutMask.tif compactness=1 binary calculate use");
		run("Morphological Filters", "operation=Closing element=Square radius=20");
		saveAs("Tiff", userChosenDirectory+saveFolder+savename+"_GutWatershed.tif");
		rename("GutWatershed.tif");
		close("GutMask.tif");
		close("points.tif");
		close("GutMask-watershed.tif");
		selectWindow("GutWatershed.tif");
		roiManager("Select", 0);
		roiManager("Measure");
		selectWindow("Results");
		saveAs("Results", userChosenDirectory+saveFolder+savename+"_GutRegionValues.csv");
		close("Results");
		roiManager("deselect");
		roiManager("delete");
		roiManager("open", userChosenDirectory +saveFolder+ savename + "_centreline.roi");
		run("Clear Results");
		roiManager("Select", 0);
		profile = getProfile();
		for (i=0; i<profile.length; i++){
  			setResult("Value", i, profile[i]);
			updateResults();}
		saveAs("Measurements", userChosenDirectory +saveFolder+savename+ "_gutMaskWatershedProfile.txt");
		roiManager("deselect");
		roiManager("delete");
		roiManager("reset");
	
		//get R3 value
		allText = File.openAsString(userChosenDirectory+saveFolder+savename+"_R3coords.csv");
		text = split(allText, "\n");
		text = split(text[1], ",");
		x=text[1];
		y=text[2];
		getVoxelSize(width, height, depth, unit);
		makePoint(x/width,y/width);
		run("Set Measurements...", "modal display redirect=None decimal=3");
		run("Measure");
		saveAs("Results", userChosenDirectory+saveFolder+savename+"_R3MaskVal.csv");
		close("Results");
	
	
		close("GutWatershed.tif");
		close("*");
	}
}

fileList2 = getFileList(userChosenDirectory+saveFolder); //list of processed files
for (file2 = 0; file2 < fileList2.length; file2++) {
	if(endsWith(fileList2[file2], "original.tif")) {
		savename = split(fileList2[file2], "(\_o)");
		savename = savename[0];
		
		open(userChosenDirectory+saveFolder+savename+"_processedDTlabelledSpots.tif");
		rename("labelledSpots");
		
		// Check if the image is all black
    	selectImage("labelledSpots");
    	getDimensions(width, height, channels, slices, frames);
    	setSlice(slices/2);
    	getStatistics(area, mean, min, max, stdDev);
   	 	if (min == 0 && max == 0) {
      	  print(savename+" contains all black pixels.");
      	  close("labelledSpots");
      	  continue;
    	} 
			
		open(userChosenDirectory+saveFolder+savename+"_original.tif");
		rename("original");	
		open(userChosenDirectory+saveFolder+savename+"_GutWatershed.tif");
		rename("GutWatershed.tif");	
		
		//make watershed 3d
		selectImage("original");
		getDimensions(width, height, channels, slices, frames);
		selectImage("GutWatershed.tif");
		run("16-bit");
		run("Concatenate...", "keep open image1=GutWatershed.tif image2=GutWatershed.tif image3=[-- None --]");
		run("Scale...", "x=1.0 y=1.0 z="+slices/2+" width=355 height=363 depth=78 interpolation=None average process create");
		rename("GutWatershed_z");
		close("Untitled*");
		close("GutWatershed.tif");
			
		selectWindow("original");
		if (width*height*slices >= 2147483647){
			close("original");
			open(userChosenDirectory+saveFolder+savename+"_original_crop.tif");
			rename("original");
		}

		//find spots and positions along gut length
		run("3D Manager Options", "volume surface compactness integrated_density mean_grey_value std_dev_grey_value mode_grey_value minimum_grey_value maximum_grey_value centroid_(pix) distance_between_centers=10 distance_max_contact=1.80 drawing=Contour");
		run("3D Manager");
		selectImage("labelledSpots");
		Ext.Manager3D_AddImage();
		Ext.Manager3D_Measure();
		Ext.Manager3D_SaveResult("M", userChosenDirectory+saveFolder+savename+"_DTspots_measurement_Results.csv");
		Ext.Manager3D_CloseResult("M");
		selectImage("original");
		Stack.setChannel(2);
		Ext.Manager3D_Quantif();
		Ext.Manager3D_SaveResult("Q", userChosenDirectory+saveFolder+savename+"_DTspots_originalImage_quantification_Results.csv");
		Ext.Manager3D_CloseResult("Q");
		selectImage("GutWatershed_z");
		Ext.Manager3D_Quantif();
		Ext.Manager3D_SaveResult("Q", userChosenDirectory+saveFolder+savename+"_DTspots_GutMaskMode_Results.csv");
		Ext.Manager3D_CloseResult("Q");
		Ext.Manager3D_Delete();
		Ext.Manager3D_Close();
		close("*");
	}
}
selectWindow("Log");
saveAs("Text", userChosenDirectory+saveFolder+"skippedFileLog.txt");
