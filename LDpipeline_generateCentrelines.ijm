//Lipid droplet counting pipeline
//Written by Laura Blackie, Miguel-Aliaga lab, The Francis Crick Institute, 2025

args = getArgument();
tokens = split(args, ",");
userChosenDirectory = tokens[0] + File.separator;
file = tokens[1];
saveFolder = "LipidDropletCount/";
saveFolder2 = "IntensityAlongLength/"

run("Bio-Formats Macro Extensions");
run("Set Measurements...", "  redirect=None decimal=3");

Ext.setId(userChosenDirectory + file);
Ext.getSeriesCount(seriesCount);

print(seriesCount);

for (series = 0; series <=seriesCount-1 ; series++) {
	//Following commands allow you to open only lightning images inside the lif file by name:
	Ext.setSeries(series);
	Ext.getSeriesName(seriesName);
	if(indexOf(seriesName, "Lng")<=0){
	}
	else {
		run("Bio-Formats Importer", "open=[" + userChosenDirectory + file + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+(series+1));
		name = getTitle();
		savename = split(seriesName, "(\/)");
		savename = file+savename[2];
		//savename = seriesName;

		run("Z Project...", "projection=[Max Intensity]");
		projection = getTitle();
		Stack.setChannel(1);
		run("Duplicate...", " ");
		rename("projection_ch2");

		//optional crop step for large images
		//3D object counter doesn't work for images above a certain size so we can crop to gut size
		selectImage(name);
		getDimensions(width, height, channels, slices, frames);
		if (width*height*slices >= 2147483647){
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
			close(projection);
			close("* projection_ch2");
			close("projection_ch2_*");
			close("projection_lpd");				
			close("projection_min");

			selectWindow("GutMask");
			run("Create Selection");
			getSelectionBounds(x, y, width, height);
			makeRectangle(x, y, width, height);
			run("Crop");
			selectWindow("projection_ch2");
			makeRectangle(x, y, width, height);
			run("Crop");
			}

		//get user drawn centreline and R3 coodinate
		selectWindow("projection_ch2");
		run("Enhance Contrast", "saturated=25");
		setTool("freeline");
		waitForUser("Before pressing OK, draw a freehand line along centreline of gut, then press OK");
		run("Fit Spline");
		run("Interpolate", "interval=1 smooth");
		roiManager("Add");
		roiManager("Save", userChosenDirectory + saveFolder+ savename + "_centreline.roi");
		roiManager("Save", userChosenDirectory + saveFolder2+ savename + "_centreline.roi");
		roiManager("delete");
		roiManager("reset");

		setTool("multipoint");
		waitForUser("Before pressing OK, place a point at R3, then press OK");
		run("Measure");
		saveAs("Results", userChosenDirectory + saveFolder+ savename + "_R3coords.csv");
		saveAs("Results", userChosenDirectory + saveFolder2+ savename + "_R3coords.csv");
		close("Results");

		close("*");
	}
}
