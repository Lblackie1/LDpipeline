//Lipid droplet counting pipeline
//Written by Laura Blackie, Miguel-Aliaga lab, The Francis Crick Institute, 2025

args = getArgument();
tokens = split(args, ",");
userChosenDirectory = tokens[0] + File.separator;
file = tokens[1];
saveFolder = "IntensityAlongLength/";

run("Bio-Formats Macro Extensions");
Ext.setId(userChosenDirectory + file);
Ext.getSeriesCount(seriesCount);
for (series = 0; series <=seriesCount-1 ; series++) {
	//Following commands allow you to open only lightning images inside the lif file by name:
	Ext.setSeries(series);
	Ext.getSeriesName(seriesName);
	if(indexOf(seriesName, "Lng")<=0){
	}
	else {
	run("Bio-Formats Importer", "open=[" + userChosenDirectory + file + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+(series+1));
	name = getTitle();
	//savename = file+seriesName;
	savename = split(seriesName, "(\/)");
	savename = file+savename[2];
	//savename = seriesName;
	
	run("Z Project...", "projection=[Max Intensity]");
	projection = getTitle();
	Stack.setChannel(1);
	run("Duplicate...", " ");
	rename("projection_ch2");
	
	//optional crop step for large images
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
		close("* projection_ch2"); 
		close("projection_ch2_*");
		close("projection_lpd");
		close("projection_min");
				
		selectWindow("GutMask");
		run("Create Selection");
		getSelectionBounds(a, b, width2, height2);
		selectWindow(projection);
		makeRectangle(a, b, width2, height2);
		run("Crop");
		close("GutMask");
		}
	close("projection_ch2");
	
	selectWindow(name);
	close();
	selectWindow(projection);
	run("Split Channels");
	selectWindow("C1-"+projection);
	close();
	selectWindow("C2-"+projection);
	
	//open user drawn centreline
	roiManager("open", userChosenDirectory + saveFolder + savename + "_centreline.roi");
	roiManager("Select", 0);

	getVoxelSize(width, height, depth, unit);
	lineWidth = 113.55/(round(width*1000)/1000);
	roiManager("Set Line Width", lineWidth);
	run("Plot Profile");
	plotname = getTitle();
  	Plot.getValues(x, y);	
	for (j=0; j<x.length; j++){
   		  setResult("Position", j, x[j]);
   		  setResult("Intensity", j, y[j]);
   		  updateResults();}
	saveAs("Measurements", userChosenDirectory + saveFolder  + savename + "_IntensityAlongCentreline.csv");
	selectWindow(plotname);
	run("Close");

  	selectWindow("C2-"+projection);
  	roiManager("Select", 0);
  	getSelectionCoordinates(x, y);
 	getPixelSize(unit, pixelWidth, pixelHeight);
 	run("Clear Results");
  	count = 0;
  	n = x.length;
  	for (k=0; k<n; k++) {
      	setResult("X", k, x[k]);
      	setResult("Y", k, y[k]);
      	dx = (x[(k+1)%n] - x[k])*pixelWidth;
      	dy = (y[(k+1)%n] - y[k])*pixelHeight;
      	length = sqrt(dx*dx+dy*dy);
      	if (k<n-1 || selectionType<=4) {
        		 setResult("Length", k, length); 
          	sumLength += length;
          	count++;
      	}      
  	}
  	selectWindow("Results"); 
	saveAs("Results", userChosenDirectory + saveFolder  + savename + "_XYCoordinatesOfLengthandDistances.csv");
	run("Close");
	selectWindow("C2-"+projection);
	run("Select None");
	close();
	roiManager("Deselect");
	roiManager("Delete");
	}
}
		