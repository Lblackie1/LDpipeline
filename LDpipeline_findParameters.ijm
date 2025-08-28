//Lipid droplet counting pipeline
//Written by Laura Blackie, Miguel-Aliaga lab, The Francis Crick Institute, 2025

//Macro to set parameters for lipid counting for each new experiment
args = getArgument();
tokens = split(args, ",");
userChosenDirectory = tokens[0] + File.separator; 
file = tokens[1];
macroPath = tokens[2];

run("Bio-Formats Macro Extensions");
Ext.setId(userChosenDirectory + file);
Ext.getSeriesCount(seriesCount);
testno = floor(seriesCount/4*3);
Ext.setSeries(testno);
Ext.getSeriesName(seriesName);

Dialog.create("Choose Test Image");
Dialog.addRadioButtonGroup("Use default image "+testno+"/"+seriesCount+" : "+seriesName+"?", newArray("Yes", "No"), 1, 2, "Yes");
Dialog.addNumber("Enter test number:", testno);
Dialog.addNumber("Enter initial threshold value:", 10200);

Dialog.show();

useDefault = Dialog.getRadioButton();
testNumber = Dialog.getNumber();
lowThresh = Dialog.getNumber();

if (useDefault == "No") {
    testno = testNumber;
}

Ext.setSeries(testno);
Ext.getSeriesName(seriesName);
print(seriesName);
if(indexOf(seriesName, "Lng")<=0){
	print("Not lightning version");
}
else {
	run("Bio-Formats Importer", "open=[" + userChosenDirectory + file + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+(testno+1));

	name = getTitle();
	print(name);
	//savename = split(seriesName, "(\/)");
	//savename = file+savename[2];
	savename = seriesName;
	savename = file+savename;

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
	
	selectImage(name);
	run("Duplicate...", "duplicate channels=2");
	rename("OGchannel2");
	run("Gaussian Blur 3D...", "x=1 y=1 z=1");
	
	doneAdjusting = false;
	lowThresh = 10200;
	slices = 75;
	while (!doneAdjusting) {
		selectImage("OGchannel2");
		run("Duplicate...", "duplicate ");
		rename("Thresholded");
		setThreshold(lowThresh, 65535);
		run("Convert to Mask", "background=Dark black");
		
		selectImage("OGchannel2");
		setSlice(slices/2);
		selectImage("Thresholded");
		setSlice(slices/2);
	
		//user review the images
		waitForUser("Review the images. Use Image>Adjust>Threshold on OGchannel2 to pick a new threshold if needed.\nIf a new threshold value is better, make a note of it to input in next dialog box.\nAlso check gutMask but if this needs adjustment will have to change code directly.\n Once finished checking threshold, click OK."); 
		
		Dialog.create("Threshold Adjustment");
    	Dialog.addMessage("After image review, was new threshold needed?\nIf so, add new threshold value below. If not, leave as default.");
    	Dialog.addRadioButtonGroup("Happy with current threshold?", newArray("Yes", "No"), 1, 2, "Yes");
    	Dialog.addNumber("New threshold if needed:", lowThresh);

   	 	Dialog.show();

    	lowThresh = Dialog.getNumber();
    	userChoice = Dialog.getRadioButton();
    
    	if (userChoice == "Yes") {
        	doneAdjusting = true; // Exit loop
    	} else {
    		selectImage("OGchannel2");
			close("\\Others");
    	}
	}

	thresholdFile = userChosenDirectory + "selected_threshold.txt";

	File.saveString(lowThresh, thresholdFile);
	close('*');
}
	
	
	
	