//Lipid droplet counting pipeline
//Written by Laura Blackie, Miguel-Aliaga lab, The Francis Crick Institute, 2025

//macroPath = getDirectory("Choose a Directory");
macroPath = "path\\to\\Macro_pipeline\\folder\\";

// Create a selection dialog with checkboxes and a file input field
Dialog.create("Select Processing Steps");
Dialog.addFile("Select Image:", "");
Dialog.addCheckbox("Step 1: Threshold Selection", true);
Dialog.addCheckbox("Step 2: Draw centrelines", true);
Dialog.addCheckbox("Step 3: Find spots and calculate measurements", true);
Dialog.addCheckbox("Step 4: QC", true);
Dialog.addCheckbox("Step 5: Calculate intensity along line", true);
//Dialog.addCheckbox("Step 6: Count nuclei", true);
Dialog.show();

// Get user inputs
imagePath = Dialog.getString();
macroChoice = Dialog.getCheckbox();
macro2Choice = Dialog.getCheckbox();
macro3Choice = Dialog.getCheckbox();
macro4Choice = Dialog.getCheckbox();
macro5Choice = Dialog.getCheckbox();
//macro6Choice = Dialog.getCheckbox();

imageDir = File.getParent(imagePath);
imageName = File.getName(imagePath);

//Create processing folders
Dir1 = imageDir+File.separator+"LipidDropletCount"+File.separator;
if (!File.exists(Dir1)){
	File.makeDirectory(Dir1);
	}
Dir2 = imageDir+File.separator+"IntensityAlongLength"+File.separator;
if (!File.exists(Dir2)){
	File.makeDirectory(Dir2);
	}

if (macroChoice) {
    runMacro(macroPath+"LDpipeline_findParameters.ijm", imageDir + "," + imageName + "," + macroPath);
}

if (macro2Choice) {
    runMacro(macroPath+"LDpipeline_generateCentrelines.ijm", imageDir + "," + imageName);
}
if (macro3Choice) {
	thresholdFile = imageDir+File.separator+"selected_threshold.txt";
	threshold = File.openAsString(thresholdFile);
    runMacro(macroPath+"LDpipeline_3DspotAnalysis.ijm", imageDir + "," + imageName + "," + threshold);
}
if (macro4Choice) {
    runMacro(macroPath+"LDpipeline_QC.ijm", imageDir + "," + imageName);
}
if (macro5Choice) {
    runMacro(macroPath+"LDpipeline_intensityAlongGutLength.ijm", imageDir + "," + imageName);
}
//if (macro6Choice) {
//    runMacro(macroPath+"Macro_openCentrelines_3DDAPIAlongGutLength_threshold_pipelineVer2.ijm", imageDir + "," + imageName);
//}

print("Processing complete!");
