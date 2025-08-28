//Lipid droplet counting pipeline
//Written by Laura Blackie, Miguel-Aliaga lab, The Francis Crick Institute, 2025

args = getArgument();
tokens = split(args, ",");
userChosenDirectory = tokens[0] + File.separator;
file = tokens[1];
saveFolder = "LipidDropletCount/";

folderPath = userChosenDirectory + saveFolder;

run("Bio-Formats Macro Extensions");
list = getFileList(folderPath);
numImages = 0;
width = height = 200;
FileTypeList = newArray("GutWatershed","GutMask");

for (j=0; j<2; j++){
	FileType = FileTypeList[j];
	numImages = 0;

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], FileType+".tif")) {
        open(folderPath + list[i]);
        name = getTitle();
        if (indexOf(name, " ") >= 0) {
        	name = replace(name, " ", "_");
        	rename(name);
        }
        savename = split(name, "(_Position)");
		savename = savename[1];
        
        run("Size...", "width="+width+" height="+height+" depth=1 constrain average interpolation=Bilinear");
        
        if (numImages == 0) {
            rename("Stack");
            setFont("SansSerif", 20);
            drawString(savename, 10, 25);

        } else {
            run("Concatenate...", "  title=Stack open image1=Stack image2="+name+"");
            setSlice(numImages+1);
            drawString(savename, 10, 25);
        }
        numImages++;
    }
}

if (numImages == 0) {
    showMessage("No images found matching "+FileType+"'_GutMask.tif'.");
} else {
    run("Make Montage...", "columns=" + floor(sqrt(numImages)) + " rows=" + (floor(numImages / floor(sqrt(numImages)))+1)+ " scale=1");
}

close("Stack");

run("8-bit");
saveAs("Tiff", folderPath+"QC_"+FileType+"_Montage.tif");
close("*");
}

FileType3 = "Spots"; 
numImages = 0;

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], FileType3+".tif")) {

		Ext.setId(folderPath + list[i]);
		Ext.getSizeZ(sizeZ);
		run("Bio-Formats Importer", "open=["+folderPath + list[i]+"] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1 z_begin="+sizeZ/2+" z_end="+sizeZ/2+" z_step=1");

        name = getTitle();
        if (indexOf(name, " ") >= 0) {
        	name = replace(name, " ", "_");
        	rename(name);
        }
        savename = split(name, "(_Position)");
		savename = savename[1];
        
        run("Size...", "width="+width+" height="+height+" depth=1 constrain average interpolation=Bilinear");
        
        if (numImages == 0) {
            rename("Stack");
            setFont("SansSerif", 20);
            drawString(savename, 10, 25);

        } else {
            run("Concatenate...", "  title=Stack open image1=Stack image2="+name+"");
            setSlice(numImages+1);
            drawString(savename, 10, 25);
        }
        numImages++;
    }
}

if (numImages == 0) {
    showMessage("No images found matching "+FileType3+"'_GutMask.tif'.");
} else {
    run("Make Montage...", "columns=" + floor(sqrt(numImages)) + " rows=" + (floor(numImages / floor(sqrt(numImages)))+1)+ " scale=1");
}

close("Stack");

run("8-bit");
saveAs("Tiff", folderPath+"QC_"+FileType3+"_Montage.tif");

close("*");
