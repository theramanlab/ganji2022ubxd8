/*
 *                                                  Mitochondria - ER Contact Points Macro
 *                                     
 *                         This macro was designed to work with ND2 image stacks. It may work with other image stacks
 *                         (using the "open()" function) however some parts may not function as expected. A new folder
 *                         will be created for each new image analyzed due to the large number of files generated per 
 *                         analysis. Analyzing the same image or an image of the same name in the same directory will
 *                         overwrite previous results. It is highly recommended to move folders with results to a new
 *                         location after finishing to prevent accidental data loss. 
 *                         
 *                         Plugins:
 *                         ND2Reader - https://imagej.nih.gov/ij/plugins/docs/ND2Reader.pdf
 *                         Colocalization Highlighter - Received from Rakesh Ganji
 *                         IsoPhotContour2 - http://www.mecourse.com/landinig/software/software.html
 */
//run("ND to Image6D");
open();
nd2path = File.directory+File.name;
nd2name = File.nameWithoutExtension;
savepath = File.directory+"\\"+nd2name+" analysis\\"; 
if(File.exists(savepath))
{
	if(getBoolean("File \""+nd2name+" analysis\"\nalready exists in that dir.\nOverwrite existing file?"))
	{
		File.makeDirectory(savepath);   //creates folder to save all analysis data
	}
	else
	{
		exit();
	}
}
else
{
	File.makeDirectory(savepath);   //creates folder to save all analysis data
}

//--------------------PRE-PROCESSING-----------------------------------
if(getBoolean("Begin pre-processing image?"))
{
	waitForUser("Select MITO channel"); //begin with the MITO image
	run("Duplicate...", "title=[mito image]"); //manipulating ND2 files with CLAHE does not work, so create a duplicate TIF file
	run("Subtract Background...", "Rolling = 50 disable");
	run("Despeckle");
	run("Enhance Local Contrast (CLAHE)", "blocksize=9 histogram=256 maximum=4 mask=*None* fast_(less_accurate)");
	run("Despeckle");
	run("Tubeness"); //tubeness creates a new image
	run("8-bit");
	mitotubename = "mito tubeness "+nd2name+".tif";
	save(savepath+mitotubename); //save the tubeness image to the save folder
	close();
	selectWindow("mito image");
	close();

	waitForUser("Select ER channel");
	run("Duplicate...", "title=[er image]"); //creating a new image isn't necessary but editing the original image may create irreversible changes
	run("Subtract Background...", "Rolling = 50 disable");
	run("Tubeness");
	run("8-bit");
	ertubename = "er tubeness "+nd2name+".tif";
	save(savepath+ertubename); //save the ER tubeness image to the save folder
	close();
	selectWindow("er image");
	close("\\Others");
	close();
	waitForUser("Pre-processing complete");
}
else
{
	close();
	waitForUser("Open mito tubeness image");
	mitotubenesspath = File.openDialog("Open pre-processed mito tubeness image");
	mitotubename = File.name;
	waitForUser("Open ER tubeness image");
	ertubenesspath = File.openDialog("Open pre-processed ER tubeness image");	
	ertubename = File.name;
}
//MAIN ANALYSIS MACRO
/*
 * In short, this macro takes two image and calculates the colocalization of one image (ER)
 * with the perimeter of the structures in the second image (mitochondria).
 * This protocol was adapted from Rakesh Ganji.
 * 
 * This macro is able to analyze multiple cells (or ROIs) in the same image. The ROIs will be numbered
 * starting from 1. The ROIs selected as well as the ROI restricted images will be saved if further analysis 
 * is needed.
 * 
 * This macro was written by Jacob Klickstein in the Raman lab at Tufts University (7/2017)
 */

if (isOpen("ROI Manager")) //prevents old ROIs from interfering by resetting the ROI manager
{
     selectWindow("ROI Manager");
     run("Close");
}

run("Clear Results"); //clear results in case there is old data
open(savepath+ertubename);

for(cellnum=1;getBoolean("Continue with analysis of "+nd2name+"?\nSelect 'No' to finish");cellnum++) //main for loop that iterates new cells in the same image
{

//ROI Processing step
waitForUser("Select ROI\nor click 'OK' to analyze\nentire image");
if(selectionType()==-1) {run("Select All");}
roiManager("Add");
roiManager("Show None");
roiManager("Select",cellnum-1);
run("Clear Outside");

//saves the ROI subsection of the ER image
ertubecell = replace(ertubename,".tif"," roi"+cellnum+".tif");
save(savepath+ertubecell);

getHistogram(erhisto,ercounts,256);
erweightedmean = ercountsum = 0;
for(k=10,ercountsum=0;k<ercounts.length;k++)
{
	ercountsum = ercountsum+ercounts[k];
}
for(k=10;k<ercounts.length;k++)
{
	erweightedmean = erweightedmean+(k*(ercounts[k]/ercountsum));
}

open(savepath+mitotubename); //open the er tubeness fill from before
roiManager("Select", cellnum-1); //use the same ROI as selected above
run("Clear Outside");


mitotubecell = replace(mitotubename,".tif"," roi"+cellnum+".tif");
save(savepath+mitotubecell); //save the Mito tubeness image of the ROI

getHistogram(mitohisto,mitocounts,256); //This returns an array (0-255) of the histogram data
mitoweightedmean = mitocountsum = 0;
for(k=10,mitocountsum=0;k<mitocounts.length;k++) //the first ten bins are considered background noise
{
	mitocountsum = mitocountsum+mitocounts[k]; //this calculates the denominator of the weighted average
}
for(k=10;k<mitocounts.length;k++)
{
	mitoweightedmean = mitoweightedmean+(k*(mitocounts[k]/mitocountsum)); //this calculates the weighted average
}

mitoweightedmean75 = 0.75*mitoweightedmean; //we will use 75% of the weighted mean intensity of the mitochondria

/*colocalization highlighter will ask for which image is which color with
 * the options being green or red . As far as I can tell, it does not matter
 * which is which as they have been converted to 8-bit images earlier. If you want
 * the merge image colors to match these colors, you may want to change this around
 * or do this comman manually.
 * This step is where the weightedmeans are used
 * This is also where we will get our total colocalization data from
 */
run("Colocalization Highligter", "channel_1=["+mitotubename+"] channel_2=["+ertubename+"] ratio=75 threshold_channel_1="+mitoweightedmean75+" threshold_channel_2="+erweightedmean+" display=255 colocalized");
colocalizename = "colocalized points "+nd2name+".tif";
windowfound=0;
while(windowfound==0) //Cycle through the windows, closing them until we get to the 8-bit image, then save it
{
	if(matches(getTitle(),".*8.bit.*"))
	{
		save(savepath+colocalizename);
		windowfound=1;
	}
	close();
}
open(savepath+mitotubecell); //open the mito tubeness ROI image created earlier

setAutoThreshold("Default"); //threshold using the 75% weighted mean
//call("ij.plugin.frame.ThresholdAdjuster.setMode", "B&W");
setThreshold(mitoweightedmean75, 255);
setAutoThreshold("Default dark");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Invert");

waitForUser("Check background");

run("IsoPhotContour2 ", "contours"); //creates a new IsoPhot image
run("8-bit");
run("Invert");
//run("Mean...", "radius=0"); for thicker borders, uncomment this command
setAutoThreshold("Default dark"); //threshold again, this time creating a binary image
setThreshold(1, 255);
setOption("BlackBackground", false);
run("Convert to Mask");

/*calcs the histo of the IsoPhot but there should be only two values for the pixels,
 * either 256 (mito perimeter) or 1 (background). If there are more values than this,
 * then the thresholding from the step before did not work properly
 */

getHistogram(mitohisto,mitocounts,256);

open(savepath+colocalizename); //opens the colocalized image created earler

//This AND function finds the points of colocalization that are also the perimeter
imageCalculator("AND create", colocalizename,"IsoPhot");

getHistogram(perihisto,pericounts,256);
ratio = pericounts[255]/mitocounts[255];
close("*tubeness*");

/*Looks for windows of interest and saves them. The IF statements prevent the macro from crashing and data loss if the user has closed one
 * of the windows earlier, or if the macro was unable to perform one of the stteps
 */
if(isOpen("Result of colocalized points "+nd2name+".tif"))
{
	selectWindow("Result of colocalized points "+nd2name+".tif");
	save(savepath+nd2name+"_perim_cont_points_roi"+cellnum+".tif");
}
if(isOpen("IsoPhot"))
{
	selectWindow("IsoPhot");
	save(savepath+nd2name+"_mito_perim_roi"+cellnum+".tif");
}

//outputs all the results to the results window which will be saved every iteration for these images
setResult("ROI num",cellnum-1,cellnum);
setResult("ER W Mean",cellnum-1,erweightedmean);
setResult("Mito W Mean",cellnum-1,mitoweightedmean);
setResult("Mito Peri",cellnum-1,mitocounts[255]);
setResult("Peri coloc",cellnum-1,pericounts[255]);
setResult("Ratio",cellnum-1,ratio);
updateResults();
selectWindow("Results");
saveAs("Results",savepath+nd2name+" results.csv");
roiManager("Deselect");
roiManager("save",savepath+nd2name+"_rois.zip");

waitForUser("Click 'okay' to continue");
//Now we save our images and clean up some of the windows
selectWindow("colocalized points "+nd2name+".tif");
close();
selectWindow("IsoPhot");
close();
selectWindow("Result of colocalized points "+nd2name+".tif");
close();
open(savepath+ertubename);
roiManager("Show All");
}
run("Clear Results");
if (isOpen(ertubename))
{
     selectWindow(ertubename);
     run("Close");
}
if (isOpen("ROI Manager")) //prevents old ROIs from interfering by resetting the ROI manager
{
     selectWindow("ROI Manager");
     run("Close");
}
