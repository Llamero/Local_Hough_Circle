package Hough_Package;

/** Hough_Circle.java:
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 @author Ben Smith (benjamin.smith@berkeley.edu)
 @Based on original plugin implementation by Hemerson Pistori (pistori@ec.ucdb.br) and Eduardo Rocha Costa
 @created February 4, 2017
 
 The Hough Transform implementation was based on 
 Mark A. Schulze applet (http://www.markschulze.net/)
 
*/

import ij.*;
import ij.process.*;
import java.awt.*;
import ij.plugin.HyperStackConverter;
import ij.plugin.frame.*;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.*;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import java.util.Arrays;

/**
 *
 * @author Ben
 */
public class Hough_Circle extends SwingWorker<Integer, String>{
    
    private int radiusMin;  // Find circles with radius grater or equal radiusMin - argument syntax: "min=#"
    private int radiusMax;  // Find circles with radius less or equal radiusMax - argument syntax: "max=#"
    private int radiusInc;  // Increment used to go from radiusMin to radiusMax - argument syntax: "inc=#"
    private int minCircles; // Minumum number of circles to be found - argument syntax: "minCircles=#"    
    private int maxCircles; // Maximum number of circles to be found - argument syntax: "maxCircles=#"
    private int threshold = -1; // An alternative to maxCircles. All circles with a value in the hough space greater then threshold are marked. Higher thresholds result in fewer circles. - argument syntax: "threshold=#"
    private double thresholdRatio; //Ratio input from GUI that expresses threshold as ratio of resolution (highest possible # of votes)
    private int resolution; //The number of steps to use per transform (i.e. number of voting rounds)
    private double ratio; // Ratio of found circle radius to clear out surrounding neighbors
    private int searchBand = 0; //The +/- range of radii to search for relative to the last found radius - argument syntax: "bandwidth=#"
    private int searchRadius = 0; //The search radius to look for the next centroid relative to the last found centroid - argument syntax: "radius=#"
    private boolean reduce = false; //Cap the transform resolution by removeing redundant steps
    private boolean local = false; //Whether or not the search is going to be local
    
    //Output parameters
    private boolean houghSeries = false; //Contains whether the user wants the Hough series stack as an output - argument syntax: "show_raw"
    private boolean showCircles = false; //Contains whether the user wants the circles found as an output - argument syntax: "show_mask"
    private boolean showRadius = false; //Contains whether the user wants a map of centroids and radii outputed from search - argument syntax: "show_centroids"
    private boolean showScores = false; //Contains whether the user wants a map of centroids and Hough scores outputed from search - argument syntax: "show_scores"
    private boolean results = false; //Contains whether the user wants to export the measurements to a reuslts table 
            
    //Hough transform variables
    public ImagePlus imp; //Initalize the variable to hold the image
    private boolean isStack = false; //True if there is more than one slice in the input data
    private int stackSlices; //number of slices in the stack
    private int maxHough; //Contains the brights pixel in the entire Hough array
    private Point maxPoint; //Stores the location of the brightest pixel in the Hough array
    private int maxRadius; //Stores the radius of the brightest pixel in the Hough array
    private Rectangle r; //Stores the ROI on the original image
    private float imageValues[]; // Raw image (returned by ip.getPixels()) - float is used to allow 8, 16 or 32 bit images
    private int houghValues[][][]; // Hough Space Values [X coord][Y coord][radius index]
    private int localHoughValues[][][][]; //Local Hough space [circle#][X coord][Y coord][radius index]
    private int localHoughParameters[][]; //Array to pass local Hough space parameters to centroid search [circle#][parameter vector]
    private int width; // ROI width
    private int height;  // ROI height
    private int depth;  // Number of slices
    private int fullWidth; // Image Width
    private int fullHeight; //Image Height
    private int offx;   // ROI x origin position
    private int offy;   // ROI y origin position
    private Point centerPoint[]; // Center Points of the Circles Found.
    private int centerRadii[]; //Corresponding radii of the cricles marked by the center points
    private int houghScores[]; //Corresponding Hough scores for each centroid
    private int circleID[]; //Corresponding ID # for each centroid
    private int circleIDcounter; //Counter for keeping track of current max ID #
    private int lut[][][]; // LookUp Table for x and y tranform shifts in an octahedral manner
    private int lutSize; //Stores the actual number of transforms performed (<=selected resolution)
    private int nCircles; //Number of circles found during search - <= maxCircles
    private int nCirlcesPrev; //Stores nCirlces from last iteration
    private boolean localSearch = false; //Record whether a local-only search was done for the frame
    
    //Variables for storing the results and exporting result images
    public ImagePlus houghPlus;
    public ImageStack houghStack;
    public ImagePlus circlePlus;
    public ImageStack circleStack;
    public ImagePlus radiusPlus;
    public ImageStack radiusStack;
    public ImagePlus scorePlus;
    public ImageStack scoreStack;
    public ImageProcessor circlesip;
    public ImageProcessor radiusip;
    public ImageProcessor scoresip;
    public ResultsTable rt; 
    private String method;
    
    //Variables for max Hough score search
    private int maxHoughArray[][]; //Matrix to store hough scores, radii, and points from multi-threaded max search
    private int ithread;
    private int totalTime = 0; //Variable to test beenfits of multithreading
    // </editor-fold>
    
   

    @Override
     protected Integer doInBackground() throws Exception{
        Hough_GUI parameters = new Hough_GUI();
        this.radiusMin = parameters.getParam_radiusMin();
        this.radiusMax = parameters.getParam_radiusMax();
        this.radiusInc = parameters.getParam_radiusInc();
        this.minCircles = parameters.getParam_minCircles();
        this.maxCircles = parameters.getParam_maxCircles();
        this.threshold = parameters.getParam_threshold();
        this.thresholdRatio = parameters.getParam_thresholdRatio();
        this.resolution = parameters.getParam_resolution();
        this.ratio = parameters.getParam_ratio();
        this.searchBand = parameters.getParam_searchBand();
        this.searchRadius = parameters.getParam_searchRadius();
        this.reduce = parameters.getParam_reduce();
        this.local = parameters.getParam_local();
        this.houghSeries = parameters.getParam_houghSeries();
        this.showCircles = parameters.getParam_showCircles();
        this.showRadius = parameters.getParam_showRadius();
        this.showScores = parameters.getParam_showScores();
        this.results = parameters.getParam_results(); 
IJ.log("I am in the analysis class.");

        startTransform();
        return 0;
    }
    
    public void startTransform(){
        //Initialize the results table
        rt = Analyzer.getResultsTable();
        if (rt == null) {
                rt = new ResultsTable();
                Analyzer.setResultsTable(rt);
        }
        
        //Initialize variables
        nCircles = 0;
        circleIDcounter = 0;
IJ.log("I am doing a transform." + radiusMax);
        //Calculate Hough parameters
        this.depth = ((this.radiusMax-this.radiusMin)/this.radiusInc)+1;
        
        //If radiusInc is not a divisor of radiusMax-radiusMin, return error
        if((radiusMax-radiusMin)%radiusInc != 0){
            IJ.showMessage("Radius increment must be a divisor of maximum radius - minimum radius.");
            IJ.showMessage("radiusMin=" + radiusMin + ", radiusMax=" + radiusMax + ", radiusInc=" + radiusInc);
            return;          
        }
        
        //Build the transform LUT (all necessary translations in Cartesian coordinates)
        //NOTE: This step must precede the calculatation of the threshold, as it can change the resolution from the original input value
IJ.log("Hi");  
        lutSize = buildLookUpTable();
IJ.log("Hi2");          
        //Calculate the threshold based off the the actual resolution
        threshold = (int) Math.round(thresholdRatio * resolution);
        
        //If the threshold rounds to 0, set threshold to 1 as 0 is nonsensical
        if(threshold <= 0) {
            threshold = 1;
        }        
        // <editor-fold desc="Send arguments to record">
        //If the macro was being recorded, return the set values to the recorder
        if (Recorder.record){
            String Command  = "run(\"Local Hough Circle Transform\",\"minRadius=" + radiusMin + ", maxRadius=" + radiusMax + ", inc=" + radiusInc + 
                    ", minCircles=" + minCircles + ", maxCircles=" + maxCircles + ", threshold=" + thresholdRatio + ", resolution=" + resolution +
                     ", ratio=" + ratio + ", bandwidth=" + searchBand + ", local_radius=" + searchRadius + ", ";
            if(reduce) Command += " reduce";
            if(local) Command += " local_search";
            if(houghSeries) Command += " show_raw";
            if(showCircles) Command += " show_mask";
            if(showRadius) Command += " show_centroids";
            if(showScores) Command += " show_scores";
            if(results) Command += " results_table";
            Command += "\");\r\n";
            Recorder.recordString(Command);
        }
        // </editor-fold>

        //Create an ImagePlus instance of the currently active image
        imp = WindowManager.getCurrentImage();
        
        //Import the ImagePlus as a stack
        ImageStack stack = imp.getStack();
        
        //Get the ROI dimensions - no ROI = full image
        r = stack.getRoi();
        offx = r.x;
        offy = r.y;
        width = r.width;
        height = r.height;
        fullWidth = stack.getWidth();
        fullHeight = stack.getHeight();      
IJ.log("There are " + WindowManager.getNonImageWindows().length + " windows and " + WindowManager.getImageCount() + " images!");

        //Convert the stack to float (allows 8, 16, and 32 stacks to all be processed as one type)
        ImageStack stackCopy = stack.duplicate();
        ImageStack floatStack = stackCopy.convertToFloat();

        if(houghSeries){
            //Frames is transform dimension
            houghStack = new ImageStack(width, height, depth*imp.getNSlices());
        }
        if(showCircles){
            circleStack = new ImageStack(width, height, imp.getNSlices());
        }
        if(showRadius){
            radiusStack = new ImageStack(width, height, imp.getNSlices());
        }
        if(showScores){
            scoreStack = new ImageStack(width, height, imp.getNSlices());
        }

        //See if input is stack
        stackSlices = imp.getStackSize();
        if(stackSlices > 1) isStack = true;
        
        //Build arrays for storing the circle result parameters
        houghScores = new int [maxCircles];
        centerRadii = new int [maxCircles];
        centerPoint = new Point[maxCircles];
        circleID = new int [maxCircles];
        Arrays.fill(circleID, -1); //Flag indeces as -1 = no circle found for this index
        Arrays.fill(centerRadii, -1);
        Arrays.fill(centerPoint, new Point(-1,-1));
        Arrays.fill(houghScores, -1);
        
        //Start timer to delay log outputs
        long startTime = System.currentTimeMillis();
        
        //Retrieve the current slice in the stack as an ImageProcessor
        for(int slice = 1; slice <= stackSlices; slice++){
            //Store and then reset nCircles
            nCirlcesPrev = nCircles;
            nCircles = 0;
            
            //Set houghMaximum to -1 so that it can be determined whether the maximum has already been found
            maxHough = -1;
            
            //Initialize local search to false to record whether a local search was performed
            localSearch = false;
         
            //Show tranform status if there is more than one slice
            if(isStack){
               if((System.currentTimeMillis() - startTime) > 500){
                   IJ.showStatus("Hough tranform - processing frame: " + slice + " of " + stackSlices + "");
                   
                   //Reset timer
                   startTime = System.currentTimeMillis();
               }
            }
            else{
               IJ.log("\\Clear");
               IJ.log("Running Hough tansfrom...\r\n");
            }
 
            ImageProcessor ip = floatStack.getProcessor(slice);
            imageValues = (float[]) ip.getPixels();

            // <editor-fold desc="Local search">
            //If the number of circles is greater than min circles - speed up the search by looking locally for each next circle
            //Multithread by giving each core a subset of the total number of circles
            if(local & slice > 1){
                //If there are still a sufficient number of circles, perform an exclusively local search
                if(nCirlcesPrev>=minCircles){
                    method = "Local";
                    
                    //Set local search to true
                    localSearch = true;
                    startLocalTransform();
                   
                    //Rebuild the local searches into a full transform if the Hough Series is desired
                    if(houghSeries){
                        convertLocaltoFullHough(slice, houghStack);

                    }
                }
                
                //Otherwise, if there still are valid circles, perform the full Hough transform, but perform a local search for the valid circles first before continuing onto a full search 
                else if(nCirlcesPrev > 0){
                    method = "Partial Local";
                    houghTransform();
                    if (houghSeries){
                        //Create the hyperstach to put into the result if needed
                        HoughSpaceSeries(slice, houghStack);
                    }
                    startPartialLocalSearch();
                    getCenterPoints();
                }
                 //If there are no circles, and minCircles is 0, then return blank results
                else if(nCirlcesPrev == 0 & minCircles == 0){
                    method = "N/A";
                    //Create an empty Hough array
                    houghValues = new int[width][height][depth];
                    if (houghSeries){
                        //Create the hyperstach to put into the result if needed
                        HoughSpaceSeries(slice, houghStack);
                    }
                    
                    //Set maxHough to an arbitrary value to flag that there is no need to find this value
                    maxHough = 1;
                } 

                //If no circles are found, then revert to a full Hough
                else{
                    method = "Full";
                    houghTransform();
                    if (houghSeries){
                        //Create the hyperstach to put into the result if needed
                        HoughSpaceSeries(slice, houghStack);
                    }
                    // Mark the center of the found circles in a new image if user wants to find centers
                    if(showCircles || showRadius || showScores || results) getCenterPoints();
                }
            }
            //Otherwise, perform the full transform
            else{
                method = "Full";
                houghTransform();
                if (houghSeries){
                    //Create the hyperstach to put into the result if needed
                    HoughSpaceSeries(slice, houghStack);
                }
                // Mark the center of the found circles in a new image if user wants to find centers
                if(showCircles || showRadius || showScores || results) getCenterPoints();
            }
             // </editor-fold>
            
            if(!isStack) IJ.log("Searching Hough space for circles...\r\n");

            // Create image View for Marked Circles.
            if(showCircles) drawCircles(slice, circleStack, width, height, offx, offy, fullWidth);
            
            //Create map of centroids where the intensity is the radius
            if (showRadius || showScores) drawCentroids(slice, radiusStack, scoreStack);
            
            //Export measurements to the results table
            if (results) resultsTable(slice);
            
        }
IJ.log("" + totalTime);
        //Draw the resulting stacks
         if(houghSeries){
            houghPlus = new ImagePlus("Hough Transform Series", houghStack);
            if(depth>1) houghPlus = HyperStackConverter.toHyperStack(houghPlus, 1, depth, imp.getNSlices(), "default", "grayscale");
            houghPlus.show();
         }
         if(showCircles) new ImagePlus("Centroid overlay", circleStack).show();
         if(showRadius) new ImagePlus("Centroid map", radiusStack).show();
         if(showScores) new ImagePlus("Score map", scoreStack).show(); 
         if(results) rt.show("Results");
         
         IJ.log("Hough transform complete!");
    }
    
    //OPTMIZED
    private void startLocalTransform(){
         //Initialize an array to store the local Hough transform from each thread
        localHoughValues = new int [nCirlcesPrev][2*searchRadius + 1][2*searchRadius + 1][2*searchBand + 1];
        
        //Initialize an array to store the final dimensions and location of each local Hough Transform
        localHoughParameters = new int [nCirlcesPrev][9];

        //Setup two processing pipelines to multithread only if there are more than 1 circles
        if(nCirlcesPrev < 2){
            localHoughTransform(0);
            localGetCenterPoint(0); //Make sure this returns -1 to index if no circle was found
            maxHough = houghScores[0];
        }
        
        //Otherwise, multithread the local search
        else{
            //Build an array to store the result from each thread
            final Thread[] threads = newThreadArray();

            //Create an atomic integer counter that each thread can use to count through the radii
            final AtomicInteger ai = new AtomicInteger(0);  

            //Build a thread for as many CPUs as are available to the JVM 
            for (ithread = 0; ithread < threads.length; ithread++) {    

                // Concurrently run in as many threads as CPUs  
                threads[ithread] = new Thread() {  

                    { setPriority(Thread.NORM_PRIORITY); }  

                    @Override
                    public void run() {  

                        //Divide the task so that each core works on a subset of circles
                        for(int circleNum = ai.getAndAdd(1); circleNum < nCirlcesPrev; circleNum = ai.getAndAdd(1)){
                            localHoughTransform(circleNum);
                            localGetCenterPoint(circleNum); //Make sure this returns -1 to index if no circle was found
                        }
                    }               
                };  
            }    
            startAndJoin(threads);
            
            //Retrieve the maximum hough value if the raw Hough is desired
            if(houghSeries){
                for(int i=0; i<nCirlcesPrev; i++){
                    if(houghScores[i] > maxHough) maxHough = houghScores[i]; 
                }
            }
        }

        //Flag move all circles to the starting indexes, so there are no gaps between circles (i.e. if there were 3 circles, then they occupy indeces 0-2)
        collapseLocalResult();
    }
    
    //UNTESTED---------------------
    private void startPartialLocalSearch(){
        //Build an array to store the result from each thread
        final Thread[] threads = newThreadArray();

        //Create an atomic integer counter that each thread can use to count through the radii
        final AtomicInteger ai = new AtomicInteger(0);  

        //Build a thread for as many CPUs as are available to the JVM 
        for (ithread = 0; ithread < threads.length; ithread++) {    

            // Concurrently run in as many threads as CPUs  
            threads[ithread] = new Thread() {  

                { setPriority(Thread.NORM_PRIORITY); }  

                @Override
                public void run() {  

                    //Divide the task so that each core works on a subset of circles
                    for(int circleNum = ai.getAndAdd(1); circleNum < nCirlcesPrev; circleNum = ai.getAndAdd(1)){
                        getLocalCenterPoint2(circleNum);
                    }
                }               
            };  
        }    
        startAndJoin(threads);
        
        //Clear out the remainder of the previous circles in the circle information arrays
        for(int a = nCirlcesPrev; a<circleID.length; a++){
            circleID[a] = -1;
            centerPoint[a] = new Point(-1,-1);
            centerRadii[a] = -1;
            houghScores[a] = -1;
        }
        
        //Flag move all circles to the starting indexes, so there are no gaps between circles (i.e. if there were 3 circles, then they occupy indeces 0-2)
        collapseLocalResult();
    }
    
    //OPTMIZED - not time limiting
    //Build an array of the cartesion transformations necessary for each increment at each radius
    private int buildLookUpTable() {       
        //Build an array to store the X and Y coordinates for each angle increment (resolution) across each radius (depth)
        lut = new int[2][resolution][depth];
        
        //Initialize a variable that will allow for measuring the maximium LUT size
        int maxLUT = 0;
        
        //Step through all radii to be sampled
        for(int radius = radiusMax; radius>=radiusMin; radius -= radiusInc) {
            
            //Index counter that also tracks the largest actual LUT size (may be <= resolution)
            int i = 0; 
            for(int resStep = 0; resStep < resolution; resStep++) {
                
                //Calcualte the angle and corresponding X and Y displacements for the specified radius
                double angle = (2D*Math.PI * (double)resStep) / (double)resolution;
                int indexR = (radius-radiusMin)/radiusInc;
                int rcos = (int)Math.round ((double)radius * Math.cos (angle));
                int rsin = (int)Math.round ((double)radius * Math.sin (angle));
                
                //Test to make sure that the coordinate is a new coordinate
                //NOTE: A continuous circle is being discretized into pixels, it is possible that small angle steps for small radii circles will occupy the same pixel.
                //Since there is no point in making redundant calculations, these points are excluded from the LUT
                //NOTE: Using the minRadius as the transform cutoff results in a strong harmonic of the image forming near the min radius
                //threfore, using the max will push this harmonic outside the search range.
                if(radius == radiusMax && reduce){
                    if(i == 0) {
                        lut[0][i][indexR] = rcos;
                        lut[1][i][indexR] = rsin;
                        i++;
                    }
                    else if ((rcos != lut[0][i-1][indexR]) | (rsin != lut[1][i-1][indexR])){
                        lut[0][i][indexR] = rcos;
                        lut[1][i][indexR] = rsin;
                        i++;
                    }
                }
                else{
                    lut[0][i][indexR] = rcos;
                    lut[1][i][indexR] = rsin;
                    i++;
                }
            }
            
            //If this is the smallest radius, see how many transforms could be done, and set this as the new resolution
            if(radius == radiusMax){
                maxLUT = i;
                resolution = maxLUT;             
            }
        }
        return maxLUT;
    }

    //OPTIMIZED
    private void houghTransform () {
        //Build an array to store the result from each thread
        final Thread[] threads = newThreadArray();
        
        //Create an atomic integer counter that each thread can use to count through the radii
        final AtomicInteger ai = new AtomicInteger(radiusMin);  
        
        //Create an array to store the Hough values
        houghValues = new int[width][height][depth];
        
        //Build a thread for as many CPUs as are available to the JVM 
        for (ithread = 0; ithread < threads.length; ithread++) {    
            
            // Concurrently run in as many threads as CPUs  
            threads[ithread] = new Thread() {  
                          
                { setPriority(Thread.NORM_PRIORITY); }  
  
                @Override
                public void run() {  
  
                //Divide the radius tasks across the cores available
                for (int radius = ai.getAndAdd(radiusInc); radius <= radiusMax; radius = ai.getAndAdd(radiusInc)) {  
                    int indexR=(radius-radiusMin)/radiusInc;
                    //For a given radius, transform each pixel in a circle, and add-up the votes 
                    for(int y = 1; y < height-1; y++) {
                        for(int x = 1; x < width-1; x++) {
                                if( imageValues[(x+offx)+(y+offy)*fullWidth] != 0 )  {// Edge pixel found                                    
                                    for(int i = 0; i < lutSize; i++) {
                                        int a = x + lut[1][i][indexR]; 
                                        int b = y + lut[0][i][indexR]; 
                                        if((b >= 0) & (b < height) & (a >= 0) & (a < width)) {
                                            houghValues[a][b][indexR] += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }  
            };  
        }    
        startAndJoin(threads);      
    }
    
    //The local Hough is an inversion of the inertial reference frame used in the full Hough
    //In the full Hough the image is seach for pixels with a value > 1, and if this is true
    //then the pixel is projected in a circle about that point.
    
    //OPTMIZED
    //To reduce the necessary transform space,
    private void localHoughTransform (int index) {
        //Initialize local search variables
        int startWidth = centerPoint[index].x - searchRadius;
        if (startWidth < 1) startWidth = 1; //Keep search area inside the image
        int endWidth = centerPoint[index].x + searchRadius;
        if (endWidth > width) endWidth = width; //Keep search area inside the image
        int lwidth = endWidth - startWidth + 1;
        
        int startHeight = centerPoint[index].y - searchRadius;
        if (startHeight < 1) startHeight = 1; //Keep search area inside the image
        int endHeight = centerPoint[index].y + searchRadius;
        if (endHeight > height) endHeight = height; //Keep search area inside the image
        int lheight = endHeight - startHeight + 1;
        
        int lradiusMin = centerRadii[index] - searchBand;
        if(lradiusMin < radiusMin) lradiusMin = radiusMin;
        int lradiusMax = centerRadii[index] + searchBand;
        if(lradiusMax > radiusMax) lradiusMax = radiusMax;
        int ldepth = ((lradiusMax-lradiusMin)/radiusInc)+1;
        
        //Store local Hough parameters into array
        localHoughParameters[index][0] = startWidth;
        localHoughParameters[index][1] = endWidth;
        localHoughParameters[index][2] = lwidth;
        localHoughParameters[index][3] = startHeight;
        localHoughParameters[index][4] = endHeight;
        localHoughParameters[index][5] = lheight;
        localHoughParameters[index][6] = lradiusMin;
        localHoughParameters[index][7] = lradiusMax;
        localHoughParameters[index][8] = ldepth;
        
        //Divide the radius tasks across the cores available
        for (int radius = lradiusMin; radius <= lradiusMax; radius += radiusInc) { 
            int indexR=(radius-lradiusMin)/radiusInc;
            //For a given radius, transform each pixel in a circle, and add-up the votes 
            for(int y = startHeight; y < endHeight-1; y++) {
                for(int x = startWidth; x < endWidth-1; x++) {
                    for(int i = 0; i < lutSize; i++){
                        int a = x + offx + lut[1][i][indexR+((lradiusMin-radiusMin)/radiusInc)]; 
                        int b = y + offy + lut[0][i][indexR+((lradiusMin-radiusMin)/radiusInc)];
                        
                        //Check to make sure pixel is within the image
                        if((a>1) & (a<fullWidth-1) & (b>1) & (b<fullHeight-1)){
                            //See if the pixel has an intensity > 1, if so, add 1 to the vote
                            if(imageValues[a + (b*fullWidth)] > 0 ){
                                localHoughValues[index][x-startWidth][y-startHeight][indexR] += 1;
                            } 
                        }
                    }
                }
            }
        }
    }  

    //OPTMIZED
    //Find the largest Hough pixel in the 3D Hough transform array to scale the 8-bit conversion
    private void houghMaximum () {
        //long startTime = System.currentTimeMillis(); //1337ms without multi, 319ms with multi, 175ms by writing to private variable per thread
        //Build an array to store the result from each thread
        final Thread[] threads = newThreadArray();
        
        //Build an array to store the max from each thread
        maxHoughArray = new int[threads.length][4];
        
        //Create an atomic integer counter that each thread can use to count through the radii
        final AtomicInteger ai = new AtomicInteger(0);
        
        //Create an integer for indexing the results array (one result per thread
        final AtomicInteger Index = new AtomicInteger(0);

        maxHough = 0;
        
        //Build a thread for as many CPUs as are available to the JVM 
        for (ithread = 0; ithread < threads.length; ithread++) {    
            
            // Concurrently run in as many threads as CPUs  
            threads[ithread] = new Thread() {  
                          
                { setPriority(Thread.NORM_PRIORITY); }  

                //Search for the largest score with each thread
                @Override
                public void run() {
                    int maxHoughThread = -1;
                    int maxRadiusThread = -1;
                    Point maxPointThread = new Point (-1,-1);
                    for(int a=ai.getAndIncrement(); a<depth; a=ai.getAndIncrement()){
                        for(int j = 0; j < height; j++) {
                            for(int k = 0; k < width; k++){
                                if(houghValues[k][j][a] > maxHoughThread) {
                                    maxHoughThread = houghValues[k][j][a];
                                    maxPointThread = new Point(k,j);
                                    maxRadiusThread = a*radiusInc + radiusMin;
                                }
                            }
                        }
                    }
                    //Have each thread report the score to a common array
                    maxHoughArray[Index.getAndIncrement()] = new int[]{maxHoughThread, maxRadiusThread, maxPointThread.x, maxPointThread.y};                    
                }
            };
        }
        startAndJoin(threads);
        
        //Search common array for highest score
        for (int[] maxHoughArray1 : maxHoughArray) {
            if (maxHoughArray1[0] > maxHough) {
                maxHough = maxHoughArray1[0];
                maxRadius = maxHoughArray1[1];
                maxPoint = new Point((int) maxHoughArray1[2], (int) maxHoughArray1[3]);    
            }
        }
        //long stopTime = System.currentTimeMillis();
        //IJ.log("Elapsed time was " + (stopTime - startTime) + " miliseconds.");
    }
    
    //OPTMIZED
    //Create a Hough stack series using the local transforms
    private void convertLocaltoFullHough(int slice, ImageStack houghStack){
       
        
        int startFrame = ((slice-1)*depth);

        //Create an array to store the Hough values
        byte localHoughPixels[][] = new byte[depth][width*height];
        
        //NOTE: Only do multithreading if the problem is sufficiently complex
        if(localHoughParameters.length < 2 & localHoughParameters.length*searchBand*searchRadius*searchRadius < 15000){
            //Extract the local search parameters
            //Divide the task so that each core works on a subset of circles
            for(int circleNum = 0; circleNum < localHoughParameters.length; circleNum++){
                //Extract the local search parameters
                int startWidth = localHoughParameters[circleNum][0];
                int endWidth = localHoughParameters[circleNum][1];
                int lwidth = localHoughParameters[circleNum][2];
                int startHeight = localHoughParameters[circleNum][3];
                int endHeight = localHoughParameters[circleNum][4];
                int lheight = localHoughParameters[circleNum][5];
                int lradiusMin = localHoughParameters[circleNum][6];
                int lradiusMax = localHoughParameters[circleNum][7];
                int ldepth = localHoughParameters[circleNum][8];
                for(int h=0; h<lheight; h++){    
                    for(int w=0; w<lwidth; w++){
                        for(int i=0; i<ldepth; i++){
                           localHoughPixels[i+(lradiusMin-radiusMin)/radiusInc][(w+startWidth)+(h+startHeight)*width] = (byte) Math.round ((localHoughValues[circleNum][w][h][i] * 255D) / maxHough);
                        }
                    }
                }
            }
        }
        
        //Otherwise, multithread the depositing of pixels
        else{
            //Build an array to store the result from each thread
            final Thread[] threads = newThreadArray();

            //Create an atomic integer counter that each thread can use to count through the radii
            final AtomicInteger ai = new AtomicInteger(0);  

            //Build a thread for as many CPUs as are available to the JVM 
            for (ithread = 0; ithread < threads.length; ithread++) {    

                // Concurrently run in as many threads as CPUs  
                threads[ithread] = new Thread() {  

                    { setPriority(Thread.NORM_PRIORITY); }  

                    @Override
                    public void run() {  

                        //Divide the task so that each core works on a subset of circles
                        for(int circleNum = ai.getAndAdd(1); circleNum < localHoughParameters.length; circleNum = ai.getAndAdd(1)){
                            //Extract the local search parameters
                            int startWidth = localHoughParameters[circleNum][0];
                            int endWidth = localHoughParameters[circleNum][1];
                            int lwidth = localHoughParameters[circleNum][2];
                            int startHeight = localHoughParameters[circleNum][3];
                            int endHeight = localHoughParameters[circleNum][4];
                            int lheight = localHoughParameters[circleNum][5];
                            int lradiusMin = localHoughParameters[circleNum][6];
                            int lradiusMax = localHoughParameters[circleNum][7];
                            int ldepth = localHoughParameters[circleNum][8];
                            for(int h=0; h<lheight; h++){    
                                for(int w=0; w<lwidth; w++){
                                    for(int i=0; i<ldepth; i++){
                                       localHoughPixels[i+(lradiusMin-radiusMin)/radiusInc][(w+startWidth)+(h+startHeight)*width] = (byte) Math.round ((localHoughValues[circleNum][w][h][i] * 255D) / maxHough);
                                    }
                                }
                            }
                        }
                    }               
                };  
            }    
            startAndJoin(threads);
        }
        
        //Deposit the array into the Hough Series stack
        //Not time limiting, even at fairly high resolutions
        for(int radius = radiusMin; radius<=radiusMax; radius += radiusInc) {
            //Calculate the corresponding index
            int houghIndex = (radius-radiusMin)/radiusInc;
            
            //Deposit the array image into the corresponding slice in the stack
            houghStack.setPixels(localHoughPixels[houghIndex], houghIndex+1+startFrame);

            //Give the current slice the appropriate radius label
            houghStack.setSliceLabel("Hough Space [r="+radius+", resolution="+resolution+"]", houghIndex+1+startFrame);
        }
       
    }
    
    //OPTIMIZED
    //Add transform series data to the hyperstack
    private void HoughSpaceSeries(int slice, ImageStack houghStack){  
        //If the maximum Hough value has not yet been assigned, search the whole Hough transform for the maximum value
        if(maxHough == -1){
            houghMaximum();  
        }

        //Create an array to store the Hough values
        byte localHoughPixels[][] = new byte[depth][width*height];        
        int startFrame = ((slice-1)*depth);

        //If a full transform was done, calcualte the Hough Pixels
        if(minCircles > 0 | nCirlcesPrev > 0){
            //Build an array to store the result from each thread
            final Thread[] threads = newThreadArray();

            //Create an atomic integer counter that each thread can use to count through the radii
            final AtomicInteger ai = new AtomicInteger(radiusMin);  

            //Build a thread for as many CPUs as are available to the JVM 
            for (ithread = 0; ithread < threads.length; ithread++) {    

                // Concurrently run in as many threads as CPUs  
                threads[ithread] = new Thread() {  

                    { setPriority(Thread.NORM_PRIORITY); }  

                    @Override
                    public void run() {  

                    //Divide the radius tasks across the cores available
                    for (int radius = ai.getAndAdd(radiusInc); radius <= radiusMax; radius = ai.getAndAdd(radiusInc)) {

                        //Calculate the corresponding index
                        int houghIndex = (radius-radiusMin)/radiusInc;

                        //If full tansform was performed, retrieve the pixel array for the current Hough radius image
                        createHoughPixels(localHoughPixels[houghIndex], houghIndex);

                        //Deposit the array image into the corresponding slice in the stack
                        houghStack.setPixels(localHoughPixels[houghIndex], houghIndex+1+startFrame);

                        //Give the current slice the appropriate radius label
                        houghStack.setSliceLabel("Hough Space [r="+radius+", resolution="+resolution+"]", houghIndex+1+startFrame);
                    }
                }  
            };
            }
            startAndJoin(threads);
        }
    
        //Deposit the array into the Hough Series stack
        //Not time limiting, even at fairly high resolutions
        for(int radius = radiusMin; radius<=radiusMax; radius += radiusInc) {
            //Calculate the corresponding index
            int houghIndex = (radius-radiusMin)/radiusInc;
            
            //Deposit the array image into the corresponding slice in the stack
            houghStack.setPixels(localHoughPixels[houghIndex], houghIndex+1+startFrame);

            //Give the current slice the appropriate radius label
            houghStack.setSliceLabel("Hough Space [r="+radius+", resolution="+resolution+"]", houghIndex+1+startFrame);
        }
    }
    
    //OPTMIZED
    // Convert Values in Hough Space to an 8-Bit Image Space.
    private void createHoughPixels (byte houghPixels[], int index) {
	//Rescale all the Hough values to 8-bit to create the Hough image - 47ms to complete - single threading okay
        for(int l = 0; l < height; l++) {
            for(int i = 0; i < width; i++) {
                houghPixels[i + l * width] = (byte) Math.round ((houghValues[i][l][index] * 255D) / maxHough);
            }

        }
    }

    // Draw the circles found in the original image.
    private void drawCircles(int slice, ImageStack circleStack, int widthROI, int heightROI, int fullWidthX, int fullWidthY, int fullWidthROI) {
		
            // Copy original input pixels into output
            // circle location display image and
            // combine with saturation at 100
            byte[] circlespixels = new byte [widthROI*heightROI];

            int roiaddr=0;
            for( int y = fullWidthY; y < fullWidthY+heightROI; y++) {
                    for(int x = fullWidthX; x < fullWidthX+widthROI; x++) {
                            // Copy;
                            circlespixels[roiaddr] = (byte) imageValues[x+fullWidthROI*y];
                            // Saturate
                            if(circlespixels[roiaddr] != 0 )
                                    circlespixels[roiaddr] = 100;
                            else
                                    circlespixels[roiaddr] = 0;
                            roiaddr++;
                    }
            }
            // Copy original image to the circlespixels image.
            // Changing pixels values to 100, so that the marked
            // circles appears more clear. Must be improved in
            // the future to show the resuls in a colored image.
            //for(int i = 0; i < width*height ;++i ) {
            //if(imageValues[i] != 0 )
            //if(circlespixels[i] != 0 )
            //circlespixels[i] = 100;
            //else
            //circlespixels[i] = 0;
            //}
            if(centerPoint == null) getCenterPoints();

            byte cor = -1;
            // Redefine these so refer to ROI coordinates exclusively
            fullWidthROI= widthROI;
            fullWidthX=0;
            fullWidthY=0;

            for(int l = 0; l < nCircles; l++) {
                    int i = centerPoint[l].x;
                    int j = centerPoint[l].y;
                    // Draw a gray cross marking the center of each circle.
                    for( int k = -10 ; k <= 10 ; ++k ) {
                            int p = (j+k+fullWidthY)*fullWidthROI+ (i+fullWidthX);
                            if(!outOfBounds(j+k+fullWidthY,i+fullWidthX))
                                    circlespixels[(j+k+fullWidthY)*fullWidthROI+ (i+fullWidthX)] = cor;
                            if(!outOfBounds(j+fullWidthY,i+k+fullWidthX))
                                    circlespixels[(j+fullWidthY)*fullWidthROI  + (i+k+fullWidthX)] = cor;
                    }
                    for( int k = -2 ; k <= 2 ; ++k ) {
                            if(!outOfBounds(j-2+fullWidthY,i+k+fullWidthX))
                                    circlespixels[(j-2+fullWidthY)*fullWidthROI+ (i+k+fullWidthX)] = cor;
                            if(!outOfBounds(j+2+fullWidthY,i+k+fullWidthX))
                                    circlespixels[(j+2+fullWidthY)*fullWidthROI+ (i+k+fullWidthX)] = cor;
                            if(!outOfBounds(j+k+fullWidthY,i-2+fullWidthX))
                                    circlespixels[(j+k+fullWidthY)*fullWidthROI+ (i-2+fullWidthX)] = cor;
                            if(!outOfBounds(j+k+fullWidthY,i+2+fullWidthX))
                                    circlespixels[(j+k+fullWidthY)*fullWidthROI+ (i+2+fullWidthX)] = cor;
                    }
            }
            circleStack.setPixels(circlespixels, slice);
            circleStack.setSliceLabel(nCircles + " circles found", slice);
    }
    
    // Draw the centroids found in the original image where intensity = radius.
    private void drawCentroids(int slice, ImageStack radiusStack, ImageStack scoreStack) {
            
        //Create arrays the same size as the images
        radiusip = new ShortProcessor(width, height);
        short[] radiuspixels = (short[])radiusip.getPixels();

        scoresip = new ShortProcessor(width, height);
        short[] scorepixels = (short[])scoresip.getPixels();
        
        for(int l = 0; l < nCircles; l++) {
                int i = centerPoint[l].x;
                int j = centerPoint[l].y;
                int radius = centerRadii[l];

                //Draw a point at the centroid and set the int to = radius
                if(showRadius) radiuspixels[j*fullWidth + i] = (short) radius;

                //Record the Hough score on a 16-bit scale
                if(showScores) scorepixels[j*fullWidth + i] = (short) Math.round(houghScores[l]/resolution*65535);
        }
        
        if(showRadius){
            radiusStack.setPixels(radiuspixels, slice);
            radiusStack.setSliceLabel(nCircles + " circles found", slice);
        }
        if(showScores){
            scoreStack.setPixels(scorepixels, slice);
            scoreStack.setSliceLabel(nCircles + " circles found", slice);
        }
    }
    
    //Export the results to the results table
    public void resultsTable(int frame){
        for(int a = 0; a < nCircles; a++){
            rt.incrementCounter();
            rt.addValue("ID", circleID[a]);
            rt.addValue("X", centerPoint[a].x);
            rt.addValue("Y", centerPoint[a].y);
            rt.addValue("Radius", centerRadii[a]);
            rt.addValue("Score", ((float) houghScores[a]/resolution));
            rt.addValue("nCircles", nCircles);
            rt.addValue("Resolution", lutSize);
            rt.addValue("Frame", frame);
            rt.addValue("Method", method);
        }
    }
    

    private boolean outOfBounds(int y,int x) {
        if(x >= width)
            return(true);
        if(x <= 0)
            return(true);
        if(y >= height)
            return(true);
        if(y <= 0)
            return(true);
        return(false);
    }

    //OPTMIZED
    /** Search for a fixed number of circles.
    @param maxCircles The number of circles that should be found.  
    */
    private void getCenterPoints () {
   
        int countCircles = nCircles;
        maxHough = threshold;

        while(countCircles < maxCircles && maxHough >= threshold) {

            //Search for the highest remaining Hough score in the matrix
            houghMaximum();

            if(maxHough>=threshold){
                circleIDcounter++;
                circleID[countCircles] = circleIDcounter;
                centerPoint[countCircles] = maxPoint;
                centerRadii[countCircles] = maxRadius;
                houghScores[countCircles] = maxHough;
                countCircles++;
                clearNeighbours(maxPoint.x,maxPoint.y, maxRadius);
            }
        }

        nCircles = countCircles;
        
        //Clear the remainder of result arrays, since these results are from the previous round
        for(int i = nCircles; i<maxCircles; i++){
            circleID[i] = -1;
            centerPoint[i] = new Point(-1,-1);
            centerRadii[i] = -1;
            houghScores[i] = -1;
        }
    }
    
    //UNTESTED------------------------
    private void getLocalCenterPoint2(int index){
        //Initialize local search variables
        //Initialize local search variables
        int startWidth = centerPoint[index].x - searchRadius;
        if (startWidth < 1) startWidth = 1; //Keep search area inside the image
        int endWidth = centerPoint[index].x + searchRadius;
        if (endWidth > width) endWidth = width; //Keep search area inside the image
        
        int startHeight = centerPoint[index].y - searchRadius;
        if (startHeight < 1) startHeight = 1; //Keep search area inside the image
        int endHeight = centerPoint[index].y + searchRadius;
        if (endHeight > height) endHeight = height; //Keep search area inside the image
        
        int lradiusMin = centerRadii[index] - searchBand;
        if(lradiusMin < radiusMin) lradiusMin = radiusMin;
        int lradiusMax = centerRadii[index] + searchBand;
        if(lradiusMax > radiusMax) lradiusMax = radiusMax;

        
        //Search locally for the highest hough score in the full Hough Space
        int maxLocalHough = -1;
        int maxLocalRadius = -1;
        Point maxLocalPoint = new Point (-1,-1);
        for(int radius = lradiusMin; radius<=lradiusMax; radius += radiusInc){
            int indexR=(radius-radiusMin)/radiusInc;
            for(int j = startHeight; j < endHeight; j++) {
                for(int k = startWidth; k < endWidth; k++){
                    if(houghValues[k][j][indexR] > maxLocalHough) {
                        maxLocalHough = houghValues[k][j][indexR];
                        maxLocalPoint = new Point(k,j);
                        maxLocalRadius = radius;
                    }
                }
            }
        }       
        //If the highest score is above the threshold, record the new circle
        if(maxLocalHough>=threshold){
            centerPoint[index] = maxLocalPoint;
            centerRadii[index] = maxLocalRadius;
            houghScores[index] = maxLocalHough;
            clearNeighbours(maxLocalPoint.x,maxLocalPoint.y, maxLocalRadius);
        }
        //Otherwise, record that the circle was lost
        else{
            circleID[index] = -1;
            centerPoint[index] = new Point(-1,-1);
            centerRadii[index] = -1;
            houghScores[index] = -1;
        }  
    }
    
    //OPTMIZED
    private void localGetCenterPoint(int index){
        
        //Extract the local search parameters
        int startWidth = localHoughParameters[index][0];
        int endWidth = localHoughParameters[index][1];
        int lwidth = localHoughParameters[index][2];
        int startHeight = localHoughParameters[index][3];
        int endHeight = localHoughParameters[index][4];
        int lheight = localHoughParameters[index][5];
        int lradiusMin = localHoughParameters[index][6];
        int lradiusMax = localHoughParameters[index][7];
        int ldepth = localHoughParameters[index][8];
        
        //Search for the highest hough score in the local Hough Space
        int maxLocalHough = -1;
        int maxLocalRadius = -1;
        Point maxLocalPoint = new Point (-1,-1);
        for(int a=0; a<ldepth; a++){
            for(int j = 0; j < lheight; j++) {
                for(int k = 0; k < lwidth; k++){
                    if(localHoughValues[index][k][j][a] > maxLocalHough) {
                        maxLocalHough = localHoughValues[index][k][j][a];
                        maxLocalPoint = new Point(k + startWidth,j + startHeight);
                        maxLocalRadius = a*radiusInc + lradiusMin;
                    }
                }
            }
        }
        
        //If the highest score is above the threshold, record the new circle
        if(maxLocalHough>=threshold){            
            centerPoint[index] = maxLocalPoint;
            centerRadii[index] = maxLocalRadius;
            houghScores[index] = maxLocalHough;
        }
        //Otherwise, record that the circle was lost
        else{
            circleID[index] = -1;
            centerPoint[index] = new Point(-1,-1);
            centerRadii[index] = -1;
            houghScores[index] = -1;
        }     
    }

    //OPTMIZED
    private void collapseLocalResult(){
        //search for indeces containing circles (i.e. index != -1)
        int[] idIndeces = new int [circleID.length];
        int indexCounter = 0;
        for(int i=0; i<circleID.length; i++){
            //If a valid ID is found, record the index
            if(circleID[i] != -1){
                idIndeces[indexCounter] = i;
                indexCounter++;
            }
        }
       
        //Move all the found results to the starting indeces of the arrays if needed
        if(indexCounter < maxCircles){
            for(int i=0; i<indexCounter; i++){
                //If index is empty, then move the result
                if(circleID[i] == -1){
                    //Check to see index is empty
                    //Move results
                    circleID[i] = circleID[idIndeces[i]];
                    houghScores[i] = houghScores[idIndeces[i]];
                    centerRadii[i] = centerRadii[idIndeces[i]];
                    centerPoint[i] = centerPoint[idIndeces[i]];

                    //Clear original index
                    circleID[idIndeces[i]] = -1;
                    houghScores[idIndeces[i]] = -1;
                    centerRadii[idIndeces[i]] = -1;
                    centerPoint[idIndeces[i]] = new Point(-1,-1);
                }
            }
        }
        
        //Update nCircles to reflect the new found number of circles
        nCircles = indexCounter;
    };
    
    //OPTIMIZED - Not time limiting, even with large circles
    /** Clear, from the Hough Space, all the counter that are near (radius/2) a previously found circle C.    
    @param x The x coordinate of the circle C found.
    @param x The y coordinate of the circle C found.
    @param x The radius of the circle C found.
    */
    private void clearNeighbours(int x,int y, int radius) {
        // The following code just clean the points around the center of the circle found.
	int radiusSquared = (int) Math.round(Math.pow(radius*ratio, 2D));
        //int radiusSquared = radius*radius;

        int y1 = (int)Math.floor (y - radius);
        int y2 = (int)Math.ceil (y + radius) + 1;
        int x1 = (int)Math.floor (x - radius);
        int x2 = (int)Math.ceil (x + radius) + 1;

        if(y1 < 0)
            y1 = 0;
        if(y2 > height)
            y2 = height;
        if(x1 < 0)
            x1 = 0;
        if(x2 > width)
            x2 = width;

        for(int indexR = 0; indexR<depth; indexR++){
            for(int i = y1; i < y2; i++) {
                for(int j = x1; j < x2; j++) {	      	     
                    if((int) (Math.pow (j - x, 2D) + Math.pow (i - y, 2D)) < radiusSquared) {
                        houghValues[j][i][indexR] = 0;
                    }
                }
            }
        }
    }
    
    /** Create a Thread[] array as large as the number of processors available. 
    * From Stephan Preibisch's Multithreading.java class. See: 
    * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD 
    */  
    private Thread[] newThreadArray() {  
        int n_cpus = Runtime.getRuntime().availableProcessors();  
        return new Thread[n_cpus];  
    }
    
    /** Start all given threads and wait on each of them until all are done. 
    * From Stephan Preibisch's Multithreading.java class. See: 
    * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD 
     * @param threads
    */  
    public static void startAndJoin(Thread[] threads)  
    {  
        for (int ithread = 0; ithread < threads.length; ++ithread)  
        {  
            threads[ithread].setPriority(Thread.NORM_PRIORITY);  
            threads[ithread].start();  
        }  
  
        try  
        {     
            for (int ithread = 0; ithread < threads.length; ++ithread)  
                threads[ithread].join();  
        } catch (InterruptedException ie)  
        {  
            throw new RuntimeException(ie);  
        }  
    } 
}