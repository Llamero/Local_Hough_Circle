/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Hough_Package;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import static ij.plugin.filter.PlugInFilter.DOES_16;
import static ij.plugin.filter.PlugInFilter.DOES_32;
import static ij.plugin.filter.PlugInFilter.DOES_8G;
import static ij.plugin.filter.PlugInFilter.DONE;
import static ij.plugin.filter.PlugInFilter.SUPPORTS_MASKING;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.KeyEventDispatcher;
import java.awt.KeyboardFocusManager;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeEvent;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SwingWorker.StateValue;
import javax.swing.UIManager;

/**
 *
 * @author Ben
 */
public class Hough_GUI implements PlugInFilter {
    //***GUI input variables***
    // <editor-fold desc="Initialize variables">  
    final private JPanel guiPanel = new JPanel(); //Initialize a panel (needed to house a frame)
    final private JFrame guiFrame = new JFrame();    //Initlialize a frame to hold the gui
    final private JLabel guiTitle = new javax.swing.JLabel();
    final private JLabel guiIntro1 = new javax.swing.JLabel();
    final private JLabel guiIntro2 = new javax.swing.JLabel();
    final private JRadioButton guiEasyModeButton = new javax.swing.JRadioButton();
    final private JRadioButton guiAdvancedModeButton = new javax.swing.JRadioButton();
    final private ButtonGroup modeButtonGroup = new javax.swing.ButtonGroup();
    final private JLabel guiSearchLabel = new javax.swing.JLabel();
    final private JLabel guiMinLabel = new javax.swing.JLabel();
    final private JTextField guiMinText = new javax.swing.JTextField();
    final private JLabel guiMaxLabel = new javax.swing.JLabel();
    final private JTextField guiMaxText = new javax.swing.JTextField();
    final private JLabel guiIncLabel = new javax.swing.JLabel();
    final private JTextField guiIncText = new javax.swing.JTextField();
    final private JLabel  guiMinNumLabel = new javax.swing.JLabel();
    final private JTextField guiMinNumText = new javax.swing.JTextField();
    final private JLabel guiMaxNumLabel = new javax.swing.JLabel();
    final private JTextField guiMaxNumText = new javax.swing.JTextField();       
    final private JLabel guiThreshLabel = new javax.swing.JLabel();
    final private JTextField guiThreshText = new javax.swing.JTextField();
    final private JLabel guiResLabel = new javax.swing.JLabel();
    final private JTextField guiResText = new javax.swing.JTextField();
    final private JLabel guiClearLabel = new javax.swing.JLabel();
    final private JTextField guiClearText = new javax.swing.JTextField();
    final private JLabel guiRadiusBandLabel = new javax.swing.JLabel();
    final private JTextField guiRadiusBandText = new javax.swing.JTextField();
    final private JLabel guiSearchRadLabel = new javax.swing.JLabel();
    final private JTextField guiSearchRadText = new javax.swing.JTextField();
    final private JCheckBox guiReduceBox = new javax.swing.JCheckBox();
    final private JCheckBox guiLocalBox = new javax.swing.JCheckBox();
    final private JLabel guiOutputLabel = new javax.swing.JLabel();
    final private JCheckBox guiRawBox = new javax.swing.JCheckBox();
    final private JCheckBox guiPointBox = new javax.swing.JCheckBox();
    final private JCheckBox guiIDBox = new javax.swing.JCheckBox();
    final private JCheckBox guiHoughBox = new javax.swing.JCheckBox();
    final private JCheckBox guiResultsBox = new javax.swing.JCheckBox();
    private JProgressBar guiProgressBar; //Do not initialize so that color of font can be changed
    final private JButton guiOKButton = new javax.swing.JButton();
    private String barString = ""; 
    
    //Search parameters
    private int radiusMin; // Find circles with radius grater or equal radiusMin - argument syntax: "min=#"
    private int radiusMax; // Find circles with radius less or equal radiusMax - argument syntax: "max=#"
    private int radiusInc; // Increment used to go from radiusMin to radiusMax - argument syntax: "inc=#"
    private int minCircles;// Minumum number of circles to be found - argument syntax: "minCircles=#"    
    private int maxCircles;// Maximum number of circles to be found - argument syntax: "maxCircles=#"
    private double thresholdRatio;//Ratio input from GUI that expresses threshold as ratio of resolution (highest possible # of votes)
    private int resolution;//The number of steps to use per transform (i.e. number of voting rounds)
    private double ratio;// Ratio of found circle radius to clear out surrounding neighbors
    private int searchBand = 0;//The +/- range of radii to search for relative to the last found radius - argument syntax: "bandwidth=#"
    private int searchRadius = 0;//The search radius to look for the next centroid relative to the last found centroid - argument syntax: "radius=#"
    private boolean reduce = false;//Cap the transform resolution by removeing redundant steps
    private boolean local = false;//Whether or not the search is going to be local
    
    //Output parameters
    private boolean houghSeries = false;//Contains whether the user wants the Hough series stack as an output - argument syntax: "show_raw"
    private boolean showCircles = false;//Contains whether the user wants the circles found as an output - argument syntax: "show_mask"
    private boolean showID = false;//Contains whether the user wants a map of centroids and radii outputed from search - argument syntax: "show_centroids"
    private boolean showScores = false;//Contains whether the user wants a map of centroids and Hough scores outputed from search - argument syntax: "show_scores"
    private boolean results = false;//Contains whether the user wants to export the measurements to a reuslts table 
    
    //Keep track of whether analysis has started, and update GUI accordingly
    private boolean analysisStarted = false;
    private boolean isGUI;
    private List<String> status;
    
    //Start instance of analysis class
    Hough_Circle guiInput;
    
    @Override
    public int setup(String arg, ImagePlus imp) {        
        
        if (arg.equals("about")) {
            showAbout();
            return DONE;
        }
        
        //Sends arduments to ImageJ that tells it how to run the plugin - tells it to accept all grayscale and supports selections
        return DOES_8G+DOES_16+DOES_32+SUPPORTS_MASKING;//Do not include DOES_STACKS, as this will call the GUI once for each slice
    }
    
    void showAbout() {
        IJ.showMessage("Hough Circle Transform v1.0.0",
                       "This plugin performs a Hough circle transform \n" +
                       "on an image or stack.  Hough circle transforms\n" +
                       "can be used to find the centroid and radius of\n" +
                       "circles embedded in complex images.  This plugin\n"+
                       "was inspired by the transform implementation\n"+
                       "of Hemerson Pistori (pistori@ec.ucdb.br)"
                      );
    }
    
    @Override
    public void run(ImageProcessor ip) {
        //Show a Dialog Window for user input of parameters
        readParameters();
    } 
    
    void readParameters() {
        // <editor-fold desc="Retrieve macro arguments">
        //See if any arguments have been passed to the plugin via a macro
        if (IJ.isMacro() && ij.Macro.getOptions() != null && !ij.Macro.getOptions().trim().isEmpty()){
            
            isGUI = false;
            String[] arguments = ij.Macro.getOptions().split(",");

            //remove all spaces from array components
            for(int a = 0; a<arguments.length; a++){
                arguments[a] = arguments[a].trim();
            }
            
            //Parse the arguments to the corresponding variables
            for (String argument : arguments) { //passes each argment in the array to the variable "argument"
                if (argument.matches(".*minRadius.*=.*")) {
                    //Retrieve min radius
                    radiusMin = Integer.parseInt(argument.replaceAll("\\D+", ""));//Remove all non digits
                }
                else if (argument.matches(".*maxRadius.*=.*")) {
                    //Retrieve max radius
                    radiusMax = Integer.parseInt(argument.replaceAll("\\D+", ""));//Remove all non digits
                }
                else if (argument.matches(".*inc.*=.*")) {
                    //Retrieve radius increment
                    radiusInc = Integer.parseInt(argument.replaceAll("\\D+", ""));//Remove all non digits
                }
                else if (argument.matches(".*minCircles.*=.*")) {
                    //Retrieve number of circles
                    minCircles = Integer.parseInt(argument.replaceAll("\\D+", ""));//Remove all non digits
                }
                else if (argument.matches(".*maxCircles.*=.*")) {
                    //Retrieve number of circles
                    maxCircles = Integer.parseInt(argument.replaceAll("\\D+", ""));//Remove all non digits
                }
                else if (argument.matches(".*threshold.*=.*")) {
                    //Retrieve Hough score threshold
                    //Code from: http://www.itgo.me/a/x8683194173055835262/get-float-or-integer-value-from-the-string-in-java
                    Pattern pattern = Pattern.compile("\\d+(?:\\.\\d+)?");// Match int or float
                    Matcher matcher = pattern.matcher(argument);
                    if(matcher.find()) thresholdRatio = Double.parseDouble(matcher.group());
                }
                else if (argument.matches(".*resolution.*=.*")) {
                    //Retrieve Hough score threshold
                    resolution = Integer.parseInt(argument.replaceAll("\\D+", ""));//Remove all non digits
                }
                else if (argument.matches(".*ratio.*=.*")) {
                    //Retrieve clearing neighbor radius ratio
                    //Code from: http://www.itgo.me/a/x8683194173055835262/get-float-or-integer-value-from-the-string-in-java
                    Pattern pattern = Pattern.compile("\\d+(?:\\.\\d+)?");// Match int or float
                    Matcher matcher = pattern.matcher(argument);
                    if(matcher.find()) ratio = Double.parseDouble(matcher.group());
                }
                else if (argument.matches(".*bandwidth.*=.*")) {
                    //Retrieve Hough score threshold
                    searchBand = Integer.parseInt(argument.replaceAll("\\D+", ""));//Remove all non digits
                }
                else if (argument.matches(".*local_radius.*=.*")) {
                    //Retrieve Hough score threshold
                    searchRadius = Integer.parseInt(argument.replaceAll("\\D+", ""));//Remove all non digits
                }

                //Retrieve checkbox status
                if (argument.matches(".*reduce.*")) reduce = true;
                if (argument.matches(".*local_search.*")) local = true;
                if (argument.matches(".*show_raw.*")) houghSeries = true;
                if (argument.matches(".*show_mask.*")) showCircles = true;
                if (argument.matches(".*show_centroids.*")) showID = true;
                if (argument.matches(".*show_scores.*")) showScores = true;
                if (argument.matches(".*results_table.*")) results = true;
            }
            
            //Start the Hough Transform
            startTransform();
          
        }
        // </editor-fold>
        //If arguments are not coming from macro, present user interface
        else{
            // <editor-fold desc="Swing GUI part 1.">
            //***Build GUI using Swing***
            //Set the frame to close when the window is closed
            isGUI = true;
            guiFrame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
            
            //Format GUI text
            guiTitle.setFont(new java.awt.Font("Tahoma", 0, 24));// NOI18N
            guiTitle.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
            guiTitle.setText("Hough Circle Transform");

            guiIntro1.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
            guiIntro1.setText("This plugin performs a Hough circle transform on an image or stack. ");

            guiIntro2.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
            guiIntro2.setText("It can be used to find and measure circular objects within an image.");

            modeButtonGroup.add(guiEasyModeButton);
            guiEasyModeButton.setSelected(true);
            guiEasyModeButton.setText("Easy Mode");
            guiEasyModeButton.addActionListener((java.awt.event.ActionEvent evt) -> {
                //Setup local easy GUI
                if (guiLocalBox.isSelected()){
                    easyLocalGUI();
                }
                //Setup full easy GUI
                else{
                    easyFullGUI();
                }
            });

            modeButtonGroup.add(guiAdvancedModeButton);
            guiAdvancedModeButton.setText("Advanced Mode");
            guiAdvancedModeButton.addActionListener((java.awt.event.ActionEvent evt) -> {
                if (guiLocalBox.isSelected()){
                    advancedLocalGUI();
                }
                //Setup full advanced GUI
                else{
                    advancedFullGUI();
                }
            });

            guiSearchLabel.setFont(new java.awt.Font("Tahoma", 0, 14));// NOI18N
            guiSearchLabel.setText("Search parameters:");

            guiMinLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiMinLabel.setText("Minimum search radius (in pixels):");

            guiMinText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiMinText.setText("10");
            guiMinText.setToolTipText("Radius of smallest circle you expect to find");

            guiMaxLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiMaxLabel.setText("Maximum search radius (in pixels):");

            guiMaxText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiMaxText.setText("100");
            guiMaxText.setToolTipText("Radius of largest circle you expect to find");

            guiIncLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiIncLabel.setText("Radius search increment (in pixels):");
            guiIncLabel.setVisible(false);

            guiIncText.setVisible(false);
            guiIncText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiIncText.setText("1");
            guiIncText.setToolTipText("How much to increase the radius at each increment");

            guiMinNumLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiMinNumLabel.setText("Minimum number of circles to be found:");
            guiMinNumLabel.setVisible(false);

            guiMinNumText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiMinNumText.setText("0");
            guiMinNumText.setToolTipText("The minimum number of circles you expect to find in the image");
            guiMinNumText.setVisible(false);

            guiMaxNumLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiMaxNumLabel.setText("Maximum number of circles to be found:");
            guiMaxNumLabel.setVisible(false);

            guiMaxNumText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiMaxNumText.setText("65535");
            guiMaxNumText.setToolTipText("The maximum number of circles you expect to find in the image");
            guiMaxNumText.setVisible(false);

            guiThreshLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiThreshLabel.setText("Hough score threshold (0.0 - 1.0):");
            guiThreshLabel.setToolTipText("This is the ratio of votes/resolution");

            guiThreshText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiThreshText.setText("0.5");
            guiThreshText.setToolTipText("The Hough score threshold above which circles are valid");
            guiThreshText.setPreferredSize(new java.awt.Dimension(36, 20));

            guiResLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiResLabel.setText("Transform resolution (# of steps per transform):");
            guiResLabel.setToolTipText("The maximum number of transform steps per radius");
            guiResLabel.setVisible(false);

            guiResText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiResText.setText("1000");
            guiResText.setToolTipText("The number of steps to use per transform");
            guiResText.setVisible(false);

            guiClearLabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
            guiClearLabel.setText("Clear neighbors radius ratio:");
            guiClearLabel.setToolTipText("The ratio relative to the circle's radius");
            guiClearLabel.setVisible(false);

            guiClearText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiClearText.setText("1.0");
            guiClearText.setToolTipText("The radius around a found circle to clear out the Hough search space");
            guiClearText.setVisible(false);

            guiRadiusBandLabel.setVisible(false);
            guiRadiusBandLabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
            guiRadiusBandLabel.setText("Local radius search bandwidth (+/- previous radius):");

            guiRadiusBandText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiRadiusBandText.setText("10");
            guiRadiusBandText.setToolTipText("The range of radii to be used in the next frame");
            guiRadiusBandText.setVisible(false);

            guiSearchRadLabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
            guiSearchRadLabel.setText("Local radius search area (search radius near last centroid):");
            guiSearchRadLabel.setVisible(false);

            guiSearchRadText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiSearchRadText.setText("10");
            guiSearchRadText.setToolTipText("The area around a centroid to look for the next centroid");
            guiSearchRadText.setVisible(false);

            guiReduceBox.setSelected(true);
            guiReduceBox.setText("Reduce transform resolution by removing redundant transform steps");
            guiReduceBox.setVisible(false);

            guiLocalBox.setSelected(false);
            guiLocalBox.setText("Use a local search space to speed up search (only applies to movies)");
            guiLocalBox.addActionListener((java.awt.event.ActionEvent evt) -> {
                if(guiLocalBox.isSelected()){
                    //Setup local easy
                    if(guiEasyModeButton.isSelected()){
                        easyLocalGUI();
                    }
                    //Setup local advanced
                    else{
                        advancedLocalGUI();
                    }
                }
                else{
                    //Setup full easy
                    if(guiEasyModeButton.isSelected()){
                        easyFullGUI();
                    }    
                    //setup full advanced
                    else{
                        advancedFullGUI();
                    }
                }
            });

            guiOutputLabel.setFont(new java.awt.Font("Tahoma", 0, 14));// NOI18N
            guiOutputLabel.setText("Output options:");

            guiRawBox.setText("Raw Hough transform series");
            guiRawBox.setVisible(false);

            guiPointBox.setSelected(true);
            guiPointBox.setText("Circle outlines overlaid on the original image mask");

            guiIDBox.setText("Filled circles marked by ID number.");
            //guiIDBox.setVisible(false);               
                            
            guiHoughBox.setText("Filled circles marked by Hough score.");
            //guiHoughBox.setVisible(false);
            
            guiResultsBox.setSelected(true);
            guiResultsBox.setText("Export measurements to the results table");
            
            //UIManager.put("ProgressBar.background", Color.ORANGE); //Bar background color
            //UIManager.put("ProgressBar.foreground", Color.BLUE); //Bar color
            UIManager.put("ProgressBar.selectionBackground", Color.BLACK); //Font color above background
            UIManager.put("ProgressBar.selectionForeground", Color.BLACK); //Font color above bar
            guiProgressBar = new javax.swing.JProgressBar();
            guiProgressBar.setFont(new java.awt.Font("Courier New", 0, 11)); // NOI18N
            guiProgressBar.setVisible(false);
            guiProgressBar.setMaximum(100);
            guiProgressBar.setMinimum(0);
            guiProgressBar.setStringPainted(true); //Allow string to be written on progress bar


            guiOKButton.setText("OK");
            guiOKButton.addActionListener((java.awt.event.ActionEvent evt) -> { 
                // </editor-fold>
                // <editor-fold desc="Retrieve GUI arguments and calculate Hough parameters">
                
                //If the analysis is not started, begin the analysis
                if(!analysisStarted){
                    //Retrive the numbers from the text boxes and combobox
                    radiusMin = Integer.parseInt(guiMinText.getText());              
                    radiusMax = Integer.parseInt(guiMaxText.getText());
                    radiusInc = Integer.parseInt(guiIncText.getText());
                    minCircles = Integer.parseInt(guiMinNumText.getText());
                    maxCircles = Integer.parseInt(guiMaxNumText.getText());
                    thresholdRatio = Double.parseDouble(guiThreshText.getText());
                    resolution = Integer.parseInt(guiResText.getText());
                    ratio = Double.parseDouble(guiClearText.getText());
                    searchBand = Integer.parseInt(guiRadiusBandText.getText());
                    searchRadius = Integer.parseInt(guiSearchRadText.getText());
                    reduce = guiReduceBox.isSelected();
                    local = guiLocalBox.isSelected();

                    //Retrieve the check box status
                    houghSeries = guiRawBox.isSelected();
                    showCircles = guiPointBox.isSelected();
                    showID = guiIDBox.isSelected();
                    showScores = guiHoughBox.isSelected();
                    results = guiResultsBox.isSelected();

                    //Override searchBand and searchRad if in local easy mode
                    if(guiEasyModeButton.isSelected() & guiLocalBox.isSelected()){
                        searchRadius = radiusMin;
                        searchBand = radiusMax-radiusMin;
                    }

                    //Override impossible inputs
                    if (maxCircles > 65535) maxCircles = 65535;
                    if (minCircles > maxCircles) minCircles = maxCircles;
                    if (searchBand < 1) searchBand = 1;
                    if (searchRadius < 1) searchRadius = 1;

                    //Remove the GUI frame now that it is no longer needed
                    //guiFrame.dispose();

                    startTransform();
                    
                    //Set button to cancel
                    guiOKButton.setText("Cancel");
                    guiFrame.pack();
                    
                    analysisStarted = true;
                }
                
                //If analysis is running, cancel the analysis
                else{
                    guiInput.interruptThreads(true);
                    guiInput.cancel(true); //Stop analysis thread
                    IJ.showProgress(0); //Reset progress bar
                    guiOKButton.setText("OK");
                    guiFrame.pack(); 
                    analysisStarted = false;
                    IJ.showStatus("Analysis cancelled..."); //Update IJ status
                    guiInput = new Hough_Circle(); //Create new instance of analysis worker, since last worker thread was cancelled
                }    
            });         
            // </editor-fold>
            // <editor-fold desc="Swing GUI Part 2">
            javax.swing.GroupLayout layout = new javax.swing.GroupLayout(guiFrame.getContentPane());
            guiFrame.getContentPane().setLayout(layout);
            layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(guiTitle, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiIntro1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiIntro2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiHoughBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiIDBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiRawBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiPointBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiResultsBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiLocalBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiReduceBox, javax.swing.GroupLayout.DEFAULT_SIZE, 391, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(guiProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(guiOKButton))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(guiRadiusBandLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiSearchRadLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(guiSearchRadText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(guiRadiusBandText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(guiMaxNumLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiMaxLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiMinNumLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiMinLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiIncLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiResLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiThreshLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiClearLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(guiMinText)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                    .addComponent(guiMaxText, javax.swing.GroupLayout.DEFAULT_SIZE, 104, Short.MAX_VALUE)
                                    .addComponent(guiIncText)
                                    .addComponent(guiMinNumText, javax.swing.GroupLayout.Alignment.LEADING))
                                .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                    .addComponent(guiResText, javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addComponent(guiClearText, javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addComponent(guiMaxNumText)
                                    .addComponent(guiThreshText, javax.swing.GroupLayout.DEFAULT_SIZE, 104, Short.MAX_VALUE)))))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(88, 88, 88)
                                .addComponent(guiEasyModeButton)
                                .addGap(43, 43, 43)
                                .addComponent(guiAdvancedModeButton))
                            .addComponent(guiOutputLabel)
                            .addComponent(guiSearchLabel))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(guiTitle)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiIntro1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiIntro2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(guiAdvancedModeButton)
                    .addComponent(guiEasyModeButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiSearchLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiMinLabel)
                    .addComponent(guiMinText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiMaxLabel)
                    .addComponent(guiMaxText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiIncText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(guiIncLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiMinNumText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(guiMinNumLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(guiMaxNumText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(guiMaxNumLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 12, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(3, 3, 3)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiThreshText, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiThreshLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiResLabel)
                    .addComponent(guiResText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiClearLabel)
                    .addComponent(guiClearText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiRadiusBandLabel)
                    .addComponent(guiRadiusBandText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiSearchRadLabel)
                    .addComponent(guiSearchRadText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiReduceBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiLocalBox)
                .addGap(18, 18, 18)
                .addComponent(guiOutputLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(guiRawBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiPointBox, javax.swing.GroupLayout.PREFERRED_SIZE, 23, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiIDBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiHoughBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiResultsBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(guiOKButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(guiProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
                       
            //Use the Pack function to make the GUI frame shrink to the smallest size necassary
            guiFrame.pack();

            //Show the GUI
            guiFrame.setVisible(true);
            // </editor-fold>
        }
    }
    
    //Send the GUI values to the analysis class, and then run the analysis on a separate thread
    void startTransform(){
        //Create and instance of the analysis class
        guiInput = new Hough_Circle();
        
        KeyboardFocusManager manager = KeyboardFocusManager.getCurrentKeyboardFocusManager();
        manager.addKeyEventDispatcher(new MyDispatcher());
        
        //Add an action listener to the status to allow for GUI to be updated
        //Code modified from: http://www.javacreed.com/swing-worker-example/

        guiInput.addPropertyChangeListener((final PropertyChangeEvent event) -> {
            switch (event.getPropertyName()) {
                
                //If event has progress flag, update progress bar
                case "progress":
                    barString = guiInput.getStatus();
                    if(isGUI){
                        guiProgressBar.setIndeterminate(false);
                        guiProgressBar.setValue((Integer) event.getNewValue()); 
                        guiProgressBar.setString(barString);
                    }
                    break;
                    
                //if event has state flag, it indicates thread completion status
                case "state":
                    switch ((StateValue) event.getNewValue()) {
                        
                        //If worker thread is done, then clear out progress indicators and set button to "OK"
                        case DONE:
                            IJ.showProgress(0);
                            if(isGUI){
                                guiProgressBar.setVisible(false);
                                guiOKButton.setText("OK");
                            }
                            analysisStarted = false;
                            IJ.showStatus("Analysis complete...");
                            break;
                            
                        //If worker has just started, set progress to indetertminate, to let user know plugin is active
                        case STARTED:
                            if(isGUI){
                                guiProgressBar.setVisible(true);
                                guiProgressBar.setIndeterminate(true);
                                guiProgressBar.setString("Starting Transform...");
                            }
                            break;
                        case PENDING:
                            if(isGUI){
                                guiProgressBar.setVisible(true);
                                guiProgressBar.setIndeterminate(true);
                            }
                            break;
                    }
                    break;
            }
        });

        //Start the background transform by sending the GUI variables to the transform
        guiInput.setParameters(radiusMin, radiusMax, radiusInc, minCircles, maxCircles, thresholdRatio, resolution, ratio, searchBand, 
                searchRadius, reduce, local, houghSeries, showCircles, showID, showScores, results, isGUI);

        //Start the analysis on a separate thread so the GUI stays free.
        guiInput.execute();
    }
    
    void easyFullGUI(){
        guiTitle.setText("Hough Circle Transform");
        guiMinLabel.setVisible(true);
        guiMinText.setVisible(true);
        guiMaxLabel.setVisible(true);
        guiMaxText.setVisible(true);
        guiIncLabel.setVisible(false);
        guiIncText.setVisible(false);
        guiMinNumLabel.setVisible(false);
        guiMinNumText.setVisible(false);
        guiMaxNumLabel.setVisible(false);
        guiMaxNumText.setVisible(false);
        guiThreshLabel.setVisible(true);
        guiThreshText.setVisible(true);
        guiResLabel.setVisible(false);
        guiResText.setVisible(false);
        guiClearLabel.setVisible(false);
        guiClearText.setVisible(false);
        guiRadiusBandLabel.setVisible(false);
        guiRadiusBandText.setVisible(false);
        guiSearchRadLabel.setVisible(false);
        guiSearchRadText.setVisible(false);
        guiReduceBox.setVisible(false);
        guiRawBox.setVisible(false);
        //guiIDBox.setVisible(false);
        //guiHoughBox.setVisible(false);

        guiMaxNumText.setText("65535");
        //guiHoughBox.setSelected(false);
        //guiIDBox.setSelected(false); 
        guiRawBox.setSelected(false);
        guiFrame.pack();
    }
    void easyLocalGUI(){
        guiTitle.setText("Local Hough Circle Transform");
        guiMinLabel.setVisible(true);
        guiMinText.setVisible(true);
        guiMaxLabel.setVisible(true);
        guiMaxText.setVisible(true);
        guiIncLabel.setVisible(false);
        guiIncText.setVisible(false);
        guiMinNumLabel.setVisible(false);
        guiMinNumText.setVisible(false);
        guiMaxNumLabel.setVisible(true);
        guiMaxNumText.setVisible(true);
        guiThreshLabel.setVisible(true);
        guiThreshText.setVisible(true);
        guiResLabel.setVisible(false);
        guiResText.setVisible(false);
        guiClearLabel.setVisible(false);
        guiClearText.setVisible(false);
        guiRadiusBandLabel.setVisible(false);
        guiRadiusBandText.setVisible(false);
        guiSearchRadLabel.setVisible(false);
        guiSearchRadText.setVisible(false);
        guiReduceBox.setVisible(false);
        guiRawBox.setVisible(false);
        //guiIDBox.setVisible(false);
        //guiHoughBox.setVisible(false);

        //guiHoughBox.setSelected(false);
        //guiIDBox.setSelected(false); 
        guiRawBox.setSelected(false);
        guiMaxNumText.setText("10");
        guiMinNumText.setText("65535");
        guiFrame.pack();
    }
    void advancedFullGUI(){
        guiTitle.setText("Hough Circle Transform");
        guiMinLabel.setVisible(true);
        guiMinText.setVisible(true);
        guiMaxLabel.setVisible(true);
        guiMaxText.setVisible(true);
        guiIncLabel.setVisible(true);
        guiIncText.setVisible(true);
        guiMinNumLabel.setVisible(false);
        guiMinNumText.setVisible(false);
        guiMaxNumLabel.setVisible(true);
        guiMaxNumText.setVisible(true);
        guiThreshLabel.setVisible(true);
        guiThreshText.setVisible(true);
        guiResLabel.setVisible(true);
        guiResText.setVisible(true);
        guiClearLabel.setVisible(true);
        guiClearText.setVisible(true);
        guiRadiusBandLabel.setVisible(false);
        guiRadiusBandText.setVisible(false);
        guiSearchRadLabel.setVisible(false);
        guiSearchRadText.setVisible(false);
        guiReduceBox.setVisible(true);
        guiRawBox.setVisible(true);
        //guiIDBox.setVisible(true);
        //guiHoughBox.setVisible(true);

        guiMinNumText.setText("1");
        guiMaxNumText.setText("1");
        guiFrame.pack();
    }
    void advancedLocalGUI(){
        guiTitle.setText("Local Hough Circle Transform");
        guiMinLabel.setVisible(true);
        guiMinText.setVisible(true);
        guiMaxLabel.setVisible(true);
        guiMaxText.setVisible(true);
        guiIncLabel.setVisible(true);
        guiIncText.setVisible(true);
        guiMinNumLabel.setVisible(true);
        guiMinNumText.setVisible(true);
        guiMaxNumLabel.setVisible(true);
        guiMaxNumText.setVisible(true);
        guiThreshLabel.setVisible(true);
        guiThreshText.setVisible(true);
        guiResLabel.setVisible(true);
        guiResText.setVisible(true);
        guiClearLabel.setVisible(true);
        guiClearText.setVisible(true);
        guiRadiusBandLabel.setVisible(true);
        guiRadiusBandText.setVisible(true);
        guiSearchRadLabel.setVisible(true);
        guiSearchRadText.setVisible(true);
        guiReduceBox.setVisible(true);
        guiRawBox.setVisible(true);
        //guiIDBox.setVisible(true);
        //guiHoughBox.setVisible(true);

        guiMinNumText.setText(guiMaxNumText.getText());
        guiFrame.pack();
    }
    
    //If hey was pressed, then record if it was the escape key
    private class MyDispatcher implements KeyEventDispatcher {
        @Override
        public boolean dispatchKeyEvent(KeyEvent e) {
            switch (e.getID()) {
                case KeyEvent.KEY_PRESSED:
                    //Do nothing
                    break;
                case KeyEvent.KEY_RELEASED:
                    if(e.getKeyCode() == KeyEvent.VK_ESCAPE) guiInput.interruptThreads(true);
                    IJ.showProgress(0);
                    break;
                case KeyEvent.KEY_TYPED:
                    //Do nothing
                    break;
                default:
                    break;
            }
            return false;
        }
    }
}
