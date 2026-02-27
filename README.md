# Amsterdam Local Field Potential Analysis (ALFA) toolbox
The ALFA toolbox is a MATLAB tool developed for the extraction, processing, restructuring and visualization of data from JSON files containing local field potential data recorded by the Medtronic PerceptTM neurostimulator.

This is an open research tool that is not intended for clinical purposes.

## About
### Developed by:
* B.J. (Bart) Keulen       -   Amsterdam UMC, Amsterdam, The Netherlands
* M.J. (Mariëlle) Stam     -   Amsterdam UMC, Amsterdam, The Netherlands
* J.T. (Jackson) Boonstra  -   Vrije Universiteit Amsterdam, The Netherlands

### GitHub repository:
https://github.com/Data-Driven-Brain-Stimulation/ALFA-toolbox

### Please cite this article when using our toolbox:
Keulen BJ et al. Amsterdam Local Field potential Analysis (ALFA) toolbox: an open source software package for deep brain stimulation research 

## Outline
### Main scripts
The ALFA toolbox contains three main scripts. Every analysis or visualization step starts with running one of these scripts
* mainExtract.m     -   For the processing of JSON data files and visualization if needed. General settings are given manually in the section Settings within the script
* mainExtractUI.m   -   For the processing of JSON data files and visualization if needed. General settings are given through a user interface which appears when running the script
* mainPlot.m        -   For the visualization of output data of the mainExtract/mainExtractUI files

After running one of these scripts (and having set the general settings in the case of mainExtractUI), a window appears in which you are asked to select a file or folder (see settings dataset and folder).

### Settings
*mainExtract*
* dataset     -   Select either a a single file (0), a single folder (1) or a set of folders (2)
* ecgMethod   -   ECG suppression of Streaming data: set to (2) for (extensive) manual suppression, (1) for automatic suppression, (0) for no suppression
* rWindow     -   Window around R-peak [before after] for calculating SVD components, in seconds. [0.25 0.4] is default, alternatively use [0.2 0.2]
* linenoise   -   Frequency and bandwidth [Hz] of line-noise to remove, usually either 50 or 60 Hz with 0.2 or 0.5 Hz bandwidth
* tZone       -   Timezone to use for correcting for UTC offset and daylight saving time. Set manually by replacing with IANA timezone name (use command timezones() for an overview)
* plotData    -   Plot data (1) or not (0)
* showFig     -   Show figures (1) or not (0)

*mainExtractUI*
* Type of dataset                       -   Select either a a single file, a folder or a set of folders
* Create figures                        -   Check box to create figures of data
* Show figures                          -   Check box to show figures of data
* Timezone                              -   Select timezone to correct UTC offset and daylight saving time. Select local timezone of system, another timezone or no correction by selecting JSON default setting
* ECG suppression method                -   ECG suppression of Streaming data: set to (2) for (extensive) manual suppression, (1) for automatic suppression, (0) for no suppression
* R-peak window                         -   Window around R-peak for calculating SVD components, specified as time [ms] before and after each R-peak. 250 ms before and 400 ms after is default, alternatively use 200 ms before and 200 ms after. 
* Line noise filter center frequency    -   Select center frequency (50 or 60 Hz) of filter to remove line-noise (mains hum) from Streaming data. Default is 50 Hz.
* Line noise filter bandwidth           -   Select bandwidth of line-noise filter, default is 0.2 Hz. Alternatively use 0.5 Hz

*mainPlot*
* folder   -   Select either a single file (0) or a folder (1)
* showFig  -   Show figures (1) or not (0)

### Output
All output files are stored within a subfolder of the location of the selected file or folder. Depending on the type of dataset, this subfolder will be named:
* single file: [filename]_unpacked
* folder: folder_[name of folder]_unpacked
* folderset: folderset_[name of folderset]_unpacked

Depending on the type of data and type of dataset processing, the datafiles given as output are saved within this created subfolder as follows:

*Datatypes Setup, Survey, Identifier and Streaming*
•	single file: [filename]_[datatype].mat
•	folder: [name of folder]_[datatype].mat
•	folderset: [datatype]\[name of folder within folderset]_[number of JSON within folder]_[datatype].mat

*Datatypes Timeline and Events*
•	single file: [filename]_[datatype].mat
•	folder: [name of folder]_[datatype].mat
•	folderset: [datatype]\[name of folder within folderset]_[datatype].mat
