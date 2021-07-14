# TIRF-SIM
Supplementary Matlab code for processing TIRF-SIM data using CME Analysis package

## Description
This Matlab code can be used to process TIRF-SIM microscopic data using cmeAnalysis package: https://github.com/DanuserLab/cmeAnalysis. This is in fact not directly possible due to signifficant data differences (resolution, amount of noise, individual CCP characteristics, etc.). These tools make it partially possible by creating (emulating) TIRF data (non SIM) which is processed by cmeAnalysis in a standard way. The results (found tracks) are then visualized together with reconstructed SIM data and it is possible to manually select *interesting* tracks (in any sense), optionally polish the data by precise centering, changing the time position of the track (begin / end) and export such tracks as separate tiff movies for further (statistical) processing.

## Prerequisities
* Matlab with working cmeAnalysis package
* TIRF-SIM data in *mrc* format (both *raw* and *reconstructed* data are needed)

## Usage workflow
### 1. Prepare data for cmeAnalysis
Use *mrcMean2tiff.m* to convert raw SIM data from *mrc* format to *tiff* while using only one of the 9 SIM images or averaging them all. This is the easiest way to generate data having TIRF like characteristics needed by cmeAnalysis.
Note, that for reading *mrc* data, *ReadMRC.m* function is used, which was taken from [MathWorks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/27021-imagic-mrc-dm-and-star-file-i-o).
### 2. Run cmeAnalysis
Refer to cmeAnalysis documentation on how to process the data, only detection and tracking need to be run. Nevertheless, analysis could be run to see whether the CCP characteristics are meaningfull (i.e. that the cmeAnalysis didn't fail).
### 3. View and select interesting tracks
Use *tirfSimGui.m* function to visualize the detected tracks on the reconstructed SIM data. For each track meant for further processing click the 'Export' button at the top. This also saves *tiff* movie for this track as is detected from cmeAnalysis.
### 4. Improve CCP centering and track time bounds, export final tracks
For tracks exported by *tirfSimGui.m* use *tirfSimCentering.m* for precise manual adjustment of CCP centers and track time bounds (begin / end). You can see two graphs for CCP evaluatin: CCP intensity time behaviour (top graph) and mean intensity based on distance from the center ("radial intensity", bottom graph). Finally export all the tracks to *tiff* movies using the 'Resave all cut-outs' button.
