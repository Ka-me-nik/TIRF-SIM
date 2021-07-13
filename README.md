# TIRF-SIM
Supplementary Matlab code for processing TIRF-SIM data using CME Analysis package

## Description
This Matlab code can be used to process TIRF-SIM microscopic data using cmeAnalysis package: https://github.com/DanuserLab/cmeAnalysis. This is in fact not directly possible due to signifficant data differences (resolution, amount of noise, individual CCP characteristics, etc.). These tools make it partially possible by creating (emulating) TIRF data (non SIM) which is processed by cmeAnalysis in a standard way. The results (found tracks) are then visualized together with reconstructed SIM data and it is possible to manually select *interesting* tracks (in any sense), optionally polish the data by precise centering, changing the time position of the track (begin / end) and export such tracks as separate tiff movies for further (statistical) processing.

## Prerequisities
* Matlab with working cmeAnalysis package
* TIRF-SIM data in *mrc* format (both *raw* and *reconstructed* data are needed)

## Usage workflow
1. Prepare data for cmeAnalysis
2. Run cmeAnalysis
3. View and select interesting tracks
4. Improve CCP centering and track time bounds
5. Export final tracks
