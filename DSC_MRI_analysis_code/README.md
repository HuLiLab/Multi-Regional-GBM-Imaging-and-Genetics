Author: Nate Semmineh
This zip folder consists of 4 files
	* dR2s_NEU_GPM.m - MATLAB code
	* README.md - A readme file
	* ROI_Coordinates.xlsx - Excel sheet
	* Sample data - Demo data


# dR2s_NEU_GPM.m

	A MATLAB code to compute the T2* weighted relativity changes (delta R2*) on tumor ROI and whole brain non enhancing (WBNE) voxels.

	# System Requirements

		* Any operating system having MATLAB R2018a or higher versions

	# Data Requirements

		* DSC MRI perfusion data - NIFTI format
		* Excel sheet containing the x,y and z coordinates of the tumor ROI

	# Code Description

		The code takes out the DSC MRI perfusion data from the subject folder. For each subject, a brain mask is created to perform skull stripping. 
		
		The steady-state and gradient timepoints were computed along with the peak of the transient dip on the T2* signal curve.

		The WBNE mask is computed by considering those voxels which has a transient dip (more than 6SD from mean baseline) and return to baseline (mean of final 10 seconds of data is +/- 2SD from mean baseline)

		The tumor coordinates (x,y,z) for the current subject are extracted from the excel sheet given, creating a 10x10 tumor ROI. The mean T2* relaxivity curve (dR2* curve) is then calculated for the tumor ROI and WBNE map.

		For each glycolytic and neuronal group, the mean dR2* for the tumor and WBNE for all the subject from each group was plotted.

	# Output

		The code outputs a plot with the dR2s curve for the tumor ROI and WBNE. The runtime for the demo is ~4 minutes.

# ROI_Coordinates.xls

	This excel sheet consists of the tumor ROI coordinates for each subject.

# Sample_data

	This folder consists of a sample data for the demo. The MCH* folder has the DSC MRI perfusion data in NIFTI format  (2 files: .nii.gz and .json) 

