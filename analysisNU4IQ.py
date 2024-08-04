#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to analyze the NEMA NU4 image(s) quality phantom. The analysis follows the NEMA NU4 protocol and the reported 
results are
1) Uniformity
2) Spill-over ratios (water and air cavities)
3) Recovery coefficients.

The script assumes that the VOIs were created and placed in AMIDE beforehand, and that VOI (volume of interest) 
statistics were saved to a .tsv file. A .xif template is normally available together with the current script so that 
the VOI can be easily placed in AMIDE. The script can handle multiple images, for instance reconstructed at different 
numbers of iterations. The VOIs statistics in the .tsv file can therefore be for multiple datasets. The final report 
will give the results for all the datasets.

Requirements:
	- VOIs are named: Unif, Water, Air, Rod_{1..5} from AMIDE (use the template!)
	- The VOI statistics are saved in a .tsv file in AMIDE (use the template!)

Summary of the NEMA NU4 protocol:
	- Uniformity: VOI with a diameter of 22.5 mm and a length of 10 mm
		- Report: Mean, Max, Min and %STD
	- Spill-over ratios: Two VOIs with a diameter of 4 mm and a length of 7.5 mm (air and water)
		- Report: ratio of the mean activity of the VOIs and the mean activity in the uniform VOI.
	- Recovery coefficients: Five VOIs with a diameter twice the rod diameter and length of 10 mm.
		- Averaging the voxels along the axial axis and determining the maximum average activity inside each VOI.
		- Report: 
			- Mean: ratio of the maximum average and the mean activity in the hot uniform VOI
			- Stdev: standard deviation of the activity along the axial axis at the position of the maximum average 
			  activity.

Notes:
	- SOR: Spill-over ratio
	- RC: Recovery coefficient
	- Coefficient of variation: Stdev / Mean

Usage example:
	python analysisNU4IQ.py -i roi_results.tsv

TODO: 
	- Clean up and check the zxc 
	- Add some checks
"""



########################################################################################################################
# Imports
########################################################################################################################
import re
import argparse
import numpy as np 
from collections import defaultdict
from typing import Dict, Union, List


########################################################################################################################
# Constants
########################################################################################################################
ROD_NAMES = ["Rod_1", "Rod_2", "Rod_3", "Rod_4", "Rod_5"]
CAVITY_NAMES = ["Water", "Air"]



########################################################################################################################
# Methods created to build the features.
########################################################################################################################
def parseAmideRawMeasurementFile(_file_path: str):
	'''
	Def: Parses the raw measurement file saved by Amide of the IQ phantom image. The function saves the information in 
	     a two-level dictionary where each key is the image name and the second level is the VOI name and its 
	     corresponding raw measurement.
	@_file_path: Name of the file where the Amide raw measurements of the VOIs were saved. 
	Returns:
		{images}{VOI}: Two-level dictionary containing the image(s) and VOIs.
	'''
	# Initialize the two-level dictionary
	data_dict = defaultdict(lambda: defaultdict(list))
		
	id = 0
	with open(_file_path, 'r') as file:
		current_data_set = None
		current_roi = None
		current_table = []
		for line in file:
			if line.startswith("#"):
				# Attempt to extract the data set name in the current line 
				data_set_match = re.match(r"#\s+Data Set:\s+(.+?)\s+Scaling Factor:", line)
				# Attempt to extract the ROI name in the current line 
				roi_match = re.match(r"#\s+ROI:\s+(\S+)\s+Type:", line)
				if data_set_match:
					# The data is provided as outer loop of VOI and inner loop of dataset.
					# Thus, each time a dataset name is detected, data follow and we need to prepare for new data.
					# Store the data accumulated in the dictionary (if not the first in both loop)
					if (current_data_set is not None) and (id != 0):
						data_dict[current_data_set][current_roi] = np.array(current_table, dtype=float)
					# Start a new one
					current_table = []
					# zxc provide other choice of IDs
					current_data_set = str(id) 
					id +=1
				elif roi_match:
					# Store the current table in the dictionary before starting a new one
					# zxc revoir si nec!!! aussi manque le current_table = [] might be okay 
					# due to data_set_match always triggering after
					if current_roi is not None:
						data_dict[current_data_set][current_roi] = np.array(current_table, dtype=float)
						id = 0
					current_roi = roi_match.group(1)
			elif line.strip():
				# Add non-comment lines to the current table
				current_table.append(line.strip().split('\t'))
		
		# Store the last table after exiting the loop
		if current_data_set and current_roi:
			data_dict[current_data_set][current_roi] = np.array(current_table, dtype=float)

	return data_dict


def extractBasicMetrics(_voiInfo: np.ndarray) -> Union[float, float, float, float]:
	'''
	Def: Extract the basic metrics (mean, std dev, min and max) from the VOIs.
	@_voiInfo[:, 5]: Array describing the pixels inside the current VOI. The columns are the pixel value, the fraction 
	                 inside the VOI and the (x, y, z) position of the pixel.
	Returns: 
	    Mean, StdDev, Min and Max of the VOI.
	'''
	# Can handle fractions of pixels
	pixelVal = _voiInfo[:, 0]
	pixelFrac = _voiInfo[:, 1]

    # Extract the basic statistics
	mean = np.sum(pixelVal * pixelFrac) / np.sum(pixelFrac)
	stdDev = np.sqrt(np.sum(pixelFrac * (pixelVal - mean)**2) / np.sum(pixelFrac))
	min = np.min(pixelVal)
	max = np.max(pixelVal)

	return mean, stdDev, min, max


def extractRodLineProfile(_voiInfo: np.ndarray) -> np.ndarray:
	'''
	Def: According to NU4 protocol, we have to average the voxels along the axial axis for the VOIs of the five hot 
	     rods. We then determine the maximum average activity. At the voxel position having the maximum average 
	     intensity, we extract the intensity values at each slice along the axial profile. These values can then be 
	     used to compute the standard deviation of the recovery coefficient in another function. This function extracts 
	     the line profile at the maximum average voxel position.
	@_voiInfo[:, 5]: Array describing the pixels inside the current VOI. The column are the pixel value, the fraction 
	                 inside the VOI and the (x, y, z) position of the pixel.
	Returns:
		Axial line profile at the maximum average voxel position.
	'''
	# Extract the X, Y position then the intensity
	voxelXPos = _voiInfo[:, 2]
	voxelYPos = _voiInfo[:, 3]
	voxelVal = _voiInfo[:, 0]

	# Create axial line profile at each X,Y position
	axialLineProfile = {}
	for i in range(_voiInfo.shape[0]):
		xyKey = str(voxelXPos[i]) + "_" + str(voxelYPos[i])
		if xyKey in axialLineProfile:
			axialLineProfile[xyKey].append(voxelVal[i])
		else:
			axialLineProfile[xyKey] = [voxelVal[i]]

	# Find the axial profile at maximum average pixel position
	maxAvgKey = None 
	maxAvg = 0.0
	for key in axialLineProfile:
		currAvg = np.mean(axialLineProfile[key])
		if maxAvg < currAvg:
			maxAvgKey = key 
			maxAvg = currAvg 

	return axialLineProfile[maxAvgKey]


def evalUniformity(_vois: Dict[str, np.ndarray]) -> Dict[str, float]:
	''' 
	Def: Evaluate uniformity from the uniformity VOI.
	@_vois: Dictionary of VOIs.
	Returns: 
		Dictionary of uniformity results
	'''
	# Evaluate mean, min and max of uniformity VOI, then coefficient of variation
	unifMean, unifStdDev, unifMin, unifMax = extractBasicMetrics(_vois['Unif'])
	unifCoeffOfVar = unifStdDev / unifMean
	unifResults = {"mean": unifMean, "min": unifMin, "max": unifMax, "unifCoeffOfVar": unifCoeffOfVar}
	return unifResults


def evalSpillOverRatios(_vois: Dict[str, np.ndarray], _unifResults: Dict[str, float]) -> Dict[str, Dict[str, float]]:
	''' 
	Def: Evaluate spill-over ratios of water and air cavities.
	@_vois: dictionary of VOIs
	@_unifResults: dictionary of uniformity results
	Returns: 
		Dictionary of spill-over ratio results
	'''
	# SOR: Spill Over Ratio
	sorResults = {}

	# Loop on the two cavities
	for cavity in CAVITY_NAMES:
		# Extract mean and stdev of current material VOI
		mean, stdev, _, _ = extractBasicMetrics(_vois[cavity])

		# Evaluate SOR and error
		sor = mean / _unifResults["mean"]
		sorCoeffOfVar = stdev / mean
		sorError = sor * np.sqrt(sorCoeffOfVar**2 + _unifResults["unifCoeffOfVar"]**2)
		sorResults[cavity] = {"sor": sor, "sorError": sorError}

	return sorResults


def evalRecovCoeff(_vois: Dict[str, np.ndarray], _unifResults) -> Dict[str, Dict[str, float]]:
	''' 
	Def: Evaluate recovery coefficients of the five rods.
	@_vois: Dictionary of VOIs
	@_unifResults: Dictionary of uniformity results
	Returns: 
		Dictionary of recovery coefficients results
	'''
	# RC: Recovery coefficient
	rcResults = {}

	# Loop on the five rods
	for cRod in ROD_NAMES:
		# Extract the line profile at the maximum average position
		maxAvglineProfile = extractRodLineProfile(_vois[cRod])

		# Evaluate recovery coefficient and error
		cRodRecCoeff = np.mean(maxAvglineProfile) / _unifResults["mean"]
		cRodCoeffOfVar = np.std(maxAvglineProfile) / np.mean(maxAvglineProfile)
		cRodError = cRodRecCoeff * np.sqrt(cRodCoeffOfVar**2 + _unifResults["unifCoeffOfVar"]**2)
		rcResults[cRod] = {"rc": cRodRecCoeff, "rcError": cRodError}

	return rcResults


def showReport(_unifResults: Dict[str, float], _sorResults: Dict[str, Dict[str, float]], \
	           _rcResults: Dict[str, Dict[str, float]]) -> None:
	''' 
	Def: Show a report of the analysis results.
	@_unifResults: Dictionary with uniformity results.
	@_sorResults: Dictionary with spill-over ratios results.
	@_rcResults: Dictionary with recovery coefficients results.
	'''
	print("\n*******************************")
	print("==== Uniformity ====")
	print(f"  Mean: {_unifResults['mean']:.4f}")
	print(f"   Min: {_unifResults['min']:.4f}")
	print(f"   Max: {_unifResults['max']:.4f}")
	print(f"  %STD: {100.0 * _unifResults['unifCoeffOfVar']:.2f}%")

	print("\n==== Spill-over ratios ====")
	for cavity in CAVITY_NAMES:
		align = "  " if cavity == "Air" else ""
		print(f" {align}{cavity}: {_sorResults[cavity]['sor']:.3f} ± {_sorResults[cavity]['sorError']:.3f}")

	print("\n==== Recovery coefficients ====")
	for cRod in ROD_NAMES:
		print(f"  {cRod.replace('Rod_', '')} mm: {_rcResults[cRod]['rc']:.2f} ± {_rcResults[cRod]['rcError']:.2f}")
	print("************************************************************\n")



########################################################################################################################
# List of methods used to use this module as a script.
########################################################################################################################
def main(args) -> None:
	''' Main function to execute the script. '''
	# Load .tsv data (VOI statistics) obtained from AMIDE
	data = parseAmideRawMeasurementFile(args.iFile)

	# Loop on each dataset (can be one or more)
	for data_set, vois in data.items():
		print(f"The current dataset is named {data_set}")

		# Uniformity, SOR and Recovery coefficients
		unifResults = evalUniformity(vois)
		sorResults = evalSpillOverRatios(vois, unifResults)
		rcResults = evalRecovCoeff(vois, unifResults)

		# Show report of results
		showReport(unifResults, sorResults, rcResults)


def readArguments():
	'''
	Def: Parse command line arguments.
	Returns: Provided arguments.
	'''
	parser = argparse.ArgumentParser(description=("Compute Nema Nu4 metrics from the raw measurements of the VOIs " \
	                                              "extracted with Amide."))
	
	# Basic:
	parser.add_argument('-i', '--iFile', action='store', type=str, required=True, dest='iFile',  
	                    help=("File where the measurements of Amide were saved."))

	# Features:
	# n/a 

	return parser.parse_args()



########################################################################################################################
# Main : Generate a empty config file or validate an existing configuration file.
########################################################################################################################
if __name__=='__main__':
	args = readArguments()
	main(args)