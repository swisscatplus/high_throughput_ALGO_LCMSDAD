import initialize as init
import os
import time

import data_processing
import databank_handling as dtb
import create_file
import ms_spectra_comparison
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter, find_peaks_cwt, find_peaks
import scipy.stats
import signal_peak_handling as sp_handling
import data_processing_dad as dpr_dad
import dad_spectra_comparison as dad_comp
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.datasets import make_blobs
import ms_deconvolution as msd
import pywt
from datetime import datetime
import combine_ms_dad as dadms
import runs_handling as runs


# Go back directories until we have project folder
current_directory = os.getcwd()
directory_project = os.path.abspath(os.path.join(current_directory, os.pardir))
print(directory_project)

# Go to background folder
background_folder = os.path.join(directory_project, "Data_examples", "background_spectra")

# define time stamp for peak folder name
now = datetime.now()

"""
---------------------------------------------------------
Define default settings
(maybe more to come)
---------------------------------------------------------
"""
settings = {
    "threshold_ms_spectra": .001,  # threshold MS spectra noise removal
    "ms_weighting_fct": "None",  # Weighting fct MS Data
    "threshold_ms_entropy_sim": 10.,
    "threshold_ms_dot_product_sim": 1.,
    "threshold_ms_weighted_dp_sim": 1.,
    "threshold_ms_bhat1": 1.,
    "User_Interface": False,
    "background folder": background_folder,
    "ms_algorithm": "all",  # Options: all, dot_product, weighted_dot_product, bhat1, entropy_sim
    "dad_algorithm": "all",  # Options: all,
    "threshold_dad_pearson": 1.,
    "threshold_dad_dot_product": 1.,
    "threshold_dad_ft": 1.,
    "threshold_dad_derivative": 1.,
    "FT_low_border": 0,  # Settings for cutting the FT spectrum before dot product; removes background; 0 is lowest
    "FT_upper_border": 53,  # To remove noise; 105 is highest (max wvl number)
    "nth_derivative": 1,  # Number of derivatives for comparison DAD.
    "retention_time_interval": .5,  # Time interval to accept equal retention time. To be optimised.
    "directory_project": directory_project,  # Later change this everywhere for more convenience.
    "peak_folder_time": now.strftime("%Y-%m-%d_%H-%M-%S"),
    "ion_detection_mode": "positive"
}
print("Weighting function used for MS Spectra: " + settings["ms_weighting_fct"])
method_name = "AceticAcid01"
background_method = "defaultNew"



"""
---------------------------------------------------------
This is the path to the analytical file directly from OpenLab.
---------------------------------------------------------
"""
directory_data = os.path.join(directory_project, "Data_examples", "20230907 145900")
testpath = os.path.join(directory_data, "20230907 145900_MS1 +TIC SCAN ESI Frag=100V Gain=1-0_spectra.cdf")
# Should go through a folder of multiple runs later.




"""
---------------------------------------------------------
To create new background file, enter path here, change the name and uncomment the function.
---------------------------------------------------------
"""
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   "240405_Acetic_acid_HILIC-04.JSON")
data_processing.create_new_background_spectra(background_filepath, method_name, settings)




"""
---------------------------------------------------------
To create new database entry from an existing peak file, enter peak file path here and uncomment the function. Enter
a molecule name. Should the function fail with the molecule name, enter the molecule smiles.
---------------------------------------------------------
"""
peak_path = os.path.join(directory_project, "Data_examples", "Peak_files", "2024-04-12_14-38-45", "peak_2.cdf")
molecule_name = "Indoline"
# create_file.create_dtb_entry(directory_project, peak_path, molecule_name)


"""
---------------------------------------------------------
To add method and retention time to existing database file, enter peak path and database file path here, then
uncomment the function.
---------------------------------------------------------
"""
# peak_path = os.path.join(directory_project, "Data_examples", "Peak_files", "2024-03-07_18-04-27", "peak_1.cdf")

dtb_entry_path = os.path.join(directory_project, "Data_examples", "database", "HJKGBRPNSJADMB-UHFFFAOYSA-N.cdf")
# create_file.add_method_to_dtb_entry(dtb_entry_path, peak_path)




"""
---------------------------------------------------------
Deal with peaks to be compared or whatever. Make sure to enter folder time stamp of latest change.
---------------------------------------------------------
"""
# Folder time stamp will later be automatised for a direct full comparison...

peak_directory = os.path.join(directory_project, "Data_examples", "Peak_files", "2024-03-14_15-30-45")

# dtb.compare_all_peaks_dtb(peak_directory, settings)



"""
---------------------------------------------------------
Testing
---------------------------------------------------------
"""
dad_path = os.path.join(directory_project, "Data_examples", "test.txt")
# init.import_dad_spectrum(dad_path)


"""
Testing for import .json file
"""
"""json_path = os.path.join(directory_project, "Data_examples", "testfiles", "carboxylic_acid_and_amines_methode_type_poroshell_7-01.JSON")
full_analysis = init.import_run_json(json_path, method = method_name)
peak_info = {
    "start_time": 177.15,
    "end_time": 190.325,
}

background_signal = data_processing.load_background_spectra_signal(full_analysis.info["LC Method"], settings)  # later give as input, requires preprocessed dtb
ms_chr = data_processing.MS_full_chr(full_analysis.ms_data3d, full_analysis.info)
extracted_ms = ms_chr.extract_ms_timespan(peak_info["start_time"], peak_info["end_time"])"""
# extracted_ms.full_processing_ms(background_signal, settings)

# data_processing.process_found_signal(full_analysis, peak_info, settings)
# sp_handling.signal_comparison(settings)
# dad_chr = dpr_dad.DAD_full_chr(full_analysis.dad_data3d, full_analysis.info)

# dadms.assign_peaks(full_analysis, settings)
# dad_comp.comparison_dad(extracted_dad, extracted_dad2, settings)

run_name = "240405_Acetic_acid_HILIC-05.JSON"
runs.analyse_single_run(run_name, method_name, settings)

run_folder_name = "example01"
# runs.analyse_multiple_runs(run_folder_name, method_name, settings)
