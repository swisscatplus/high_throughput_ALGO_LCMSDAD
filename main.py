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
import output as out
import optimization_script as opt_sc

starttime = time.time()
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
    "threshold_ms_spectra": .2,  # threshold MS spectra noise removal
    "ms_weighting_fct": "sin2",  # Weighting fct MS Data
    "threshold_ms_entropy_sim": 10.,
    "threshold_ms_dot_product_sim": .2,
    "threshold_ms_weighted_dp_sim": 1.,
    "threshold_ms_bhat1": 1.,
    "background folder": background_folder,
    "ms_algorithm": "dot_product",  # Options: all, dot_product, weighted_dot_product, bhat1, entropy_sim
    "dad_algorithm": "derivative",  # Options: all,
    "threshold_dad_pearson": 1.,
    "threshold_dad_dot_product": .000000000000000001,
    "threshold_dad_ft": 1.,
    "threshold_dad_derivative": .1,
    "FT_low_border": 0,  # Settings for cutting the FT spectrum before dot product; removes background; 0 is lowest
    "FT_upper_border": 53,  # To remove noise; 105 is highest (max wvl number)
    "nth_derivative": 4,  # Number of derivatives for comparison DAD.
    "retention_time_interval": 10.,  # Time interval to accept equal retention time. To be optimised.
    "directory_project": directory_project,  # Later change this everywhere for more convenience.
    "peak_folder_time": now.strftime("%Y-%m-%d_%H-%M-%S"),
    "ion_detection_mode": "positive",  # positive or negative
    "method_name": "AceticAcid01",
    "Number Columns": 8,  # Number for columns for heatmap visualization
}
print("Weighting function used for MS Spectra: " + settings["ms_weighting_fct"])
method_name = "AceticAcid01"  # Change into settings later!
background_method = "defaultNew"
print(time.time())



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
                                   "240611_synthese_Leander_2024-06-11_16-38-22+02-00-01.JSON")
# data_processing.create_new_background_spectra(background_filepath, background_method, settings)




"""
---------------------------------------------------------
To create new database entry from an existing peak file, enter peak file path here and uncomment the function. Enter
a molecule name. Should the function fail with the molecule name, enter the molecule smiles.
---------------------------------------------------------
"""
peak_path = os.path.join(directory_project, "Data_examples", "Peak_files", "2024-04-12_14-38-45", "peak_2.cdf")
# molecule_name = "Indoline"
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

# extracted_ms.full_processing_ms(background_signal, settings)

# data_processing.process_found_signal(full_analysis, peak_info, settings)
# sp_handling.signal_comparison(settings)
# dad_chr = dpr_dad.DAD_full_chr(full_analysis.dad_data3d, full_analysis.info)

# dadms.assign_peaks(full_analysis, settings)
# dad_comp.comparison_dad(extracted_dad, extracted_dad2, settings)

run_name = "240611_synthese_Leander_2024-06-11_17-29-04+02-00-04.JSON"
# runs.analyse_single_run(run_name, method_name, background_method, settings)

run_folder_name = "screening_suzuki"
# runs.analyse_multiple_runs(run_folder_name, method_name, background_method, settings)

"""
Testing superimposed ms peak deco
"""
# runs.delete_old_sgn_files(settings)
json_path = os.path.join(directory_project, "Data_examples", "testfiles", run_name)

# full_analysis = init.import_run_json(json_path, method=background_method)
# ms_chr = data_processing.MS_full_chr(full_analysis.ms_data3d, full_analysis.info)
# ms_peaks = msd.ms_create_peaks(full_analysis, ms_chr, settings)
# ms_spectra = ms_chr.extract_single_ms(330)
# ms_spectra.plot_spectra()

# out.dtb_molecule_list(settings)
# out.dtb_molecule_full_data("RCIJACVHOIKRAP-UHFFFAOYSA-M.cdf", settings)
# out.create_analysis_report(settings, run_folder_name, report_name="screening_suzuki_positive2")  # peak_folder="2024-04-22_15-43-35"
# 2024-05-22_10-00-49
"""opt_sc.comparison_dtb_named_files(settings, True)
print("Done with processed")
print(time.time())
opt_sc.comparison_dtb_named_files(settings, False)"""

# opt_sc.plot_optimization_dad("dad_optimization_peaks_dtb_compTrue.csv", settings)
opt_sc.plot_optimization_ms("ms_optimization_peaks_dtb_compFalse.csv", settings)


"""peak_file_path = os.path.join(directory_project, "Data_examples", "Peak_files", "2024-06-14_16-37-50chrysin",
                              "peak_4.cdf")
out.plot_peak_data(peak_file_path, True)

molecule_name = "Benzaldehyde"
#create_file.create_dtb_entry(directory_project, peak_file_path, molecule_name)"""




"""
Analyse molecules for dtb:
"""


"""timestamp_peakfolder = settings["peak_folder_time"]"""



"""
---------------------------------------------------------------------------------------


molecule_folder = "cloro-nitropyridine"
background_file = "240419-screening_C18_acid_acetic-C-37.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)

---------------------------------------------------------------------------------------
"""





print(time.time())
print(time.time()-starttime)
