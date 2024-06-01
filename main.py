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
    "threshold_ms_spectra": .001,  # threshold MS spectra noise removal
    "ms_weighting_fct": "None",  # Weighting fct MS Data
    "threshold_ms_entropy_sim": 10.,
    "threshold_ms_dot_product_sim": .000000000000000001,
    "threshold_ms_weighted_dp_sim": 1.,
    "threshold_ms_bhat1": 1.,
    "background folder": background_folder,
    "ms_algorithm": "dot_product",  # Options: all, dot_product, weighted_dot_product, bhat1, entropy_sim
    "dad_algorithm": "dot_product",  # Options: all,
    "threshold_dad_pearson": 1.,
    "threshold_dad_dot_product": .000000000000000001,
    "threshold_dad_ft": 1.,
    "threshold_dad_derivative": 1.,
    "FT_low_border": 0,  # Settings for cutting the FT spectrum before dot product; removes background; 0 is lowest
    "FT_upper_border": 53,  # To remove noise; 105 is highest (max wvl number)
    "nth_derivative": 1,  # Number of derivatives for comparison DAD.
    "retention_time_interval": 5.,  # Time interval to accept equal retention time. To be optimised.
    "directory_project": directory_project,  # Later change this everywhere for more convenience.
    "peak_folder_time": now.strftime("%Y-%m-%d_%H-%M-%S"),
    "ion_detection_mode": "positive",  # positive or negative
    "method_name": "AceticAcid01",
    "Number Columns": 5,  # Number for columns for heatmap visualization
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
                                   "screening-molecules_and-melange-Leander_2-27.JSON")
# data_processing.create_new_background_spectra(background_filepath, background_method, settings)




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

# extracted_ms.full_processing_ms(background_signal, settings)

# data_processing.process_found_signal(full_analysis, peak_info, settings)
# sp_handling.signal_comparison(settings)
# dad_chr = dpr_dad.DAD_full_chr(full_analysis.dad_data3d, full_analysis.info)

# dadms.assign_peaks(full_analysis, settings)
# dad_comp.comparison_dad(extracted_dad, extracted_dad2, settings)

run_name = "screening-molecules_and-melange-Leander_2-06.JSON"
# runs.analyse_single_run(run_name, method_name, background_method, settings)

run_folder_name = "testnmf"
# runs.analyse_multiple_runs(run_folder_name, method_name, background_method, settings)

"""
Testing superimposed ms peak deco
"""
"""runs.delete_old_sgn_files(settings)
json_path = os.path.join(directory_project, "Data_examples", "testfiles", run_name)

full_analysis = init.import_run_json(json_path, method=background_method)
ms_chr = data_processing.MS_full_chr(full_analysis.ms_data3d, full_analysis.info)
ms_peaks = msd.ms_create_peaks(full_analysis, ms_chr, settings)
for peak in ms_peaks:
    print(peak)"""

# out.dtb_molecule_list(settings)
# out.dtb_molecule_full_data("VAOCPAMSLUNLGC-UHFFFAOYSA-N.cdf", settings)
#out.create_analysis_report(settings, run_folder_name, report_name="testnewexport")  # peak_folder="2024-04-22_15-43-35"
# 2024-05-22_10-00-49
# opt_sc.comparison_dtb_named_files(settings)
# opt_sc.plot_optimization_dad("dad_optimization_peaks_dtb_comp.csv", settings)
# opt_sc.plot_optimization_ms("ms_optimization_peaks_dtb_comp.csv", settings)

timestamp_peakfolder = settings["peak_folder_time"]
"""
Analyse molecules for dtb:
"""
molecule_folder = "cis-dimethyl-dioxane-dione"
background_file = "240418-screening_C18_acid_acetic-C-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "pyrrolidone-dione"
background_file = "240418-screening_C18_acid_acetic-C-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethyl-glutanimide"
background_file = "240418-screening_C18_acid_acetic-C-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""

molecule_folder = "dimethoxy-benzil"
background_file = "240418-screening_C18_acid_acetic-C-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "ditrimethylolpropane"
background_file = "240418-screening_C18_acid_acetic-C-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dihydroxy-diphenyl-methane"
background_file = "240418-screening_C18_acid_acetic-C-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "carboxymethyl-cyclopentane-tricarboxylic-acid"
background_file = "240418-screening_C18_acid_acetic-C-16.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethoxy-acetophenone"
background_file = "240418-screening_C18_acid_acetic-C-16.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "tetraclorohydroquinone"
background_file = "240418-screening_C18_acid_acetic-C-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "bromo-methy-h-imidazole"
background_file = "240418-screening_C18_acid_acetic-C-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "phloroglucinol"
background_file = "240418-screening_C18_acid_acetic-C-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "bispentafluorophenylcarbonate"
background_file = "240418-screening_C18_acid_acetic-C-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "tetrachloro-nitrobenzene"
background_file = "240418-screening_C18_acid_acetic-C-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "terephthaldehyde"
background_file = "240418-screening_C18_acid_acetic-C-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzyloxy-methoxy-propenylbenzene"
background_file = "240418-screening_C18_acid_acetic-C-16.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "triphenylmethanol"
background_file = "240418-screening_C18_acid_acetic-C-16.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "bicyclohexane-diol"
background_file = "240418-screening_C18_acid_acetic-C-16.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethoxyphenol"
background_file = "240418-screening_C18_acid_acetic-C-23.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzin"
background_file = "240418-screening_C18_acid_acetic-C-23.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "pentoerythritol"
background_file = "240418-screening_C18_acid_acetic-C-23.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "pyrazole"
background_file = "240418-screening_C18_acid_acetic-C-30.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "phenyl-acetophenone"
background_file = "240418-screening_C18_acid_acetic-C-30.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "curcumin"
background_file = "240418-screening_C18_acid_acetic-C-30.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "crown-ether"
background_file = "240418-screening_C18_acid_acetic-C-37.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "cloro-nitropyridine"
background_file = "240418-screening_C18_acid_acetic-C-37.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethyl-terephthalate"
background_file = "240418-screening_C18_acid_acetic-C-37.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dioctyl-sulfo-succinate"
background_file = "240418-screening_C18_acid_acetic-C-44.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "cyclohexanediaceticanhydride"
background_file = "240418-screening_C18_acid_acetic-C-44.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "diphenanesulfane"
background_file = "240418-screening_C18_acid_acetic-C-44.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dlbenzyl"
background_file = "240418-screening_C18_acid_acetic-C-51.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethane-sulfonamide"
background_file = "240418-screening_C18_acid_acetic-C-51.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "cyclohexylphenol"
background_file = "240418-screening_C18_acid_acetic-C-51.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "butylsulfone"
background_file = "240418-screening_C18_acid_acetic-C-58.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "hydroxyquinoline"
background_file = "240418-screening_C18_acid_acetic-C-58.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "triphenylphosphoniumchloride"
background_file = "240418-screening_C18_acid_acetic-C-58.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "methyltriphenylphosphoniumbromide"
background_file = "240418-screening_C18_acid_acetic-C-65.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "picodine-n-oxide"
background_file = "240418-screening_C18_acid_acetic-C-65.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "piperidine-dione"
background_file = "240418-screening_C18_acid_acetic-C-65.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "tetraamine-palladiumchloride"
background_file = "240418-screening_C18_acid_acetic-C-72.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "methyl-benzenoicacid"
background_file = "240417-screening_C18_a.acetic-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "bis-diphenylphosphino-ethane"
background_file = "240417-screening_C18_a.acetic-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzenesulfonylchloride"
background_file = "240417-screening_C18_a.acetic-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "n-phenylbisdifluoromethanesulfonimide"
background_file = "240417-screening_C18_a.acetic-08.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "bis-oxo-oxazolidinylphonphonicacid"
background_file = "240417-screening_C18_a.acetic-08.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "cinnamylalcohol"
background_file = "240417-screening_C18_a.acetic-08.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "eugenol"
background_file = "240417-screening_C18_a.acetic-15.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "ethylbenzoylacetate"
background_file = "240417-screening_C18_a.acetic-15.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzoin"
background_file = "240417-screening_C18_a.acetic-15.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "aminobutyricacid"
background_file = "240417-screening_C18_a.acetic-22.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethoxyaniline"
background_file = "240417-screening_C18_a.acetic-22.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethoxyethane"
background_file = "240417-screening_C18_a.acetic-22.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethoxyhydroxyphenylacetic"
background_file = "240417-screening_C18_a.acetic-29.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dihydroxybiphenyl"
background_file = "240417-screening_C18_a.acetic-29.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "ditertbutylmethylphenol"
background_file = "240417-screening_C18_a.acetic-29.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethoxy-hydroxybenzaldehyde"
background_file = "240417-screening_C18_a.acetic-36.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "3-5-dimethoxyphenol"
background_file = "240417-screening_C18_a.acetic-36.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethylglutaricanhydride"
background_file = "240417-screening_C18_a.acetic-36.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "diacetoxybenzaldehyde"
background_file = "240417-screening_C18_a.acetic-43.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethoxyhydroacetophenone"
background_file = "240417-screening_C18_a.acetic-43.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "diethylamino-dimethylfluoran"
background_file = "240417-screening_C18_a.acetic-43.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "palladiumacetate"
background_file = "240412_screening_HILIC_acetic_acid-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "imidazolcarboxylicacid"
background_file = "240412_screening_HILIC_acetic_acid-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "aminodihydropyrimidine"
background_file = "240412_screening_HILIC_acetic_acid-02.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "pyridylpiperazine"
background_file = "240412_screening_HILIC_acetic_acid-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "isopropylacetoacetate"
background_file = "240412_screening_HILIC_acetic_acid-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "hexanemethyldisilane"
background_file = "240412_screening_HILIC_acetic_acid-09.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "aniline-sulfonicacid"
background_file = "240412_screening_HILIC_acetic_acid-16.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "acesulfane"
background_file = "240412_screening_HILIC_acetic_acid-16.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "hydroxynitropyridine"
background_file = "240412_screening_HILIC_acetic_acid-21.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "heptanol"
background_file = "240412_screening_HILIC_acetic_acid-21.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzaldehyde"
background_file = "240412_screening_HILIC_acetic_acid-21.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "menthol"
background_file = "240412_screening_HILIC_acetic_acid-21.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "methylbutanol"
background_file = "240412_screening_HILIC_acetic_acid-30.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "methyltriphenylphosphoniumiodide"
background_file = "240412_screening_HILIC_acetic_acid-30.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "hydroxymethyltriphenylphosphoniumchloride"
background_file = "240412_screening_HILIC_acetic_acid-35.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "aspartame"
background_file = "240412_screening_HILIC_acetic_acid-35.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "tetrazacyclotetradecane"
background_file = "240412_screening_HILIC_acetic_acid-35.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "tetracyanoquinodimethane"
background_file = "240412-screening_C18_acetic_acide_b-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzoferrocene"
background_file = "240412-screening_C18_acetic_acide_b-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "diphenyldibromide"
background_file = "240412-screening_C18_acetic_acide_b-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "labetalol"
background_file = "240412-screening_C18_acetic_acide_b-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "chrysin"
background_file = "240412-screening_C18_acetic_acide_b-10.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "diamylphtalate"
background_file = "240412-screening_C18_acetic_acide_b-10.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "cinnamaldehyde"
background_file = "240412-screening_C18_acetic_acide_b-10.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "diethylethylmaloate"
background_file = "240412-screening_C18_acetic_acide_b-17.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "tetrazacyclotetrsdecane"
background_file = "240412-screening_C18_acetic_acide_b-17.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "amino-triazole"
background_file = "240412-screening_C18_acetic_acide_d-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "dimethylethylpyrrole"
background_file = "240412-screening_C18_acetic_acide_d-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "mebendazole"
background_file = "240412-screening_C18_acetic_acide_d-01.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzophenone"
background_file = "240412-screening_C18_acetic_acide_d-08.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "clofarabine"
background_file = "240412-screening_C18_acetic_acide_d-08.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "palladiumpicinnamylcloridedimer"
background_file = "240412-screening_C18_acetic_acide_d-08.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "cisstilbeneoxide"
background_file = "240412-screening_C18_acetic_acide_d-15.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "xantphos"
background_file = "240412-screening_C18_acetic_acide_d-15.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "bisdiphenylphpsphinobinaphthyl"
background_file = "240412-screening_C18_acetic_acide_d-15.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "1-4-pyridylpiperazine"
background_file = "240412-Screening-C18-acetic_acid_b-03.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "methylantranilate"
background_file = "240412-Screening-C18-acetic_acid_b-03.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "isopropylacetoacetate2"
background_file = "240412-Screening-C18-acetic_acid_b-03.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "hexanemethyldisilane2"
background_file = "240412-Screening-C18-acetic_acid_b-03.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzo-h-quionoline"
background_file = "240412-Screening-C18-acetic_acid_b-03.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "1-heptanol"
background_file = "240412-Screening-C18-acetic_acid_b-14.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""
molecule_folder = "benzaldehyde2"
background_file = "240412-Screening-C18-acetic_acid_b-14.JSON"

print(molecule_folder)
settings["peak_folder_time"] = timestamp_peakfolder + molecule_folder
background_filepath = os.path.join(directory_project, "Data_examples", "testfiles",
                                   background_file)
data_processing.create_new_background_spectra(background_filepath, background_method, settings)
runs.analyse_multiple_runs(molecule_folder, method_name, background_method, settings)
"""
---------------------------------------------------------------------------------------
"""




print(time.time())
print(time.time()-starttime)
