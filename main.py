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

# Go back directories until we have project folder
current_directory = os.getcwd()
directory_project = os.path.abspath(os.path.join(current_directory, os.pardir))

# Go to background folder
background_folder = os.path.join(directory_project, "Data_examples", "background_files")

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
    "ion_detection_mode": "negative",  # positive or negative
    "method_name": "AceticAcid01",
    "Number Columns": 8,  # Number for columns for heatmap visualization
}
method_name = "new_method"
background_method = "defaultNew"

"""
--------------------------------
These functions need to be called directly by their origin: (write in documentation)

runs.analyse_single_run(run_name, method_name, background_method, settings)

out.dtb_molecule_list(settings)
out.dtb_molecule_full_data("RCIJACVHOIKRAP-UHFFFAOYSA-M.cdf", settings)
out.create_analysis_report(settings, run_folder_name, report_name="_")

opt_sc.comparison_dtb_named_files(settings, False)
opt_sc.plot_optimization_dad("dad_optimization_peaks_dtb_compTrue.csv", settings)
opt_sc.plot_optimization_ms("ms_optimization_peaks_dtb_compFalse.csv", settings)

create_file.add_method_to_dtb_entry(dtb_entry_path, peak_path)
create_file.create_dtb_entry(directory_project, peak_path, molecule)
"""
def change_data_directory(new_directory_path):
    """
    Function to change the directory of the Data folder. Needs to be run before all other functions.
    :param new_directory_path: Path to contain the "Data" folder, which includes the subordinate files.
    """
    directory_project = new_directory_path
    settings["directory_project"] = directory_project
    settings["background folder"] = os.path.join(directory_project, "Data_examples", "background_files")
    return

def custom_settings(new_settings):
    """
    Function to change settings and parameters. Needs to be run before the actual analysis function.
    :param new_settings: Dictionary of all settings to be changed. The key needs to be the same as in the original settings.
    """
    for key in new_settings:
        if key in settings:
            settings[key] = new_settings[key]
            print(f"Setting changed: {key} = {settings[key]}")
    return

def new_dtb_entry(peak_folder, peak_name, molecule):
    """
    Function to create a new database entry. If the function fails due to the molecule name being unknown to PubChem,
    the create_dtb_entry function needs to be called directly.
    :param peak_folder: Name of the peak folder.
    :param peak_name: Name of the peak file to be transferred into the database.
    :param molecule: Name of the molecule (IUPAC).
    """
    peak_path = os.path.join(directory_project, "Data_examples", "Peak_files", peak_folder, peak_name)
    create_file.create_dtb_entry(directory_project, peak_path, molecule)
    return

def analyse_run_folder(run_folder, background_run, method_name, create_report = True, report_name = None):
    """
    :param run_folder: Name of the run folder to be analysed
    :param background_run: Name of the background run file to be used. Needs to be in the ASM format.
    :param method_name: Name of the LC method used.
    :param create_report: Option to create the report file. Default is true.
    :param report_name: Name of the report file. Default is the name of the run + _Report.
    """

    background_filepath = os.path.join(directory_project, "Data_examples", "testfiles", background_run)
    data_processing.create_new_background_spectra(background_filepath, background_method, settings)

    run_folder_filepath = os.path.join(directory_project, "Data_examples", "testfiles", run_folder)
    runs.analyse_multiple_runs(run_folder_filepath, method_name, background_method, settings)

    if create_report:
        if report_name is None:
            report_name = run_folder + "_Report"
        out.create_analysis_report(settings, run_folder, report_name=report_name)
    return

def save_background_method(background_run, method_name):
    """
    A function to create a background run file for a new method is used.
    :param background_run: Name of the background run used (ASM file type).
    :param method_name: Name of the new method.
    """
    background_filepath = os.path.join(directory_project, "Data_examples", "testfiles", background_run)
    data_processing.create_new_background_spectra(background_filepath,method_name, settings)
    return

