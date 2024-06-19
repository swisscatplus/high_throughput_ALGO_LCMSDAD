import os
from datetime import datetime
import pandas as pd
from netCDF4 import Dataset
import numpy as np

import create_file
import data_processing
import data_processing_dad as dpr_dad
import initialize as init
import ms_spectra_comparison as ms_comp
import dad_spectra_comparison as dad_comp

"""
Handling of the signal and peak files. Comparison with each other and rewriting of files.
"""

def signal_comparison(settings):
    """
    Iterate over all signal files and compare with each other. If they are equal, note the run they occur in and delete
    the less pure spectrum. Afterwards
    :param directory_project:
    :return:
    """
    directory_project = settings["directory_project"]
    signal_files_path = os.path.join(directory_project, "Data_examples", "Signal_files")
    signal_list = os.listdir(signal_files_path)

    for signal in signal_list:
        first_signal_path = os.path.join(signal_files_path, signal)
        for signal2 in signal_list:
            if not signal == signal2: # compare
                second_signal_path = os.path.join(signal_files_path, signal2)
                try:
                    first_signal = load_signal_file(first_signal_path)
                    sec_signal = load_signal_file(second_signal_path)
                except FileNotFoundError:  # Comes up bcs equal signal files are being deleted
                    continue

                old_threshold_ms = settings["threshold_ms_dot_product_sim"]
                old_threshold_dad = settings["threshold_dad_derivative"]
                settings["threshold_ms_dot_product_sim"] = .8
                settings["threshold_dad_derivative"] = .8
                comparison_ms = ms_comp.comparison_ms(first_signal.ms_spectrum, sec_signal.ms_spectrum, settings)
                comparison_dad = dad_comp.comparison_dad(first_signal.dad_spectrum, sec_signal.dad_spectrum, settings)
                comparison_rf = dad_comp.compare_retention_time(first_signal.info["Retention Time"],
                                                                sec_signal.info["Retention Time"], settings)
                settings["threshold_ms_dot_product_sim"] = old_threshold_ms
                settings["threshold_dad_derivative"] = old_threshold_dad
                if comparison_ms and comparison_dad and comparison_rf:
                    chosen_signal = choose_better_signal(first_signal, sec_signal)
                    combined_run_nr = np.hstack((first_signal.all_runs["All runs"], sec_signal.all_runs["All runs"]))
                    combined_ret_times = np.hstack((first_signal.all_runs["Retention times"], sec_signal.all_runs["Retention times"]))
                    combined_rel_pur = np.hstack((first_signal.all_runs["Relative Purities"], sec_signal.all_runs["Relative Purities"]))
                    combined_integrals = np.hstack((first_signal.all_runs["Integrals"], sec_signal.all_runs["Integrals"]))
                    combined_pure = np.hstack((first_signal.all_runs["Pure"], sec_signal.all_runs["Pure"]))

                    if chosen_signal == 1:
                        signal_file = Dataset(first_signal_path, "r+")
                        signal_file.variables["All runs"][:] = combined_run_nr
                        signal_file.variables["All times"][:] = combined_ret_times
                        signal_file.variables["All relative purity"][:] = combined_rel_pur
                        signal_file.variables["All integrals"][:] = combined_integrals
                        signal_file.variables["All pure"][:] = combined_pure
                        signal_file.close()
                        os.remove(second_signal_path)
                    elif chosen_signal == 2:
                        signal_file = Dataset(second_signal_path, "r+")
                        signal_file.variables["All runs"][:] = combined_run_nr
                        signal_file.variables["All times"][:] = combined_ret_times
                        signal_file.variables["All relative purity"][:] = combined_rel_pur
                        signal_file.variables["All integrals"][:] = combined_integrals
                        signal_file.variables["All pure"][:] = combined_pure
                        signal_file.close()
                        os.remove(first_signal_path)
                    signal_list = os.listdir(signal_files_path)

    signal_list_updated = os.listdir(signal_files_path)

    peak_folder_path = os.path.join(directory_project, "Data_examples", "Peak_files", settings["peak_folder_time"])
    os.makedirs(peak_folder_path, exist_ok=True)

    peak_nr = 1
    for signal in signal_list_updated:  # Creates a peak file for every signal file that still exists
        signal_file_path = os.path.join(signal_files_path, signal)
        create_file.create_peak_file(signal_file_path, peak_folder_path, peak_nr)  # Not yet tested!
        # The number of all runs of that signal should be saved within the peak file
        peak_nr += 1

    return

def choose_better_signal(first_signal, sec_signal):
    """
    Select one of two signals for the creation of a peak file based on quantity and purity.
    :param first_signal:
    :param sec_signal:
    :param settings:
    :return:
    """
    compare_purity = first_signal.info["Relative Purity"] / sec_signal.info["Relative Purity"]
    compare_integral = first_signal.info["Integral"] / sec_signal.info["Integral"]
    if first_signal.info["Pure"] and sec_signal.info["Pure"]:
        if compare_integral >= 1.:
            return 1
        else:
            return 2
    elif first_signal.info["Pure"] and not sec_signal.info["Pure"]:
        return 1
    elif sec_signal.info["Pure"] and not first_signal.info["Pure"]:
        return 2
    elif (compare_purity >= 1.) and (compare_integral >= 1.):
        return 1
    elif (compare_purity < 1.) and (compare_integral < 1.):
        return 2
    elif (compare_purity >= 1.) and (compare_integral < 1.):
        if compare_purity >= (1/compare_integral):
            return 1
        else:
            return 2
    elif (compare_purity < 1.) and (compare_integral >= 1.):
        if (1/compare_purity) >= compare_integral:
            return 2
        else:
            return 1
    else:
        print("Warning: This signal file selection criteria shouldn't occur.")
        return 1


def load_signal_file(path, processed = True):
    """
    Load previously created signal for comparisons.
    :param path:
    :return: Object containing all processed spectra and their info.
    """
    signal_file = Dataset(path, "r")
    if processed:
        data_ms = pd.DataFrame({
            "m/z": signal_file.groups["MS Data"].groups["Processed MS Data"].variables["mass"][:].compressed(),
            "Intensity": signal_file.groups["MS Data"].groups["Processed MS Data"].variables["Intensity"][:].compressed(),
        })
        data_dad = pd.DataFrame({
            "Wavelength": signal_file.groups["DAD Data"].groups["Processed DAD Data"].variables["Wavelength"]\
                [:].compressed(),
            "Intensity": signal_file.groups["DAD Data"].groups["Processed DAD Data"].variables["Intensity"]\
                [:].compressed(),
            # "FT Wavelength": signal_file.groups["DAD Data"].groups["Processed DAD Data"].variables["FT Wavelength"]\
            # [:].compressed(),  # uncomment once FT is implemented
        })
    else:
        data_ms = pd.DataFrame({
            "m/z": signal_file.groups["MS Data"].groups["Raw MS Data"].variables["mass"][:].compressed(),
            "Intensity": signal_file.groups["MS Data"].groups["Raw MS Data"].variables["Intensity"][:].compressed(),
        })
        data_dad = pd.DataFrame({
            "Wavelength": signal_file.groups["DAD Data"].groups["Raw DAD Data"].variables["Wavelength"] \
                [:].compressed(),
            "Intensity": signal_file.groups["DAD Data"].groups["Raw DAD Data"].variables["Intensity"] \
                [:].compressed(),
        })

    info = {
        "Processed": processed,
        "Database": False,
        "Method": signal_file.Method,
        "Retention Time": signal_file.time,
        "Run Nr": signal_file.run_nr,
        "Relative Purity": signal_file.relative_purity,
        "Integral": signal_file.integral,
    }

    all_runs = {
        "All runs": signal_file.variables["All runs"][:].compressed(),
        "Retention times": signal_file.variables["All times"][:].compressed(),
        "Relative Purities": signal_file.variables["All relative purity"][:].compressed(),
        "Integrals": signal_file.variables["All integrals"][:].compressed(),
        "Pure": signal_file.variables["All pure"][:].compressed(),
    }

    if signal_file.pure == "True":
        info["Pure"] = True
    else:
        info["Pure"] = False

    print(all_runs["All runs"])

    ms_spectrum = data_processing.MS_Spectra(data_ms, init.import_spectra_info())
    dad_spectrum = dpr_dad.DAD_Spectra(data_dad, init.import_spectra_info())
    ms_spectrum.info["Normalized"] = processed
    dad_spectrum.info["Normalized"] = processed
    signal = Signal(ms_spectrum, dad_spectrum, info, all_runs)

    signal_file.close()
    return signal


class Signal:
    """
      Class for a Signal, containing the information on the method, and the processed DAD and MS spectra and
      the Rf value.
      Used when dealing with either peaks or signals.
      """

    def __init__(self, ms_spectrum, dad_spectrum, info, all_runs = None):
        self.ms_spectrum = ms_spectrum
        self.dad_spectrum = dad_spectrum
        self.info = info
        self.all_runs = all_runs

        # Not necessary to correct for decimals since this is in the initalization of ms_spectra..