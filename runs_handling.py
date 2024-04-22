import os
import combine_ms_dad as dadms
import initialize as init
import signal_peak_handling as sp_handling
import databank_handling as dtb

def analyse_single_run(run_name, method_name, background_name, settings):
    """
    From run name to database comparison and creation of output files.
    :param run_name: Needs to be entered as a string.
    :param method_name: Needs to be entered as a string.
    :param settings:
    :return:
    """
    directory_project = settings["directory_project"]
    json_path = os.path.join(directory_project, "Data_examples", "testfiles",
                             run_name)
    delete_old_sgn_files(settings)

    full_analysis = init.import_run_json(json_path, method=background_name)
    dadms.assign_peaks(full_analysis, method_name, settings)

    sp_handling.signal_comparison(settings)

    peak_directory = os.path.join(directory_project, "Data_examples", "Peak_files", settings["peak_folder_time"])
    dtb.compare_all_peaks_dtb(peak_directory, settings)

    del_default_method(background_name, settings)
    return

def analyse_multiple_runs(run_folder_name, method_name, background_name, settings):
    """
    From run folder to dtb comparison and creation of output files.
    :param run_folder_name: Needs to be entered as a string.
    :param method_name: Needs to be entered as a string.
    :param settings:
    :return:
    """
    directory_project = settings["directory_project"]
    run_folder_directory = os.path.join(directory_project, "Data_examples", "testfiles", run_folder_name)
    analysis_run_list = os.listdir(run_folder_directory)
    delete_old_sgn_files(settings)

    run_nr = 1
    for analysis_run_name in analysis_run_list:
        run_path = os.path.join(run_folder_directory, analysis_run_name)
        if os.path.isfile(run_path) and run_path.lower().endswith(".json"):
            full_analysis = init.import_run_json(run_path, run_nr = str(run_nr), method=background_name)
            dadms.assign_peaks(full_analysis, method_name, settings)
            run_nr += 1
        else:
            print("Warning: File " + analysis_run_name + " is not a supported analysis file type!")

    sp_handling.signal_comparison(settings)

    peak_directory = os.path.join(directory_project, "Data_examples", "Peak_files", settings["peak_folder_time"])
    dtb.compare_all_peaks_dtb(peak_directory, settings)

    del_default_method(background_name, settings)
    return

def delete_old_sgn_files(settings):
    directory_project = settings["directory_project"]
    signals_directory = os.path.join(directory_project, "Data_examples", "Signal_files")
    signal_list = os.listdir(signals_directory)

    for signal_name in signal_list:
        signal_path = os.path.join(signals_directory, signal_name)
        os.remove(signal_path)
    return

def del_default_method(background_name, settings):
    directory_project = settings["directory_project"]
    background_directory = os.path.join(directory_project, "Data_examples", "background_spectra")
    if background_name == "defaultNew":
        background_path = os.path.join(background_directory, "defaultNew.cdf")
        os.remove(background_path)
    return
