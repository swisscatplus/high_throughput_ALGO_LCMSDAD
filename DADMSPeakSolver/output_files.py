import os
from netCDF4 import Dataset
import pandas as pd
from . import data_processing as dpr
from . import data_processing_dad as dpr_dad
from . import initialize as init

def dtb_molecule_name(dtb_path):
    dtb_file = Dataset(dtb_path, "r")

    molecule_name = dtb_file.molecule_name

    dtb_file.close()
    return molecule_name

def dtb_molecule_methods(dtb_path):
    dtb_file = Dataset(dtb_path, "r")

    methods_list = []
    for name, method_name in dtb_file.groups["Methods"].groups.items():
        methods_list.append(name)

    dtb_file.close()
    return methods_list

def dtb_molecule_spectra(file_path, processed):
    dtb_file = Dataset(file_path, "r")

    if processed:
        data_ms = pd.DataFrame({
            "m/z": dtb_file.groups["MS Data"].groups["Processed MS Data"].variables["mass"][:].compressed(),
            "Intensity": dtb_file.groups["MS Data"].groups["Processed MS Data"].variables["Intensity"][:].compressed(),
        })
        data_dad = pd.DataFrame({
            "Wavelength": dtb_file.groups["DAD Data"].groups["Processed DAD Data"].variables["Wavelength"] \
                [:].compressed(),
            "Intensity": dtb_file.groups["DAD Data"].groups["Processed DAD Data"].variables["Intensity"] \
                [:].compressed(),
            # "FT Wavelength": signal_file.groups["DAD Data"].groups["Processed DAD Data"].variables["FT Wavelength"]\
            # [:].compressed(),  # uncomment once FT is implemented
        })
    else:
        data_ms = pd.DataFrame({
            "m/z": dtb_file.groups["MS Data"].groups["Raw MS Data"].variables["mass"][:].compressed(),
            "Intensity": dtb_file.groups["MS Data"].groups["Raw MS Data"].variables["Intensity"][:].compressed(),
        })
        data_dad = pd.DataFrame({
            "Wavelength": dtb_file.groups["DAD Data"].groups["Raw DAD Data"].variables["Wavelength"] \
                [:].compressed(),
            "Intensity": dtb_file.groups["DAD Data"].groups["Raw DAD Data"].variables["Intensity"] \
                [:].compressed(),
        })

    ms_spectrum = dpr.MS_Spectra(data_ms, init.import_spectra_info())
    dad_spectrum = dpr_dad.DAD_Spectra(data_dad, init.import_spectra_info())
    ms_spectrum.info["Normalized"] = processed
    dad_spectrum.info["Normalized"] = processed
    ms_spectrum.info["Molecule"] = dtb_file.molecule_name
    dad_spectrum.info["Molecule"] = dtb_file.molecule_name

    dtb_file.close()
    return ms_spectrum, dad_spectrum

def peak_spectra(peak_path, processed):
    peak_file = Dataset(peak_path, "r")

    if processed:
        data_ms = pd.DataFrame({
            "m/z": peak_file.groups["MS Data"].groups["Processed MS Data"].variables["mass"][:].compressed(),
            "Intensity": peak_file.groups["MS Data"].groups["Processed MS Data"].variables["Intensity"][:].compressed(),
        })
        data_dad = pd.DataFrame({
            "Wavelength": peak_file.groups["DAD Data"].groups["Processed DAD Data"].variables["Wavelength"] \
                [:].compressed(),
            "Intensity": peak_file.groups["DAD Data"].groups["Processed DAD Data"].variables["Intensity"] \
                [:].compressed(),
            # "FT Wavelength": signal_file.groups["DAD Data"].groups["Processed DAD Data"].variables["FT Wavelength"]\
            # [:].compressed(),  # uncomment once FT is implemented
        })
    else:
        data_ms = pd.DataFrame({
            "m/z": peak_file.groups["MS Data"].groups["Raw MS Data"].variables["mass"][:].compressed(),
            "Intensity": peak_file.groups["MS Data"].groups["Raw MS Data"].variables["Intensity"][:].compressed(),
        })
        data_dad = pd.DataFrame({
            "Wavelength": peak_file.groups["DAD Data"].groups["Raw DAD Data"].variables["Wavelength"] \
                [:].compressed(),
            "Intensity": peak_file.groups["DAD Data"].groups["Raw DAD Data"].variables["Intensity"] \
                [:].compressed(),
        })
    ms_spectrum = dpr.MS_Spectra(data_ms, init.import_spectra_info())
    dad_spectrum = dpr_dad.DAD_Spectra(data_dad, init.import_spectra_info())
    ms_spectrum.info["Normalized"] = processed
    dad_spectrum.info["Normalized"] = processed

    peak_file.close()
    return ms_spectrum, dad_spectrum

def peak_molecule_name(peak_path):
    peak_file = Dataset(peak_path, "r")

    if not peak_file.dtb_hits == "":
        molecule_name = peak_file.dtb_hits
    # implement equal rt later
    elif not (peak_file.groups["dtb hits only MS"].dtb_hits == "" and
              peak_file.groups["dtb hits only DAD"].dtb_hits == ""):
        molecule_name = "Unknown Molecule " \
                        + "Hits only MS: " + peak_file.groups["dtb hits only MS"].dtb_hits \
                        + " Hits only DAD: " + peak_file.groups["dtb hits only DAD"].dtb_hits
    else:
        molecule_name = "Unknown Molecule"

    peak_file.close()

    return molecule_name

def peak_inchi(peak_path):
    peak_file = Dataset(peak_path, "r")
    if peak_file.several_hits == "Yes" or peak_file.dtb_hits_inchi == "":
        return None
    else:
        inchi = peak_file.dtb_hits_inchi

    peak_file.close()
    return inchi

def peak_retention_time(peak_path):
    peak_file = Dataset(peak_path, "r")
    retention_time = peak_file.time
    peak_file.close()
    return retention_time

def run_chromatograms(run_path, settings):
    full_analysis = init.import_run_json(run_path, method=settings["method_name"])

    if full_analysis.info["plus_minus_acq"]:
        if settings["ion_detection_mode"] == "positive":
            ms_chrom = pd.DataFrame({
                "time": full_analysis.ms_data3d["time"][::2],
                "total intensity": full_analysis.ms_data3d["total intensity"][::2]
            })
        elif settings["ion_detection_mode"] == "negative":
            ms_chrom = pd.DataFrame({
                "time": full_analysis.ms_data3d["time"][1::2],
                "total intensity": full_analysis.ms_data3d["total intensity"][1::2]
            })
        else:
            raise NameError("Please enter a valid ion detection mode.")
    else:
        ms_chrom = pd.DataFrame({
            "time": full_analysis.ms_data3d["time"],
            "total intensity": full_analysis.ms_data3d["total intensity"]
        })
    dad_chrom = pd.DataFrame({
        "time": full_analysis.dad_data3d["time"],
        "total intensity": full_analysis.dad_data3d["total intensity"]
    })
    return ms_chrom, dad_chrom

def associate_peaks_runs(run_folder, peak_folder, settings):
    peaks_directory = os.path.join(settings["directory_project"], "Data_examples", "Peak_files", peak_folder)
    runs_directory = os.path.join(settings["directory_project"], "Data_examples", "testfiles", run_folder)
    run_list = os.listdir(runs_directory)
    peak_list = os.listdir(peaks_directory)


    run_nr = 1.
    all_peaks = pd.DataFrame()
    for run_name in run_list:
        for peak_name in peak_list:
            peak_path = os.path.join(peaks_directory, peak_name)
            all_runs = all_runs_details(peak_path)
            for peak_run_nr in all_runs["Run Nr."]:
                if peak_run_nr == run_nr:
                    all_peaks = None

        run_nr += 1
    return

def all_peaks_details(run_nr, peak_folder, settings):
    peaks_directory = os.path.join(settings["directory_project"], "Data_examples", "Peak_files", peak_folder)
    peak_list = os.listdir(peaks_directory)

    combined_peak_numbers = []
    combined_ret_times = []
    combined_rel_pur = []
    combined_integrals = []
    combined_pure = []
    for peak_name in peak_list:
        peak_path = os.path.join(peaks_directory, peak_name)
        all_runs = all_runs_details(peak_path)
        for peak_run_nr in all_runs["Run Nr."]:
            index = 0
            if peak_run_nr == float(run_nr):
                peak_nr = determine_peak_nr(peak_path)
                combined_peak_numbers.append(peak_nr)
                combined_ret_times.append(all_runs["Retention Time"][index])
                combined_rel_pur.append(all_runs["Relative Purity"][index])
                combined_integrals.append(all_runs["Integral"][index])
                combined_pure.append(all_runs["Pure"][index])
            index += 1

    if combined_peak_numbers == []:
        combined_peak_numbers = [None]
        combined_ret_times = [None]
        combined_rel_pur = [None]
        combined_integrals = [None]
        combined_pure = [None]

    all_peaks = pd.DataFrame({
        "Peak Nr.": combined_peak_numbers,
        "Retention Time": combined_ret_times,
        "Pure": combined_pure,
        "Relative Purity": combined_rel_pur,
        "Integral": combined_integrals,
    })
    return all_peaks

def determine_peak_nr(peak_path):
    peak_file = Dataset(peak_path, "r")
    peak_nr = float(peak_file.peak_nr)
    peak_file.close()
    return peak_nr

def all_runs_details(peak_path):
    peak_file = Dataset(peak_path, "r")

    pure_boolean = []
    for pure_value in peak_file.variables["All pure"][:].compressed():
        if pure_value == 1:
            pure_boolean.append(True)
        elif pure_value == 0:
            pure_boolean.append(False)
        else:
            pure_boolean.append("Unclear Value")

    all_runs = pd.DataFrame({
        "Run Nr.": peak_file.variables["All runs"][:].compressed(),
        "Retention Time": peak_file.variables["All times"][:].compressed(),
        "Pure": pure_boolean,
        "Relative Purity": peak_file.variables["All relative purity"][:].compressed(),
        "Integral": peak_file.variables["All integrals"][:].compressed(),
    })

    peak_file.close()
    return all_runs