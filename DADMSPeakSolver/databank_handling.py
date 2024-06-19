import initialize as init
import os
from netCDF4 import Dataset
import pandas as pd

import data_processing
import data_processing_dad as dpr_dad
import ms_spectra_comparison as ms_comp
import signal_peak_handling as sp_handling
import dad_spectra_comparison as dad_comp

"""
Deal with file paths and spectra processing for all databank entries (for now only ms).
"""

def compare_all_peaks_dtb(peak_directory, settings):
    peak_directory_list = os.listdir(peak_directory)
    for peak_file_name in peak_directory_list:
        peak_file_path = os.path.join(peak_directory, peak_file_name)
        if os.path.isfile(peak_file_path) and peak_file_path.lower().endswith(".cdf"):
            databank_comparison_all(peak_file_path, settings)
        else:
            print("!!! File " + peak_file_name + " is not a supported peak file type!")
    return

def databank_comparison_all(analyte_path, settings):
    """
    Compares a peak of an analysis run with database spectra for MS and DAD.
    :analyte_spectra: Already fully processed.
    :return: Annotation of matches.
    """
    directory_project = settings["directory_project"]
    directory_database = os.path.join(directory_project, "Data_examples", "database")  # Database path
    database_spectra = os.listdir(directory_database)  # List of all files in the database

    analyte_signal = sp_handling.load_signal_file(analyte_path)
    analyte_spectrum_ms = analyte_signal.ms_spectrum
    analyte_spectrum_dad = analyte_signal.dad_spectrum
    analyte_method = analyte_signal.info["Method"]
    analyte_retention_time = analyte_signal.info["Retention Time"]

    for database_entry in database_spectra:
        single_entry = os.path.join(directory_database, database_entry)
        if os.path.isfile(single_entry) and single_entry.lower().endswith(".cdf"):
            dtb_spectrum = load_dtb_entry(single_entry)  # extract file
            dtb_spectrum_ms = dtb_spectrum.ms_spectrum  # Obtains MS spectrum from signal object
            dtb_spectrum_dad = dtb_spectrum.dad_spectrum
            dtb_methods = dtb_spectrum.info["Methods"]

            equal_ms_spectra = ms_comp.comparison_ms(analyte_spectrum_ms, dtb_spectrum_ms, settings)
            equal_dad_spectra = dad_comp.comparison_dad(analyte_spectrum_dad, dtb_spectrum_dad, settings)

            if analyte_method in dtb_methods.columns:
                dtb_retention_time = dtb_methods[analyte_method][0]
                equal_method = True
                equal_retention_time = dad_comp.compare_retention_time(dtb_retention_time, analyte_retention_time, settings)
            else:
                equal_method = False

            analyte_file = Dataset(analyte_path, "r+")
            if equal_method:
                if equal_ms_spectra and equal_dad_spectra and equal_retention_time:
                    # print("Database Hit!")
                    # add_molecule_to_analyte(analyte_file, dtb_spectrum)
                    current_hits = getattr(analyte_file.groups["dtb hits with retention time"], "dtb_hits")
                    new_hits = current_hits + dtb_spectrum.info["Molecule Name"] + ", "
                    analyte_file.groups["dtb hits with retention time"].dtb_hits = new_hits
            if equal_ms_spectra and equal_dad_spectra:
                # print("Both MS and DAD are equal.")
                add_molecule_to_analyte(analyte_file, dtb_spectrum)
                # add to molecule list of peak
            elif equal_ms_spectra and not equal_dad_spectra:
                # print("Only equal MS.")
                current_hits = getattr(analyte_file.groups["dtb hits only MS"], "dtb_hits")
                new_hits = current_hits + dtb_spectrum.info["Molecule Name"] + ", "
                analyte_file.groups["dtb hits only MS"].dtb_hits = new_hits
            elif equal_dad_spectra and not equal_ms_spectra:
                # print("Only equal DAD spectra.")
                current_hits = getattr(analyte_file.groups["dtb hits only DAD"], "dtb_hits")
                new_hits = current_hits + dtb_spectrum.info["Molecule Name"] + ", "
                analyte_file.groups["dtb hits only DAD"].dtb_hits = new_hits

            analyte_file.close()
        else:  # Prints error in case there is a wrong file
            print("!!! File " + single_entry + " is not a supported database file type!")
    return

def add_molecule_to_analyte(analyte_file, dtb_spectrum):
    # print("Database Hit!")
    if not getattr(analyte_file, "dtb_hits") == "":
        # print("There are several dtb hits for this peak.")
        analyte_file.several_hits = "Yes"
        current_hits = getattr(analyte_file, "dtb_hits")
        new_hits = current_hits + " and " + dtb_spectrum.info["Molecule Name"]
        analyte_file.dtb_hits = new_hits

        current_hits_inchi = getattr(analyte_file, "dtb_hits_inchi")
        new_hits_inchi = current_hits_inchi + " and " + dtb_spectrum.info["Molecule InChI"]
        analyte_file.dtb_hits_inchi = new_hits_inchi
    else:
        analyte_file.dtb_hits = dtb_spectrum.info["Molecule Name"]
        analyte_file.dtb_hits_inchi = dtb_spectrum.info["Molecule InChI"]
    return

def load_dtb_entry(dtb_entry_path, processed = True):
    """
    Load a database entry as Signal object for easy comparison. May or may not load the pre-processed data.
    :param dtb_entry_path:
    :return:
    """
    dtb_file = Dataset(dtb_entry_path, "r")
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

    methods = pd.DataFrame({})  # pd Dataframe to store all methods and retention times.

    for name, method_name in dtb_file.groups["Methods"].groups.items():
        methods[name] = method_name.variables["Retention Time"][:]

    info = {
        "Processed": processed,
        "Database": True,
        "Molecule Name": dtb_file.molecule_name,
        "Molecule InChI": dtb_file.molecule_inchi,
        "Methods": methods,
    }

    ms_spectrum = data_processing.MS_Spectra(data_ms, init.import_spectra_info())
    dad_spectrum = dpr_dad.DAD_Spectra(data_dad, init.import_spectra_info())
    ms_spectrum.info["Normalized"] = processed
    dad_spectrum.info["Normalized"] = processed
    ms_spectrum.info["Molecule"] = dtb_file.molecule_name
    dad_spectrum.info["Molecule"] = dtb_file.molecule_name

    dtb_file.close()
    dtb_entry = sp_handling.Signal(ms_spectrum, dad_spectrum, info)

    return dtb_entry