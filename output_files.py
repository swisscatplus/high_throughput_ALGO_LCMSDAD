import os
from netCDF4 import Dataset
import pandas as pd
import data_processing as dpr
import data_processing_dad as dpr_dad
import initialize as init

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
    elif not peak_file.groups["dtb hits different method"].dtb_hits == "":
        molecule_name = "Unknown Molecule " \
                        + "Hits different method: " + peak_file.groups["dtb hits different method"].dtb_hits
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