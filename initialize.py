from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import sys
import json

def load_analysis_run_netCDF(analysis_path):
    """
    Load the cdf file
    For the moment only MS Data
    Extract the data we need
    :param analysis_path: path of run to analyse
    :return: Return dataframe with nexessary data and annotations
    """
    #np.set_printoptions(threshold=sys.maxsize)  # To print complete arrays
    """ All informations on unparsed .netCDF file
    Scan acquisition time gives time when scan was done in seconds.
    Scan duration & inter scan time & resolution are default and constant.
    So are a d sampling rate & coaddition factor.
    
    Actual scan number: Real index from 0 to 1593.
    Total intensity: Sum of all intensities, indexed on actual scan number.
    
    Scan index: indexed itself on actual scan number. Could be the measured elements per scan in "intensity_values?
    -> scan index gives first value of mass/intensities which belongs to the corresponding scan!
    
    Mass range min is always 0. Mass range max changes, always > 998 m/z.
    Time range min&max are default and constant.
    
    Point count could be the number of m/z values measured per spectra.
    Flag count is always 0, not sure what it is.
    
    Mass & intensity values: Dimension of point number, give corresponding value...
    Time value: Set to default, constant. Not sure what it should contain initially.
    """

    data_cdf = Dataset(analysis_path)  # Imports the readable .netCDF file.

    """ What data to extract?
    Each index contains the time, the total intensity and a 2D spectra array.
    The 2D spectra array contains the m/z values and their corresponding intensities.
    -> They can be obtained via the scan index of the corresponding actual scan number.
    
    Ideally also parse the name of the LC Method, but not contained.. Need to parse it externally later
    """
    # dtype = [("time", "float32"), ("total_intensity", "float32"), ("ms_spectra", "float32")]
    chromatogram_ms = np.empty((len(data_cdf.variables["actual_scan_number"]), 2)) # Initialize array for data.
    ms_data = np.empty((len(data_cdf.variables["actual_scan_number"]), 2), dtype=object)

    for i in data_cdf.variables["actual_scan_number"]:  # Iterates over all measured scans.
        chromatogram_ms[i,0] = data_cdf.variables["scan_acquisition_time"][i]
        chromatogram_ms[i,1] = data_cdf.variables["total_intensity"][i]
        begin_pn_index = data_cdf.variables["scan_index"][i]

        if (i+1) >= data_cdf.variables["actual_scan_number"][-1]: # For last scan, since end pn index isn't defined
            ms_data[i,0] = [data_cdf.variables["mass_values"][begin_pn_index:]]
            ms_data[i,1] = [data_cdf.variables["intensity_values"][begin_pn_index:]]
        else:
            end_pn_index = data_cdf.variables["scan_index"][i+1]
            ms_data[i, 0] = [data_cdf.variables["mass_values"][begin_pn_index:end_pn_index]]
            ms_data[i, 1] = [data_cdf.variables["intensity_values"][begin_pn_index:end_pn_index]]
        ms_spectrum = np.array(ms_data)

    full_chromatogram_ms = {
        "time": chromatogram_ms[:,0],
        "total intensity": chromatogram_ms[:,1],
        "MS spectra": ms_spectrum, # 0 to access mass values and 1 to access intensity values
        }
    return full_chromatogram_ms # returns chromatogram in array structure with all information needed

def import_spectra_info(database = False, merged_spectra = False, time =0, run_nr = "1", sgl_nr = "1"):
    """
    Imports all additional data for the spectra.
    Predefined for test runs, will be parsed directly once extraction is automated
    :return:
    """
    spectra_info = {
        "Molecule": "Unknown",
        "Normalized": False,
        "Analysis_peak": True,
        "Weight Fct": "",
        "From database": database,
        "Merged Spectra": merged_spectra,
        "Background correction": False,
        "LC Method": "Testmethod",
        "Database matches": [],
        "Time": time,
        "Run Nr": run_nr,
        "Signal Nr": sgl_nr,
        "Processed": False,
        }
    return spectra_info

def import_full_chr_info(background = False, run_nr = "1", method = "Testmethod"):
    """
    Imports all additional data for the spectra.
    Predefined for test runs, will be parsed directly once extraction is automated
    :return:
    """
    chr_info = {
        "Normalized": False,
        "LC Method": method,
        "Background": background,
        "Run Nr": run_nr,
        "Name": "",
    }
    return chr_info


def import_dad_spectrum(dad_path):
    """
    Fct to parse preliminary file with dad spectrum as obtained from MestreNova Plugin.
    Parse data as pd dataframe.
    :param dad_path:
    :return:
    """
    uv_spectrum = pd.read_csv(dad_path, sep = "\t", names = ["Wavelength", "Intensity"])
    print(uv_spectrum)
    plt.plot(uv_spectrum["Wavelength"], uv_spectrum["Intensity"])
    plt.show()
    return

def import_run_json(path, run_nr = str(1), method = None):
    """
    Imports the data from a .json file which contains 3D dad and ms data.
    Ideally parse the method type directly, but currently the information is not included in the file.
    :param path:
    :return:
    """
    """
    json file structure: (dictionary type in python)
    liquid chromatography aggregate document - dict type
        device system document  # a bunch of info on the device, no need to be parsed
        liquid chromatography document - list type (length 1)  # what we care about
            [0] - dict type
                analyst - str type  # just gives name of analyst
                measurement aggregate document
                    measurement document - list type (length 5)
                        0 - dict type (length 8) # DAD1B chromatogram
                            chromatogram data cube - dict type (length 4)
                                data is dict containing measures (Intensity) and dimensions (time)
                                Other elements contain additional informatio
                                Given for lambda = 254.4 nm
                        1 - dict type (length 8)  # MS+ chromatogram
                            chromatogram data cube
                                as above, contains data with MS intensity and dimension (time
                            processed data document
                                contains list of identified peaks etc
                            sample document
                                information on the sample (including name!)
                            injection document
                                information on method, could be possible to extract method type!
                        2 - dict type (length 8)  # MS- chromatogram
                            same as above but yields data in negative mass mode
                        3 - dict type (length 8)  # DAD1 Intensity 3D
                            measurement identifier
                            chromatography column document  # empty list (?)
                            device control aggregate document  # dict of device infos...
                            sample document  # sample identifier and name
                            injection document
                            detection type  # just says single channel
                            chromatogram data cube - dict type  # empty chromatogram
                            three-dimensional ultraviolet spectrum data cube - dict type (length 4)  # 3D DAD data
                                label  # DAD spectrum
                                cube-structure  # explains how 3D cube is build:
                                    dimensions: list
                                        retention time in secounds (double)
                                        wavelength in nm (double)
                                    measures: list
                                        absorbance in mAU (double)
                                data - dict type
                                    measures - list type (length 1)
                                        0 - list (length 3000)  # contains absorbance values for each point in time (3000)
                                            0 - list (length 106)  # contains absorbance values for each wavelength (106) for time[0]
                                            1 - same but for time[1]
                                            ...
                                    dimensions - list type (length 2)
                                        0 - list (length 3000)  # contains time values in secounds
                                            0...3000 contains ordered value from 0.2 s to 600.0 s in 0.2 s
                                        1 - list (length 106)  # contains the wavelength values in nm
                                            ordered from 190.0 to 400.0 in 2.0
                                identifier  # DAD1I
                        4 - dict type (length 8)  # 3D MS Data
                            measurement identifier # SQ1A and time stamp
                            rest is just as above
                            three-dimensional mass spectrum data cube - dict type (length 4)  # 3D MS Data
                                label  # SQ1A +/-ESI  Scan Frag=100.0V 
                                cube-structure  # as above
                                    dimensions for time and m/z values
                                    measures for intensity values
                                data - list type (length 3188)
                                    0 - dict type (length 3)
                                        time  # numeric time value as float [secounds]
                                        measures  # list (length 1)
                                            0 - list (length 784)
                                                intensity values as floats
                                        dimensions  # list (length 1)
                                            0 - list (length 784)  # needs to be same length as measures [0]
                                                m/z values as floats (only those actually measured)
                                    1 - same as above
                                        but measures[0] has length 753! only recorded m/z are saved...
                                identifier  # SQ1A 20230907 145900          
                @index - int type # returns 1
                
    The only information that is missing is the aquisition method...
    """
    with open(path, "r") as json_file:
        data = json_file.read()
        data_deserialized = json.loads(data)
    measurement_document = data_deserialized["liquid chromatography aggregate document"]\
        ["liquid chromatography document"][0]["measurement aggregate document"]["measurement document"]
    dad_data3d = measurement_document[3]["three-dimensional ultraviolet spectrum data cube"]["data"]  # path for 3d dad
    ms_data3d = measurement_document[4]["three-dimensional mass spectrum data cube"]["data"]  # path for 3d ms data

    """for key in measurement_document[1]:
        print(key)
    print(measurement_document[1]["measurement identifier"])"""
    #print(data_deserialized["liquid chromatography aggregate document"]["device system document"])

    # first parse the 3D ms data
    chromatogram_ms = np.empty((len(ms_data3d), 2)) # Initialize array for data.
    ms_data = np.empty((len(ms_data3d), 2), dtype=object)

    for i in range(len(ms_data3d)):
        chromatogram_ms[i, 0] = ms_data3d[i]["time"]  # save time value
        chromatogram_ms[i, 1] = sum(ms_data3d[i]["measures"][0])  # save total intensity per time
        ms_mass_values = np.empty(len(ms_data3d[i]["measures"][0]))
        ms_intensity_values = np.empty(len(ms_data3d[i]["measures"][0]))
        for index in range(len(ms_data3d[i]["measures"][0])):
            ms_mass_values[index] = ms_data3d[i]["dimensions"][0][index]  # save m/z values
            ms_intensity_values[index] = ms_data3d[i]["measures"][0][index]  # save intensity values
        ms_data[i, 0] = ms_mass_values  # saving data points in time dependency
        ms_data[i, 1] = ms_intensity_values
    # load all data into pd dataframe for processing
    ms_spectrum = np.array(ms_data)
    full_chromatogram_ms = {
        "time": chromatogram_ms[:, 0],
        "total intensity": chromatogram_ms[:, 1],
        "MS spectra": ms_spectrum,  # 0 to access m/z values and 1 to access intensity values
    }


    # parse 3D UV data
    # PRELIMINARY, maybe change of mocca needs it different.. just for testing of all data is here
    chromatogram_dad = np.empty((len(dad_data3d["dimensions"][0]), 2))
    dad_data = np.empty((len(dad_data3d["dimensions"][0]), 2), dtype = object)
    dad_wavelengths = np.array(dad_data3d["dimensions"][1]) # save wavelengths (no time dependency)

    for i in range(len(dad_data3d["dimensions"][0])):
        chromatogram_dad[i, 0] = dad_data3d["dimensions"][0][i]  # save time value
        chromatogram_dad[i, 1] = sum(dad_data3d["measures"][0][i])  # save total intensity per time
        dad_intensity_values = np.empty(len(dad_data3d["measures"][0][i]))
        for index in range(len(dad_data3d["measures"][0][i])):
            dad_intensity_values[index] = dad_data3d["measures"][0][i][index]  # save intensity values
        dad_data[i, 0] = dad_wavelengths # saving data points in time dependency
        dad_data[i, 1] = dad_intensity_values
    # load all data into pd dataframe for processing
    dad_spectrum = np.array(dad_data)
    full_chromatogram_dad = {
        "time": chromatogram_dad[:, 0],
        "total intensity": chromatogram_dad[:, 1],
        "DAD spectra": dad_spectrum,  # 0 to access wavelengths and 1 to access intensity values
    }

    if method == None:
        raise ValueError("Please enter a name for the method used.")
    info = import_full_chr_info(run_nr=run_nr, method = method)
    info["Name"] = measurement_document[0]["sample document"]["written name"]

    return full_analysis(full_chromatogram_ms, full_chromatogram_dad, info)

class full_analysis:
    """
    Object to contain the complete data and information of an analysis run. Parsed by any document type.
    """
    def __init__(self, ms_data3d, dad_data3d, info):
        self.ms_data3d = ms_data3d
        self.dad_data3d = dad_data3d
        self.info = info
