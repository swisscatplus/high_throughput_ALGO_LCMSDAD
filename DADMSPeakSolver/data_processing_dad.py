import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from . import initialize as init
import os
from netCDF4 import Dataset
from scipy.signal import savgol_filter

from .mocca.user_interaction.campaign import HplcDadCampaign
from .mocca.user_interaction.user_objects import Gradient
from .mocca.user_interaction.user_objects import Compound
from .mocca.user_interaction.user_objects import InternalStandard
from .mocca.user_interaction.user_objects import HplcInput
from .mocca.user_interaction.settings import Settings
from .mocca.dad_data.models import ParafacData


current_directory = os.getcwd()
directory_code = os.path.abspath(os.path.join(current_directory, os.pardir))
directory_project = os.path.abspath(os.path.join(directory_code, os.pardir))
testpath = os.path.join(directory_project, "Data_examples", "testfiles", "carboxylic_acid_and_amines_methode_type_poroshell_7-01.JSON")
report_path = os.path.join(directory_project, "Data_examples")

def run_mocca(full_analysis, settings):
    testgradient = None
    trial_campaign = HplcDadCampaign()
    dad_object = full_analysis.dad_data3d
    testgradient = mocca_gradient(dad_object)  # Change to actual background file
    background_gradient = load_background_chr_dad(full_analysis.info["LC Method"], settings)
    exp = HplcInput(testpath, background_gradient)
    exp.custom_data = dad_object
    trial_campaign.add_hplc_input(exp)
    settings_trial = Settings('custom',  # dont give wl high or low pass, not implemented for custom
                        absorbance_threshold=1500,  # default is 500, doesn't work well for us, could optimize this (also no gradient yet..)
                        peaks_high_pass=1, peaks_low_pass=50,  # limits in retention time for peak detection.
                        spectrum_correl_thresh=0.99, relative_distance_thresh=0.0025)
    trial_campaign.process_all_hplc_input(settings_trial)
    # report(trial_campaign, report_path)  # gives error, could be from the missing gradient

    dad_peaklist = []
    for peak in trial_campaign.peak_db:
        mocca_index_to_time(peak.maximum, dad_object)
        if not peak.compound_id == None:  # To ensure that PARAFAC peaks replace the initial peak
            dad_peaklist.append(dad_peak(peak, dad_object))
    total_dad_integral = 0
    for peak in dad_peaklist:
        total_dad_integral += peak.integral  # For later calculation of relative abs

    for peak in dad_peaklist:
        integral_sum = 0
        for second_peak in dad_peaklist:
            if peak.left <= second_peak.right and peak.right >= second_peak.left:
                integral_sum += second_peak.integral
        peak.relative_purity = peak.integral/integral_sum
        if peak.integral == integral_sum:
            peak.pure = True
        else:
            peak.pure = False
        peak.relative_integral = peak.integral/total_dad_integral

    return dad_peaklist

def mocca_index_to_time(index, dad_object):
    """
    Takes the index for peaks as given from mocca and returns the corresponding time value in seconds.
    :return:
    """
    time = index * dad_object["time"][-1] /len(dad_object["time"])
    return time

"""
File for specific data processing of DAD data.
"""

class DAD_Spectra:
    def __init__(self, data, info):
        self.data = data
        self.info= info

    def normalize(self, normalization_factor = 1):
        """
        Fct to normalize the spectra.
        :param normalization_factor: Value the spectra should be normalized to. Default is set to 1.
        :return: object with normalized peak intensities
        """
        max_intensity = np.max(self.data["Intensity"])/normalization_factor
        self.data['Intensity'] = self.data['Intensity'].astype(float) # Necessary change, was integer before
        for i in range(len(self.data["Intensity"])):
            self.data.loc[i, "Intensity"] /= max_intensity
        self.info["Normalized"] = True
        return

    def full_processing_dad(self, background_signal, settings):
        """
        Apply all processing fct for dad spectra.
        :param settings:
        :return:
        """
        self.raw_data = self.data.copy()
        self.background_subtraction(background_signal, settings)
        self.smoothing_function()
        self.normalize()
        self.info["Processed"] = True
        return

    def plot_dad(self):
        """
        Simple fct to plot the spectra.
        :return:
        """
        plt.plot(self.data["Wavelength"], self.data["Intensity"])
        plt.show()
        return

    def background_subtraction(self, background_signal, settings):
        """
        Subtract background spectrum from not normalized spectrum.
        :param background_spectra: Averaged background spectra, not normalized
        :return:
        """
        background_spectra = background_signal.dad_spectrum

        merged_data = pd.merge(self.data, background_spectra.data, on="Wavelength", how="outer", suffixes=("_a", "_b"))
        merged_data["Intensity"] = merged_data["Intensity_a"] - merged_data["Intensity_b"]
        merged_data.drop(["Intensity_a", "Intensity_b"], axis=1, inplace=True)
        merged_data.loc[merged_data["Intensity"] < 0, "Intensity"] = 0
        self.info["Background correction"] = True
        self.data = merged_data
        self.info["Normalized"] = False
        return

    def smoothing_function(self):
        """
        Savitzky-Golay as standard filter.
        """
        # window length 25 means 8.3 nm (at .3nm resolution) for averaging window. Seems reasonable, uv peaks are broad anyway.
        self.data["Intensity"] = savgol_filter(self.data["Intensity"], window_length=25, polyorder=3)
        return

class DAD_full_chr:
    def __init__(self, data ,chr_info):
        self.data = data
        self.chr_info = chr_info

    def extract_single_dad(self, time):
        """
            Takes in a full chromatogram and gives out the DAD spectra for a certain point in time
            :param full_chr_dad: Needs to be array in array type object
            :param time: time value in seconds
            :return: MS object as pandas dataframe
            """
        time = float(time)
        timespan = np.array(self.data["time"])
        if time > timespan[-1]:
            print("No DAD spectra for t>=" + str(time))
            return
        index = np.argmin(np.abs(timespan - time))  # gives the index value of the closest point in time
        wavelength_values = self.data["DAD spectra"][index, 0]
        intensity_values = self.data["DAD spectra"][index, 1]
        extracted_dad = pd.DataFrame({
            "Wavelength": wavelength_values,
            "Intensity": intensity_values,
        })

        return DAD_Spectra(extracted_dad, init.import_spectra_info(time=time))

    def extract_dad_timespan(self, time_0, time_end):
        """
        Extract the averaged spectrum of a given timespan between time_0 and time_end.
        Both t values are still included.
        :param time_0:
        :param time_end:
        :return: single dad_object with the averaged spectrum. The time value is set as the middle.
        """
        time_0 = float(time_0)
        time_end = float(time_end)
        timespan = np.array(self.data["time"])
        index_0 = np.argmin(np.abs(timespan - time_0))
        index_end = np.argmin((np.abs(timespan - time_end)))

        wavelength_values = self.data["DAD spectra"][index_0, 0]
        intensity_values = self.data["DAD spectra"][index_0, 1]
        extracted_dad = pd.DataFrame({
            "Wavelength": wavelength_values,
            "Intensity": intensity_values,
        })

        for i in range(index_0+1, index_end+1, 1):
            wavelength_values = self.data["DAD spectra"][i, 0]
            intensity_values = self.data["DAD spectra"][i, 1]
            extracted_dad_next = pd.DataFrame({
                "Wavelength": wavelength_values,
                "Intensity": intensity_values,
            })
            extracted_dad = pd.merge(extracted_dad, extracted_dad_next, on="Wavelength", how="outer", suffixes=("_a", "_b"))
            extracted_dad.fillna(0, inplace=True)
            extracted_dad["Intensity"] = extracted_dad["Intensity_a"] + extracted_dad["Intensity_b"]
            extracted_dad.drop(["Intensity_a", "Intensity_b"], axis=1, inplace=True)
        extracted_dad["Intensity"] = extracted_dad["Intensity"]/(index_end - index_0 +1) # Norms the spectrum intensity

        return DAD_Spectra(extracted_dad, init.import_spectra_info(time = ((time_end - time_0)/2), run_nr=self.chr_info["Run Nr"]))

def load_background_chr_dad(method, settings):
    """
    Searches for the background file corresponding to a method, returns signal object.
    :param:
    :return:
    """
    background_folder = os.path.join(settings["directory_project"], "Data_examples", "background_spectra")
    background_list = os.listdir(background_folder)
    background_chr_dad = None

    for background_file_name in background_list:
        background_file_path = os.path.join(background_folder, background_file_name)
        if os.path.isfile(background_file_path) and background_file_path.lower().endswith(".cdf"):
            if background_file_name[:-4] == method:  # -4 to remove .cdf ending!
                background_chr_dad = extract_background_chr_dad(background_file_path)

    if background_chr_dad == None:
        print("No backgroundfile exists for the method " + method)
        return
    return background_chr_dad

def extract_background_chr_dad(background_file_path):
    background_file = Dataset(background_file_path, "r")

    dad_data = np.empty((len(
        background_file.groups["DAD Data"].groups["Chromatogram DAD Data"].variables["Intensity"][:,0].compressed()), 2),
        dtype=object)
    dad_wavelengths = np.array(
        background_file.groups["DAD Data"].groups["Chromatogram DAD Data"].variables["Wavelength"][:].compressed())  # save wavelengths (no time dependency)

    # load all data into pd dataframe for processing
    for i in range(len(
            background_file.groups["DAD Data"].groups["Chromatogram DAD Data"].variables["Intensity"][:,0].compressed())):
        dad_data[i,0] = dad_wavelengths
        dad_data[i,1] = background_file.groups["DAD Data"].groups["Chromatogram DAD Data"].variables["Intensity"][i,:].compressed()
    dad_spectrum = np.array(dad_data)
    full_chromatogram_dad = {
        "time": background_file.groups["DAD Data"].groups["Chromatogram DAD Data"].variables["Time"][:].compressed(),
        "DAD spectra": dad_spectrum,  # 0 to access wavelengths and 1 to access intensity values
    }

    background_file.close()
    background_chr_dad = mocca_gradient(full_chromatogram_dad)

    return background_chr_dad


# Friday combine everything with mocca peak picking. Then implement MS chromatogram subtraction and peak picking there.
# Finally, develop peak deconvolution for MS peaks
# Modify peak deconvolution in mocca
# Put everything together and iterate over multiple runs
# Write output files.
# Write fitting fct for dad spectra if comparison insufficient
# Write optimatization script
# Make sure MOCCA background takes in actual background file for the method indicated!!!
# Verfiy visusally that data extraction (also PARAFAC) to peak object works

class mocca_gradient:
    def __init__(self, data):
        self.data = data
        self.path = None

class dad_peak:
    def __init__(self, peak, dad_object):
        self.max = mocca_index_to_time(peak.maximum, dad_object)
        self.left = mocca_index_to_time(peak.left, dad_object)
        self.right = mocca_index_to_time(peak.right, dad_object)
        self.pure = peak.pure
        self.integral = peak.integral
        if type(peak.dataset) == ParafacData:
            self.parafac = True
        else:
            self.parafac = False

        """
        Extract spectrum and give averaged spectrum as property of the peak.
        """
        time_0 = float(self.left)
        time_end = float(self.right)
        timespan = peak.dataset.time *60  # convert minutes to secounds
        index_0 = np.argmin(np.abs(timespan - time_0))
        index_end = np.argmin((np.abs(timespan - time_end)))

        intensity_values = np.zeros(len(peak.dataset.data))
        for i in range(len(peak.dataset.data)):
            intensity_values[i] += np.sum(peak.dataset.data[i, index_0:index_end])

        extracted_dad = pd.DataFrame({
            "Wavelength": peak.dataset.wavelength,
            "Intensity": intensity_values,
        })
        extracted_dad["Intensity"] = extracted_dad["Intensity"]/(index_end - index_0 +1) # Norms the spectrum intensity

        self.data = DAD_Spectra(extracted_dad, init.import_spectra_info(time=self.max))  # Fill info later


    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))