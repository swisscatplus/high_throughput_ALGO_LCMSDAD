import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import os
from decimal import Decimal, ROUND_HALF_UP
from netCDF4 import Dataset

import initialize as init
import create_file
import data_processing_dad as dpr_dad
import ms_deconvolution as msd

def process_deconvoluted_peak(ms_spectrum, dad_spectrum, method, settings):
    """
    Take aligned spectra and create signal file based upon this.
    :param ms_spectrum:
    :param dad_spectrum:
    :param settings:
    :return:
    """
    background_signal = load_background_spectra_signal(method,
                                                       settings)  # later give as input, requires preprocessed dtb
    ms_spectrum.full_processing_ms(background_signal, settings)
    dad_spectrum.full_processing_dad(background_signal, settings)

    create_file.create_signal_file(ms_spectrum, dad_spectrum, settings["directory_project"])
    return

def process_found_signal(full_analysis, peak_info, settings):
    """
    Processing of a peak into a signal file.
    For manual entry of peaks based of times within peak_info.
    :return:
    """
    background_signal = load_background_spectra_signal(full_analysis.info["LC Method"],
                                                       settings)  # later give as input, requires preprocessed dtb
    ms_chr = MS_full_chr(full_analysis.ms_data3d, full_analysis.info)
    extracted_ms = ms_chr.extract_ms_timespan(peak_info["start_time"], peak_info["end_time"])
    extracted_ms.full_processing_ms(background_signal, settings)

    dad_chr = dpr_dad.DAD_full_chr(full_analysis.dad_data3d, full_analysis.info)

    extracted_dad = dad_chr.extract_dad_timespan(2.74 * 60, 2.84 * 60)
    extracted_dad.full_processing_dad(background_signal, settings)

    create_file.create_signal_file(extracted_ms, extracted_dad, settings["directory_project"])

    return


class MS_Spectra:
    """
    Class for processing of all MS_Spectra
    """
    def __init__(self, data, info):
        self.data = data
        self.info = info
        # fct to ensure that m/z is given with correct decimals, avoids floating-point inaccuracies
        self.data["m/z"] = self.data["m/z"].apply(lambda mass: float(Decimal(str(mass)).quantize(Decimal("0.1"))))
        # lambda is mock fct as input, conversion first to str better for Decimal, back to float bcs of pd dt frame


    def normalize(self, normalization_factor = 1):
        """
        Fct to normalize the spectra.
        :param normalization_factor: Value the spectra should be normalized to. Default is set to 1.
        :return: object with normalized peak intensities
        """
        total_abundance = self.data["Intensity"].sum() / normalization_factor
        self.data['Intensity'] = self.data['Intensity'].astype(float) # Necessary change, was integer before
        for i in range(len(self.data["Intensity"])):
            self.data.loc[i, "Intensity"] /= total_abundance
        self.info["Normalized"] = True
        return

    def noise_removal_threshold(self, settings):
        """
        Fct to remove noise from the spectra. All signals with an intensity lower than the threshold relative to the
        largest peak are being zero-filled.
        :param threshold: min intensity relative to max peak
        :return: Object with noice spectra removed.
        """

        threshold = settings["threshold_ms_spectra"]
        max_intensity = self.data["Intensity"].max()
        for i in range(len(self.data["Intensity"])):
            if self.data["Intensity"][i] < threshold*max_intensity:
                self.data.loc[i, "Intensity"] = 0
        self.info["Normalized"] = False
        self.remove_zeros()

        return

    def signal_smooving(self):
        """  Not sure if this is a good idea... Test!
        Take an average of signal in a given bandwidth to account for low precision spectra.
        :return: The spectra with newly averaged signals.
        """
        max_mass_value = self.data["m/z"].max()
        for mass in self.data["m/z"]:
            intensity = np.zeros(3)
            mass_vector = np.array([(mass-0.1), mass, (mass+0.1)])
            # Convert to correct decimal, avoiding floating-point error
            mass_vector = np.array(
                [float(Decimal(str(m)).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP)) for m in mass_vector])
            index = np.full(3, np.nan)
            for i in range(3):
                try:
                    index[i] = int(self.data.loc[self.data["m/z"] == mass_vector[i]].index[0])
                    intensity[i] = self.data.loc[index[i], "Intensity"]
                except IndexError:
                    intensity[i] = 0
                    # print("Error - no corresponding mass" + str(mass-0.1+0.1*i))
            total_intensity = np.sum(intensity)
            # now find the mass with the highest intensity
            mass_max_index = np.argmax(intensity)

            # then give it the complete intensity and set all others to zero
            for i in range(3):
                if i == mass_max_index and not np.isnan(index[i]):
                    self.data.loc[index[i], "Intensity"] = total_intensity
                else:
                    if mass_vector[i] <= max_mass_value and not np.isnan(index[i]):  # changed to max_mass_value from last df index
                         try:
                            self.data.loc[index[i], "Intensity"] = 0
                         except IndexError:  # Shouldn't actually occur
                            pass
        self.remove_zeros()
        return

    def plot_spectra(self):
        """
        Simple plot of the MS spectra.
        :return:
        """
        plt.vlines(self.data["m/z"], ymin = 0, ymax = self.data["Intensity"])
        plt.show()
        return

    def full_processing_ms(self, background_signal, settings):
        """
        Simple command for full standardised processing, returns spectrum for data base comparison.
        :param threshold_noise:
        :return:
        """
        self.raw_data = self.data.copy()
        self.noise_removal_threshold(settings)
        self.normalize()
        self.signal_smooving()
        self.background_subtraction(background_signal, settings)
        self.normalize()
        self.weighting_fct(settings)
        self.info["Processed"] = True
        return

    def weighting_fct(self, settings):
        """
        Apply a weighting function to the intensities of a given ms distributiom
        :param spectrum: The spectra as an MS_Spectra object
        :return: Spectrum with additional weighted intensity in the data
        ms_weighting_fct is predefined in settings (main.py)
        Types of weighting functions:
            exponential: Weights the (normalized) intensity by I_new = I^0.25
        """
        ms_weighting_fct = settings["ms_weighting_fct"]

        if ms_weighting_fct == "exponential":
            for i in range(len(self.data["Intensity"])):
                self.data.loc[i, "Intensity"] = self.data.loc[i, "Intensity"] ** 0.5
                self.info["Normalized"] = False
                self.info["Weight Fct"] = ms_weighting_fct
        if ms_weighting_fct == "exponential2":
            for i in range(len(self.data["Intensity"])):
                self.data.loc[i, "Intensity"] = self.data.loc[i, "Intensity"] * (math.e ** (-self.data.loc[i, "Intensity"]))
                self.info["Normalized"] = False
                self.info["Weight Fct"] = ms_weighting_fct
        if ms_weighting_fct == "logarithmic":
            for i in range(len(self.data["Intensity"])):
                if self.data.loc[i, "Intensity"] > 0:
                    self.data.loc[i, "Intensity"] = math.log(self.data.loc[i, "Intensity"])
                else:
                    pass
                self.info["Normalized"] = False
                self.info["Weight Fct"] = ms_weighting_fct
        if ms_weighting_fct == "sin":
            for i in range(len(self.data["Intensity"])):
                self.data.loc[i, "Intensity"] = math.sin(self.data.loc[i, "Intensity"] * math.pi/2)
                self.info["Normalized"] = False
                self.info["Weight Fct"] = ms_weighting_fct
        if ms_weighting_fct == "sin2":
            for i in range(len(self.data["Intensity"])):
                self.data.loc[i, "Intensity"] = math.sin(self.data.loc[i, "Intensity"] * math.pi/2) ** 2
                self.info["Normalized"] = False
                self.info["Weight Fct"] = ms_weighting_fct
        else:
            self.info["Weight Fct"] = ms_weighting_fct

        if ms_weighting_fct == "None" and self.info["Normalized"] == True:
            pass
        else:
            self.normalize()
        return

    def background_subtraction(self, background_signal, settings):
        """
        Substract background spectrum from not normalized spectrum.
        :param background_spectra: Averaged background spectra, not normalized
        :return:
        """
        background_spectra = background_signal.ms_spectrum

        merged_data = pd.merge(self.data, background_spectra.data, on="m/z", how="outer", suffixes=("_a", "_b"))
        merged_data.fillna(0, inplace=True)
        merged_data["Intensity"] = merged_data["Intensity_a"] - merged_data["Intensity_b"]
        merged_data.drop(["Intensity_a", "Intensity_b"], axis=1, inplace=True)
        merged_data.loc[merged_data["Intensity"] < 0, "Intensity"] = 0
        self.info["Background correction"] = True
        self.data = merged_data
        self.info["Normalized"] = False
        self.remove_zeros()
        return

    def remove_zeros(self):
        """
        Fct to remove all zero-valued m/z entries, to increaase comp. speed
        Resets index to start from zero again
        :return:
        """
        self.data = self.data[self.data["Intensity"] != 0]
        self.data.reset_index(drop = True, inplace = True)
        return


class MS_full_chr:
    """

    """
    def __init__(self, data, info):
        self.data = data
        self.chr_info = info

    def extract_single_ms(self, time):
        """
            Takes in a full chromatogram and gives out the MS spectra for a certain point in time
            :param full_chr_ms: Needs to be array in array type object
            :param time: time value in seconds
            :return: MS object as pandas dataframe
            """
        time = float(time)
        timespan = np.array(self.data["time"])
        if time > timespan[-1]:
            print("No MS spectra for t>=" + str(time))
            return
        index = np.argmin(np.abs(timespan - time)) # gives the index value of the closest point in time
        mass_values = self.data["MS spectra"][index, 0]
        intensity_values = self.data["MS spectra"][index, 1]
        extracted_ms = pd.DataFrame({
            "m/z": mass_values,
            "Intensity": intensity_values,
        })

        return MS_Spectra(extracted_ms, init.import_spectra_info(time = time))

    def extract_ms_timespan(self, time_0, time_end):
        """
        Extract the averaged spectrum of a given timespan between time_0 and time_end.
        Both t values are still included.
        :param time_0:
        :param time_end:
        :return: single ms_object with the averaged spectrum. The time value is set as the middle.
        """
        time_0 = float(time_0)
        time_end = float(time_end)

        timespan = np.array(self.data["time"])
        index_0 = np.argmin(np.abs(timespan - time_0))
        index_end = np.argmin((np.abs(timespan - time_end)))

        mass_values = self.data["MS spectra"][index_0, 0]
        intensity_values = self.data["MS spectra"][index_0, 1]
        extracted_ms = pd.DataFrame({
            "m/z": mass_values,
            "Intensity": intensity_values,
        })

        for i in range(index_0+1, index_end+1, 1):
            mass_values = self.data["MS spectra"][i, 0]
            intensity_values = self.data["MS spectra"][i, 1]
            extracted_ms_next = pd.DataFrame({
                "m/z": mass_values,
                "Intensity": intensity_values,
            })
            extracted_ms = pd.merge(extracted_ms, extracted_ms_next, on="m/z", how="outer", suffixes=("_a", "_b"))
            extracted_ms.fillna(0, inplace=True)
            extracted_ms["Intensity"] = extracted_ms["Intensity_a"] + extracted_ms["Intensity_b"]
            extracted_ms.drop(["Intensity_a", "Intensity_b"], axis=1, inplace=True)
        extracted_ms["Intensity"] = extracted_ms["Intensity"]/(index_end - index_0 +1)  # Norms the spectrum intensity

        return MS_Spectra(extracted_ms, init.import_spectra_info(time = ((time_end - time_0)/2), run_nr=self.chr_info["Run Nr"]))

    def extract_ms_timespan_values(self, time_0, time_end, mass_ranges_decon):
        time_0 = float(time_0)
        time_end = float(time_end)

        timespan = np.array(self.data["time"])
        index_0 = np.argmin(np.abs(timespan - time_0))
        index_end = np.argmin((np.abs(timespan - time_end)))

        accepted_masses = []
        for mass_range in mass_ranges_decon:
            for mass in np.arange(mass_range-0.25, mass_range+3.25, 0.1):  # +3.25 to account for low intensity isotope peaks
                # Low intensity isotope peaks often not detected in deconvolution but survive secondary background subt. and noise removal
                accepted_masses.append(round(mass, 1))  # round to avoid floating point error

        mass_values = self.data["MS spectra"][index_0, 0]
        intensity_values = self.data["MS spectra"][index_0, 1]

        mass_filter_mask = [mass in accepted_masses for mass in mass_values]

        extracted_ms = pd.DataFrame({
            "m/z": np.array(mass_values)[mass_filter_mask],  # changed to putting in np.array directly, shouldn't change df structure
            "Intensity": np.array(intensity_values)[mass_filter_mask],
        })

        for i in range(index_0 + 1, index_end + 1, 1):
            mass_values = self.data["MS spectra"][i, 0]
            intensity_values = self.data["MS spectra"][i, 1]

            mass_filter_mask = [mass in accepted_masses for mass in mass_values]

            extracted_ms_next = pd.DataFrame({
                "m/z": np.array(mass_values)[mass_filter_mask],
                "Intensity": np.array(intensity_values)[mass_filter_mask],
            })
            extracted_ms = pd.merge(extracted_ms, extracted_ms_next, on="m/z", how="outer", suffixes=("_a", "_b"))
            extracted_ms.fillna(0, inplace=True)
            extracted_ms["Intensity"] = extracted_ms["Intensity_a"] + extracted_ms["Intensity_b"]
            extracted_ms.drop(["Intensity_a", "Intensity_b"], axis=1, inplace=True)
        extracted_ms["Intensity"] = extracted_ms["Intensity"] / (
                    index_end - index_0 + 1)  # Norms the spectrum's intensity

        return MS_Spectra(extracted_ms,
                          init.import_spectra_info(time=((time_end - time_0) / 2), run_nr=self.chr_info["Run Nr"]))

def create_new_background_spectra(background_filepath, method_name, settings):
    """
    Goes through all background_spectra files (later). Takes the one with the correct method and loads it as an
    chr object.
    :param background_folder:
    :return: Averaged Background spectrum.
    """

    background_full_analysis = init.import_run_json(background_filepath, method=method_name)

    background_chr_ms = MS_full_chr(background_full_analysis.ms_data3d\
                                    , init.import_full_chr_info(background=True))
    background_chr_ms.chr_info["LC Method"] = method_name  # Later parse this directly is possible (?)

    if background_chr_ms.chr_info["Background"] == False:
        print("The file " + background_filepath + " is not a background spectrum.")
        return

    background_folder = os.path.join(settings["directory_project"], "Data_examples", "background_spectra")
    file_name = os.path.join(background_folder, (background_chr_ms.chr_info["LC Method"] + ".cdf"))

    background_list = os.listdir(background_folder)
    if (background_chr_ms.chr_info["LC Method"] + ".cdf") in background_list:
        print("The method " + background_chr_ms.chr_info["LC Method"] + " already exists!")
        return

    background_masses_list = msd.determine_background_masses_list(background_full_analysis, background_chr_ms, settings)

    if background_full_analysis.info["plus_minus_acq"]:
        background_chr_ms = MS_full_chr(background_full_analysis.ms_data3d_polarized \
                                        , init.import_full_chr_info(background=True))
        background_chr_ms.chr_info["LC Method"] = method_name

    begin_t_msr = 0 #  later define this as the difference in measurment time, until fraction arrives at ms (?)
    end_t_msr = background_chr_ms.data["time"][-1]
    background_ms = background_chr_ms.extract_ms_timespan(begin_t_msr, end_t_msr)  # Change data depending on polarization!
    background_ms.raw_data = background_ms.data.copy()

    background_ms.noise_removal_threshold(settings)
    background_ms.normalize()
    background_ms.signal_smooving()
    background_ms.remove_zeros()

    background_chr_dad = dpr_dad.DAD_full_chr(background_full_analysis.dad_data3d \
                                              , init.import_full_chr_info(background=True))
    begin_t_dad = 0
    end_t_dad = background_chr_dad.data["time"][-1]
    background_dad = background_chr_dad.extract_dad_timespan(begin_t_dad, end_t_dad)
    background_dad.raw_data = background_dad.data.copy()

    # possibly pre-process dad background data here

    background_file = Dataset(file_name, 'w', format='NETCDF4')

    background_file.Method = background_chr_ms.chr_info["LC Method"]

    intensity = background_file.createDimension("intensity", None)
    ft_wvl = background_file.createDimension("ft_wvl", None)
    time_dimension = background_file.createDimension("time", None)

    grp_ms = background_file.createGroup("MS Data")
    grp_dad = background_file.createGroup("DAD Data")

    mass = grp_ms.createDimension("mass", None)
    mass_list_dim = grp_ms.createDimension("Background masses dim", len(background_masses_list))
    mass_list = grp_ms.createVariable("Background masses list", "f8", mass_list_dim)
    mass_list[:] = background_masses_list[:]
    grp_raw_ms = grp_ms.createGroup("Raw MS Data")
    raw_mass = grp_raw_ms.createVariable("mass", "f8", mass)
    raw_ms_intensity = grp_raw_ms.createVariable("Intensity", "f8", intensity)
    raw_mass[:] = background_ms.raw_data.loc[:, "m/z"]
    raw_ms_intensity[:] = background_ms.raw_data.loc[:, "Intensity"]

    grp_proc_ms = grp_ms.createGroup("Processed MS Data")
    proc_mass = grp_proc_ms.createVariable("mass", "f8", mass)
    proc_ms_intensity = grp_proc_ms.createVariable("Intensity", "f8", intensity)
    proc_mass[:] = background_ms.data.loc[:, "m/z"]
    proc_ms_intensity[:] = background_ms.data.loc[:, "Intensity"]

    wavelength = grp_dad.createDimension("Wavelength", len(background_chr_dad.data["DAD spectra"][0,0]))
    ft_wvl = grp_dad.createDimension("ft_wvl", None)
    intensity_dad_chr = grp_dad.createDimension("Intensity_dad", len(background_chr_dad.data["DAD spectra"][:]))

    grp_raw_dad = grp_dad.createGroup("Raw DAD Data")
    raw_wvl = grp_raw_dad.createVariable("Wavelength", "f8", wavelength)
    raw_dad_intensity = grp_raw_dad.createVariable("Intensity", "f8", intensity)
    raw_wvl[:] = background_dad.raw_data.loc[:, "Wavelength"]
    raw_dad_intensity[:] = background_dad.raw_data.loc[:, "Intensity"]

    grp_proc_dad = grp_dad.createGroup("Processed DAD Data")
    proc_wvl = grp_proc_dad.createVariable("Wavelength", "f8", wavelength)
    proc_dad_intensity = grp_proc_dad.createVariable("Intensity", "f8", intensity)
    proc_ft_wvl = grp_proc_dad.createVariable("FT Wavelength", "f8", ft_wvl)  # delete if ft comparison not working well
    proc_wvl[:] = background_dad.data.loc[:, "Wavelength"]
    proc_dad_intensity[:] = background_dad.data.loc[:, "Intensity"]
    # proc_ft_wvl[:] = dad_spec.data.loc[:, "FT Wavelength"]  # Uncomment once FT has been constructed.

    grp_chr_dad = grp_dad.createGroup("Chromatogram DAD Data")
    chr_dad_time = grp_chr_dad.createVariable("Time", "f8", time_dimension)
    chr_intensity = grp_chr_dad.createVariable("Intensity", "f8", (intensity_dad_chr, wavelength))
    chr_wavelength = grp_chr_dad.createVariable("Wavelength", "f8", wavelength)

    chr_dad_time[:] = background_chr_dad.data["time"][:]
    chr_wavelength[:] = background_chr_dad.data["DAD spectra"][0,0]  # wavelength array shouldn't change with time
    for i in range(len(background_chr_dad.data["DAD spectra"])):
        chr_intensity[i,:] = background_chr_dad.data["DAD spectra"][i,1]  # wavelength corresponding by index with chr_wavelength

    background_file.close()
    return

def load_background_spectra_signal(method, settings, processed=True):
    """
    Searches for the background file corresponding to a method, returns signal object.
    :param:
    :return:
    """
    background_folder = os.path.join(settings["directory_project"], "Data_examples", "background_spectra")
    background_list = os.listdir(background_folder)
    background_signal = None

    for background_file_name in background_list:
        background_file_path = os.path.join(background_folder, background_file_name)
        if os.path.isfile(background_file_path) and background_file_path.lower().endswith(".cdf"):
            if background_file_name[:-4] == method:  # -4 to remove .cdf ending!
                background_signal = extract_background_signal(background_file_path, processed)
        else:
            print("!!! File " + background_file_name + " is not a supported background file type!")

    if background_signal == None:
        print("No backgroundfile exists for the method " + method)
        return
    return background_signal


def extract_background_signal(background_file_path, processed):
    background_file = Dataset(background_file_path, "r")
    if processed:
        data_ms = pd.DataFrame({
            "m/z": background_file.groups["MS Data"].groups["Processed MS Data"].variables["mass"][:].compressed(),
            "Intensity": background_file.groups["MS Data"].groups["Processed MS Data"].variables["Intensity"][
                         :].compressed(),
        })
        data_dad = pd.DataFrame({
            "Wavelength": background_file.groups["DAD Data"].groups["Processed DAD Data"].variables["Wavelength"][:]\
                .compressed(),
            "Intensity": background_file.groups["DAD Data"].groups["Processed DAD Data"].variables["Intensity"][:]\
                .compressed()
        })
    else:
        data_ms = pd.DataFrame({
            "m/z": background_file.groups["MS Data"].groups["Raw MS Data"].variables["mass"][:].compressed(),
            "Intensity": background_file.groups["MS Data"].groups["Raw MS Data"].variables["Intensity"][:].compressed(),
        })
        data_dad = pd.DataFrame({
            "Wavelength": background_file.groups["DAD Data"].groups["Raw DAD Data"].variables["Wavelength"][:] \
                .compressed(),
            "Intensity": background_file.groups["DAD Data"].groups["Raw DAD Data"].variables["Intensity"][:] \
                .compressed()
        })

    info = {
        "Processed": processed,
        "Methods": background_file.Method,
    }

    ms_spectrum = MS_Spectra(data_ms, init.import_spectra_info())
    dad_spectrum = dpr_dad.DAD_Spectra(data_dad, init.import_spectra_info())
    background_file.close()
    background_signal = Background(ms_spectrum, dad_spectrum, info)

    return background_signal

class Background:
    """
      Class for a Signal, containing the information on the method, and the processed DAD and MS spectra and
      the Rf value.
      Used when dealing with either peaks or signals.
      """
    # Consider loading more fct into this class..
    def __init__(self, ms_spectrum, dad_spectrum, info):
        self.ms_spectrum = ms_spectrum
        self.dad_spectrum = dad_spectrum
        self.info = info

        # Not necessary to correct for decimals since this is in the initalization of ms_spectra..