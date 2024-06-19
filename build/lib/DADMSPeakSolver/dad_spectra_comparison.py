import pandas as pd
import numpy as np
import scipy as sc
import math
import matplotlib.pyplot as plt

from . import initialize as init
from . import data_processing
from . import data_processing_dad as dpr_dad

def comparison_dad(spectra_1, spectra_2, settings):
    """
    Fct to use all comparison fct for dad spectra depending on the settings.
    :return:
    """
    # set all comparison booleans to False
    pearson_comp = False
    dot_product_comp = False
    ft_comp = False
    derivative_comp = False

    if not spectra_1.info["Normalized"]:
        spectra_1.normalize()
    if not spectra_2.info["Normalized"]:
        spectra_2.normalize()

    if settings["ms_algorithm"] == "all":
        if pearson_similarity(spectra_1, spectra_2) < settings["threshold_dad_pearson"]:
            print("R^2 similarity shows high DAD similarity: " +\
                  str(pearson_similarity(spectra_1, spectra_2)) +\
                  " with " + spectra_2.info["Molecule"])
            pearson_comp = True
        if dot_product_similarity(spectra_1, spectra_2) < settings["threshold_dad_dot_product"]:
            print("Dot product similarity shows high DAD similarity: " +\
                  str(dot_product_similarity(spectra_1, spectra_2)) +\
                  " with " + spectra_2.info["Molecule"])
            dot_product_comp = True
        if ft_similarity(spectra_1, spectra_2, settings) < settings["threshold_dad_ft"]:
            print("FT similarity shows high DAD similarity: " +\
                  str(ft_similarity(spectra_1, spectra_2, settings)) +\
                  " with " + spectra_2.info["Molecule"])
            ft_comp = True
        if derivative_similarity(spectra_1, spectra_2, settings) < settings["threshold_dad_derivative"]:
            print("Derivative similarity shows high DAD similarity: " +\
                  str(derivative_similarity(spectra_1, spectra_2, settings)) +\
                  " with " + spectra_2.info["Molecule"])
            derivative_comp = True

        if pearson_comp and dot_product_comp and ft_comp and derivative_comp:
            return True
        else:
            return False
    elif settings["dad_algorithm"] == "pearson":
        if pearson_similarity(spectra_1, spectra_2) < settings["threshold_dad_pearson"]:
            print("R^2 similarity shows high DAD similarity: " +\
                  str(pearson_similarity(spectra_1, spectra_2)) +\
                  " with " + spectra_2.info["Molecule"])
            pearson_comp = True
        return pearson_comp
    elif settings["dad_algorithm"] == "dot_product":
        if dot_product_similarity(spectra_1, spectra_2) < settings["threshold_dad_dot_product"]:
            print("Dot product similarity shows high DAD similarity: " +\
                  str(dot_product_similarity(spectra_1, spectra_2)) +\
                  " with " + spectra_2.info["Molecule"])
            dot_product_comp = True
        return dot_product_comp
    elif settings["dad_algorithm"] == "ft":
        if ft_similarity(spectra_1, spectra_2, settings) < settings["threshold_dad_ft"]:
            print("FT similarity shows high DAD similarity: " +\
                  str(ft_similarity(spectra_1, spectra_2, settings)) +\
                  " with " + spectra_2.info["Molecule"])
            ft_comp = True
        return ft_comp
    elif settings["dad_algorithm"] == "derivative":
        if derivative_similarity(spectra_1, spectra_2, settings) < settings["threshold_dad_derivative"]:
            print("Derivative similarity shows high DAD similarity: " +\
                  str(derivative_similarity(spectra_1, spectra_2, settings)) +\
                  " with " + spectra_2.info["Molecule"])
            derivative_comp = True
        return derivative_comp
    else:
        print(settings["dad_algorithm"] + " is an invalid DAD comparison algorithm.")
        return False


def pearson_similarity(spectra_1, spectra_2):
    """
    Comparison algorithm based on the pearson similarity of two spectra.
    If > 0.1 (in this metric) would have been impure peak on mocca.
    Known not to be super precise, typically unimodal changes over retention profile observed (mocca).
    Squared to get the R^2 value.
    :param spectra_1:
    :param spectra_2:
    :return:
    """
    pearson_coeff = np.corrcoef(spectra_1.data["Intensity"], spectra_2.data["Intensity"])[0,1]**2  # obtain off-diagonal matrix element (=coefficient)
    return 1-abs(pearson_coeff)

def dot_product_similarity(spectra_1, spectra_2):
    """
    Computes the cos of the spectral vectors in analogy to the MS algorithm.
    :param spectra_1:
    :param spectra_2:
    :return:
    """
    dotproduct = np.sum(spectra_1.data["Intensity"]*spectra_2.data["Intensity"])  # dot product
    sq_cos_score = dotproduct ** 2 / (np.sum(spectra_1.data["Intensity"] ** 2) * \
                                      np.sum(spectra_2.data["Intensity"] ** 2))
    return 1 - np.sqrt(sq_cos_score)

def ft_similarity(spectra_1, spectra_2, settings):
    """
    Comparison algorithm based on the FT of the spectrum and a dot product in the frequency domain.
    Settings to decide where to cut in the frequency domain to avoid noice and background.
    :param spectra_1:
    :param spectra_2:
    :return:
    """
    fft_intensities_1 = np.fft.fft(spectra_1.data["Intensity"])  # Obtain FFT intensity values
    fft_intensities_2 = np.fft.fft(spectra_2.data["Intensity"])

    frequency_rate = spectra_1.data["Wavelength"][1]-spectra_1.data["Wavelength"][0]  # Gets frequency sampling rate from wavelength spacing
    # We get the same number of frequencies as wavelengths [0,701]
    fft_frequencies_1 = np.fft.fftfreq(fft_intensities_1.size, frequency_rate)
    fft_frequencies_2 = np.fft.fftfreq(fft_intensities_2.size, frequency_rate)

    # Remove symmetric half
    upper_half = int(len(fft_intensities_1)/2)
    fft_intensities_1 = fft_intensities_1[settings["FT_low_border"]:upper_half]
    fft_intensities_2 = fft_intensities_2[settings["FT_low_border"]:upper_half]
    fft_frequencies_1 = fft_frequencies_1[settings["FT_low_border"]:upper_half]
    fft_frequencies_2 = fft_frequencies_2[settings["FT_low_border"]:upper_half]

    # Remove beginning and trailing in frequency domain
    fft_intensities_1 = fft_intensities_1[settings["FT_low_border"]:settings["FT_upper_border"]]
    fft_intensities_2 = fft_intensities_2[settings["FT_low_border"]:settings["FT_upper_border"]]
    fft_frequencies_1 = fft_frequencies_1[settings["FT_low_border"]:settings["FT_upper_border"]]
    fft_frequencies_2 = fft_frequencies_2[settings["FT_low_border"]:settings["FT_upper_border"]]

    """plt.vlines(fft_frequencies_2, ymax=np.abs(fft_intensities_2), ymin=0)
    plt.show()"""

    # Now calculate dot product as similarity score between the FFT spectra.
    dotproduct = np.sum(np.abs(fft_intensities_1) * np.abs(fft_intensities_2))  # dot product
    sq_cos_score = dotproduct ** 2 / (np.sum(np.abs(fft_intensities_1) ** 2) * \
                                      np.sum(np.abs(fft_intensities_2) ** 2))
    return 1 - np.sqrt(sq_cos_score)  # Try out how result changes with fitting fct!

# still write the derivative based one
def derivative_similarity(spectra_1, spectra_2, settings):
    """
    Take the nth derivative of the spectra and compare via dot product.
    To enhance fine feature as shoulders in spectrum.
    :param spectra_1:
    :param spectra_2:
    :param settings:
    :return:
    """
    intensities_1 = spectra_1.data["Intensity"]
    intensities_2 = spectra_2.data["Intensity"]
    wavelength_spacing = spectra_1.data["Wavelength"][1] - spectra_1.data["Wavelength"][0]

    for n in range(settings["nth_derivative"]):
        intensities_1 = np.gradient(intensities_1, wavelength_spacing, edge_order = 2)
        intensities_2 = np.gradient(intensities_2, wavelength_spacing, edge_order = 2)
    """plt.plot(spectra_1.data["Wavelength"], intensities_1)
    plt.show()"""

    # Now calculate dot product as similarity score between the derivative spectra.
    dotproduct = np.sum(intensities_1 * intensities_2)  # dot product
    sq_cos_score = dotproduct ** 2 / (np.sum(intensities_1 ** 2) * \
                                      np.sum(intensities_2 ** 2))
    return 1 - np.sqrt(sq_cos_score)

def compare_retention_time(dtb_rt, peak_rt, settings):
    """
    Fct. to check if the retention time of two methods has the same value within an interval.
    :param dtb_rt:
    :param peak_rt:
    :param settings:
    :return:
    """
    difference_rt = abs(dtb_rt - peak_rt)

    return difference_rt <= settings["retention_time_interval"]