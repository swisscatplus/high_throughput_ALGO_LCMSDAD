import pandas as pd
import numpy as np
import scipy as sc
import math

from . import initialize as init
from . import data_processing



def calculate_entropy(spectrum):
    """
    Calculate the Shannon Entropy for a given spectra input (MS_Spectra object)
    :param spectrum:
    :return: The spectral entropy.
    """
    return sc.stats.entropy(spectrum.data["Intensity"])

def entropy_similarity(spectrum_a, spectrum_b):
    """
    Calculates the similarity factor between two spectra based on the (Shannon) search algorithm.
    Small values give good similarity.
    Values are not normalized for the momen ~ need to define max entropy value found.
    :param spectrum_a:
    :param spectrum_b:
    :return:
    """

    merged_spectrum = merge_two_spectra(spectrum_a, spectrum_b)

    # New Spectra is normalized and single intensity coloumns removed
    merged_spectrum.data.drop(["Intensity_a", "Intensity_b"], axis=1, inplace=True)
    merged_spectrum.normalize()

    # Calculate entropy similarity
    return (2 * calculate_entropy(merged_spectrum) - calculate_entropy(spectrum_a) - calculate_entropy(spectrum_b))/math.log(4)

def merge_two_spectra(spectrum_a, spectrum_b):
    """
    Creates spectra with the intensity, fills in non-exsiting m/z in either spectra with zeros, then adds the
    intensity of both spectra.
    :param spectrum_a:
    :param spectrum_b: Background spectrum
    :return: MS Spectra Object of the merged spectra and the individual intensities.
    """
    merged_data = pd.merge(spectrum_a.data, spectrum_b.data, on="m/z", how="outer", suffixes=("_a", "_b"))
    merged_data.fillna(0, inplace=True)
    merged_data["Intensity"] = merged_data["Intensity_a"] + merged_data["Intensity_b"]

    return data_processing.MS_Spectra(merged_data, init.import_spectra_info(merged_spectra=True))


def dot_product_similarity(spectrum_a, spectrum_b):
    """
    Calculates the similarity of two spectra based on the dot product algorithm.
    Small values give good similarity.
    :param spectrum_a: Must be dataframe already
    :param spectrum_b: Must be dataframe already!
    :return:
    """
    merged_spectrum = merge_two_spectra(spectrum_a, spectrum_b)
    dotproduct = np.sum(merged_spectrum.data["Intensity_a"] * merged_spectrum.data["Intensity_b"]) # dot product
    sq_cos_score = dotproduct ** 2 / (np.sum(merged_spectrum.data["Intensity_a"] ** 2) *\
                                      np.sum(merged_spectrum.data["Intensity_b"] ** 2))

    return 1 - np.sqrt(sq_cos_score)

def weighted_dp_similarity(spectrum_a, spectrum_b):
    """
    Calculates the similarity of two spectra based on the weighted dot product algorithm:
    Math. equal to dot product but with I_i = M_i^3 * I_i^0.6
    :param spectrum_a:
    :param spectrum_b:
    :return:
    """
    merged_spectrum = merge_two_spectra(spectrum_a, spectrum_b)
    # Give weigth to intensities
    merged_spectrum.data["Intensity_a"] = merged_spectrum.data["m/z"]**3 * merged_spectrum.data["Intensity_a"]**0.6
    merged_spectrum.data["Intensity_b"] = merged_spectrum.data["m/z"]**3 * merged_spectrum.data["Intensity_b"]**0.6
    # Calculate similarity score
    dotproduct = np.sum(merged_spectrum.data["Intensity_a"] * merged_spectrum.data["Intensity_b"])  # dot product
    sq_cos_score = dotproduct ** 2 / (np.sum(merged_spectrum.data["Intensity_a"] ** 2) * \
                                      np.sum(merged_spectrum.data["Intensity_b"] ** 2))

    return 1 - np.sqrt(sq_cos_score)

def bhattacharya1_distance(spectrum_a, spectrum_b):
    """
    Calculate Bhattacharya distance:
    arccos(sum(sqrt(spectrum_a*spectrum_b)))
    :param spectrum_a:
    :param spectrum_b:
    :return:
    """
    sum = 0
    merged_spectrum = merge_two_spectra(spectrum_a, spectrum_b)
    for i in merged_spectrum.data.index:
        sum += math.sqrt(merged_spectrum.data.loc[i, "Intensity_a"] * merged_spectrum.data.loc[i, "Intensity_b"])
    sum = np.clip(sum, -1, 1) # makes sure the arccos doesn't run into its limit (happens due to approx. error)
    return (np.arccos(sum) ** 2)/(np.arccos(0) ** 2)


def comparison_ms(analyte_spectra, dtb_spectrum_ms, settings):
    """
    Iterates spectra comparison for all databank entries for the MS spectra.
    :return: Annotates any spectra match in the analytes' spectrum info dictionary.
    """
    dot_product_comp = False
    entropy_comp = False
    weighted_dp_comp = False
    bhat1_comp = False

    if not analyte_spectra.info["Normalized"]:
        analyte_spectra.normalize()

    """
    Next section compares spectra and prints matches.
    """

    if settings["ms_algorithm"] == "all":
        if dot_product_similarity(analyte_spectra, dtb_spectrum_ms) < settings["threshold_ms_dot_product_sim"]:
            print("Dot product algorithm shows high MS Spectra similarity: " +\
                  str(dot_product_similarity(analyte_spectra, dtb_spectrum_ms)) +\
                  " with " + dtb_spectrum_ms.info["Molecule"])
            dot_product_comp = True
        if entropy_similarity(analyte_spectra, dtb_spectrum_ms) < settings["threshold_ms_entropy_sim"]:
            print("Spectral Entropy algorithm shows high MS Spectra similarity: " + \
                  str(entropy_similarity(analyte_spectra, dtb_spectrum_ms)) +\
                  " with " + dtb_spectrum_ms.info["Molecule"])
            entropy_comp = True
        if weighted_dp_similarity(analyte_spectra, dtb_spectrum_ms) < settings["threshold_ms_weighted_dp_sim"]:
            print("Weighted dot product algorithm shows high MS Spectra similarity: " + \
                  str(weighted_dp_similarity(analyte_spectra, dtb_spectrum_ms)) + \
                  " with " + dtb_spectrum_ms.info["Molecule"])
            weighted_dp_comp = True
        if bhattacharya1_distance(analyte_spectra, dtb_spectrum_ms) < settings["threshold_ms_bhat1"]:
            print("Bhattacharya 1 algorithm shows high MS Spectra similarity: " + \
                  str(bhattacharya1_distance(analyte_spectra, dtb_spectrum_ms)) + \
                  " with " + dtb_spectrum_ms.info["Molecule"])
            bhat1_comp = True

        if dot_product_comp and entropy_comp and weighted_dp_comp and bhat1_comp:
            return True
        else:
            return False

    elif settings["ms_algorithm"] == "dot_product":
        if dot_product_similarity(analyte_spectra, dtb_spectrum_ms) < settings["threshold_ms_dot_product_sim"]:
            print("Dot product algorithm shows high MS Spectra similarity: " +\
                  str(dot_product_similarity(analyte_spectra, dtb_spectrum_ms)) +\
                  " with " + dtb_spectrum_ms.info["Molecule"])
            dot_product_comp = True
        return dot_product_comp

    elif settings["ms_algorithm"] == "weighted_dot_product":
        if weighted_dp_similarity(analyte_spectra, dtb_spectrum_ms) < settings["threshold_ms_weighted_dp_sim"]:
            print("Weighted dot product algorithm shows high MS Spectra similarity: " + \
                  str(weighted_dp_similarity(analyte_spectra, dtb_spectrum_ms)) + \
                  " with " + dtb_spectrum_ms.info["Molecule"])
            weighted_dp_comp = True
        return weighted_dp_comp

    elif settings["ms_algorithm"] == "bhat1":
        if bhattacharya1_distance(analyte_spectra, dtb_spectrum_ms) < settings["threshold_ms_bhat1"]:
            print("Bhattacharya 1 algorithm shows high MS Spectra similarity: " + \
                  str(bhattacharya1_distance(analyte_spectra, dtb_spectrum_ms)) + \
                  " with " + dtb_spectrum_ms.info["Molecule"])
            bhat1_comp = True
        return bhat1_comp

    elif settings["ms_algorithm"] == "entropy_sim":
        if entropy_similarity(analyte_spectra, dtb_spectrum_ms) < settings["threshold_ms_entropy_sim"]:
            print("Spectral Entropy algorithm shows high MS Spectra similarity: " + \
                  str(entropy_similarity(analyte_spectra, dtb_spectrum_ms)) +\
                  " with " + dtb_spectrum_ms.info["Molecule"])
            entropy_comp = True
        return entropy_comp

    else:
        print(settings["ms_algorithm"] + " is an invalid MS comparison algorithm.")
        return False