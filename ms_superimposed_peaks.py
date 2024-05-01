import numpy as np
from matplotlib import pyplot as plt
import hdbscan
import pandas as pd
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt
from scipy.optimize import curve_fit, root_scalar
from scipy.stats import skewnorm
from scipy.integrate import quad
import ms_spectra_comparison
import pywt
import tensorly as tl
from tensorly.decomposition import non_negative_parafac
from sklearn.decomposition import NMF


def superimposed_peak_deconvolution(intensity, time, data_sum, masses, peak_left, peak_right):
    # tensor = tl.tensor(intensity)
    range_comp = range(1, (len(masses)+1))

    intensities = pd.DataFrame()

    for mass in masses:
        for sec_mass in data_sum.index:
            if mass == sec_mass:
                new_intensity = data_sum.loc[mass].iloc[peak_left:peak_right]
                # intensities = pd.concat([intensities, new_intensity], ignore_index=True)
                intensities[str(mass)] = new_intensity
    tensor = tl.tensor(intensities)

    for n_comp in range_comp:
        print(n_comp)
        factors = non_neg_matrix_factorization(tensor, n_comp)
        for x in factors:
            print(x)
        error = reconstruction_error(tensor, factors)
        print(error)
    for n_comp in range_comp:
        print(n_comp)
        nmf = NMF(n_components=n_comp)
        W = nmf.fit_transform(intensities)
        H = nmf.components_
        print(W)
        print(H)
    return

def non_neg_matrix_factorization(tensor, n_comp):
    """
    Calculate factors based on tensorly
    """
    return non_negative_parafac(tensor, rank=n_comp)

def reconstruction_error(tensor, factors):
    """
    Calculate error to evaluate best deconstruction
    """
    reconstruction = tl.kruskal_to_tensor(factors)
    return tl.norm(tensor - reconstruction)