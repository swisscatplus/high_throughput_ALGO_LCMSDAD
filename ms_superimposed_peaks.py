import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
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
from sklearn.metrics import silhouette_score


def superimposed_peak_deconvolution(intensity, time, data_sum, masses, peak_left, peak_right):
    threshold = .5
    # tensor = tl.tensor(intensity)
    range_comp = range(2, (len(masses)+1))

    intensities = pd.DataFrame()

    for mass in masses:
        for sec_mass in data_sum.index:
            if mass == sec_mass:
                new_intensity = data_sum.loc[mass].iloc[peak_left:peak_right]
                # intensities = pd.concat([intensities, new_intensity], ignore_index=True)
                intensities[str(mass)] = new_intensity
    tensor = tl.tensor(intensities)

    """for n_comp in range_comp:
        print(n_comp)
        factors = non_neg_matrix_factorization(tensor, n_comp)
        for x in factors:
            print(x)
        error = reconstruction_error(tensor, factors)
        print(error)"""
    print("----")
    print("len time: " + str(len(intensity)))
    rec_errors = []
    sil_scores = []
    for n_comp in range_comp:
        print(n_comp)
        W, H = run_nnmf(intensities, n_comp)
        rc_error, sil_score = calculate_erros(W, H, intensities)
        rec_errors.append(rc_error)
        sil_scores.append(sil_score)

        W_df = pd.DataFrame(W, columns=[f"Component {i + 1}" for i in range(n_comp)])
        H_df = pd.DataFrame(H, columns=intensities.columns)
        plt.figure(figsize=(10, 4))
        for i in range(n_comp):
            plt.bar(np.arange(len(W_df[f'Component {i+1}'])), W_df[f'Component {i+1}'], label=f"Component {i+1}")
        plt.legend()
        plt.title("Sample Contributions to Components")
        plt.show()
        plt.figure(figsize=(6, 4))
        sns.heatmap(H_df.corr(), annot=True, cmap="coolwarm")
        plt.title("Component Correlations")
        plt.show()
        print(group_mass_values(H_df, threshold))  # later decide which n_comp is best and apply on that
    plt.figure(figsize=(12, 5))

    # Reconstruction error plot
    plt.subplot(1, 2, 1)
    plt.plot(range_comp, rec_errors, marker='o')
    plt.xlabel("Number of Components")
    plt.ylabel("Reconstruction Error")
    plt.title("NNMF Reconstruction Error")

    # Silhouette score plot
    plt.subplot(1, 2, 2)
    plt.plot(range_comp, sil_scores, marker='o')
    plt.xlabel("Number of Components")
    plt.ylabel("Silhouette Score")
    plt.title("NNMF Silhouette Score")

    plt.tight_layout()
    plt.show()
    return

def non_neg_matrix_factorization(tensor, n_comp):
    """
    Calculate factors based on tensorly
    """
    return non_negative_parafac(tensor, rank=n_comp)


def run_nnmf(intensities, n_comp):
    """
    Sklearn-based NMF decomposition
    """
    nmf = NMF(n_components=n_comp)
    W = nmf.fit_transform(intensities)
    H = nmf.components_
    return W, H

def calculate_erros(W, H, intensities):
    reconstruction = np.dot(W, H)
    reconstruction_error = np.linalg.norm(intensities - reconstruction)

    silhouette_error_score = silhouette_score(W, np.argmax(W, axis=1))
    return reconstruction_error, silhouette_error_score

def group_mass_values(H_df, threshold):
    mass_groups = []

    for mass in H_df.columns:
        group = set(mass)
        candidate_masses = set(H_df.columns)
        while candidate_masses:
            cand_mass = candidate_masses.pop()
            if all(H_df.corr().loc[cur_mass, cand_mass] > threshold for cur_mass in group):
                group.add(cand_mass)
        mass_groups.append(group)
    return mass_groups