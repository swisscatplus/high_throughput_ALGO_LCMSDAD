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


def superimposed_peak_deconvolution(time, data_sum, masses, peak_left, peak_right, plot_ = False):
    threshold = .5
    # tensor = tl.tensor(intensity)
    range_comp = range(2, (len(masses)+1))

    intensities = pd.DataFrame()

    for mass in masses:
        for sec_mass in data_sum.index:  # Not sure why but necessary, avoid key_error
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
    rec_errors = []
    sil_scores = []
    H_matrixes = []
    for n_comp in range_comp:
        W, H = run_nnmf(intensities, n_comp)
        rc_error, sil_score = calculate_erros(W, H, intensities)
        rec_errors.append(rc_error)
        sil_scores.append(sil_score)

        W_df = pd.DataFrame(W, columns=[f"Component {i + 1}" for i in range(n_comp)])
        H_df = pd.DataFrame(H, columns=intensities.columns)
        H_matrixes.append(H_df)

        if plot_:
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
    if plot_:
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

    # For now use min rec error as best n_comp (Check best metric with real superimp data).
    # Could also use Silhouette score or else.
    min_rec_error = min(rec_errors)
    min_rec_error_index = rec_errors.index(min_rec_error)
    h_df = H_matrixes[min_rec_error_index]

    mass_groups = group_mass_values(h_df, threshold)
    nnmf_peaks = fit_grouped_masses(time, mass_groups, data_sum, peak_left, peak_right)
    return nnmf_peaks

def fit_grouped_masses(time, mass_groups, data_sum, peak_left, peak_right, plot_=False):
    total_peak_intensity = np.zeros(shape=len(data_sum.columns))
    nnmf_peaks = []
    for masses in mass_groups:
        for mass in masses:
            total_peak_intensity[peak_left:peak_right] += data_sum.loc[float(mass)].iloc[peak_left:peak_right]
        limited_total_peak_intensity = total_peak_intensity[peak_left:peak_right]
        try:
            successful_fit, fit_parameters, peak_borders, r_squared, peak_integral = fit_custom_peak_fct(masses, limited_total_peak_intensity, time, plot=plot_)
            height = fit_parameters[0] * np.max(limited_total_peak_intensity)
        except RuntimeError:
            print("Peak (NNMF): " + str(masses) + " did not converge.")
            successful_fit = False
        if successful_fit:
            peak = NNMFPeak(fit_parameters, peak_borders, r_squared, peak_integral, height, masses)
            nnmf_peaks.append(peak)
    return nnmf_peaks


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
        group = set([mass])
        candidate_masses = set(list(H_df.columns)) - {mass}
        while candidate_masses:
            cand_mass = candidate_masses.pop()
            if all(H_df.corr().loc[cur_mass, cand_mass] > threshold for cur_mass in group):
                group.add(cand_mass)
            candidate_masses = candidate_masses - {cand_mass}
        mass_groups.append(group)
    mass_groups = set(frozenset(group) for group in mass_groups)  # To remove duplicate groups
    mass_groups = [set(group) for group in mass_groups]  # To convert back to list
    return mass_groups

def fit_custom_peak_fct(name, intensity, time, plot=False, print_=True):
    """
    Fit the custom fct, defined as a skewed gaussian.
    :return:
    """
    window = np.ones(5)/5
    avg_intensity = np.convolve(intensity, window, "same")
    smoothed_intensity = savgol_filter(avg_intensity, window_length=4, polyorder=3)
    maximum_intensity = np.max(smoothed_intensity)
    smoothed_intensity = smoothed_intensity/maximum_intensity

    # Define initial parameters and constrains
    amplitude_guess = 1  # np.max(smoothed_intensity)
    mean_guess = np.argmax(smoothed_intensity)
    std_guess = 10
    skewness_guess = .1
    initial_conditions = [amplitude_guess, mean_guess, std_guess, skewness_guess]
    bounds = ([0, 0, 0, 0], [np.inf, time[-1], 100, 15])

    # Skewed gaussian fit
    parameters, covariance = curve_fit(model_peak, time, smoothed_intensity, p0=initial_conditions, bounds=bounds)
    residual = smoothed_intensity - model_peak(time, *parameters)
    ss_red = np.sum(residual**2)
    ss_tot = np.sum((smoothed_intensity - np.mean(smoothed_intensity))**2)
    r_squared = 1- (ss_red/ss_tot)
    if print_:
        print("Name: " + str(name) + " R2: " + str(r_squared))
        print(*parameters)

    if plot:
        plt.figure()
        plt.scatter(time, smoothed_intensity, label = "Peak data")
        plt.plot(time, model_peak(time, *parameters), label = "Fit")
        plt.xlabel("Scan Number")
        plt.ylabel("Normalized Intensity")
        plt.legend()
        plt.show()

    relative_intensity = 0.25  # Cutoff left and right for the model peak
    def solving_model_borders(x):  # Function in fct bcs we need to access parameters and can't input them due to root_scalar
        """
        Function to determine the left and right border of the peak by their relative height.
        :return:
        """
        target_value = relative_intensity * 1  # 1 is amplitude since we normed the peak to this
        return model_peak(x, *parameters) - target_value

    if r_squared > .90:
        solution_left = root_scalar(solving_model_borders, bracket=[parameters[1] - 10*parameters[2], parameters[1]], method="brentq")
        solution_right = root_scalar(solving_model_borders, bracket=[parameters[1], parameters[1] + 10*parameters[2]], method="brentq")
        peak_borders = (solution_left.root, solution_right.root)
        peak_integral, integral_error = quad(model_peak, -np.inf, np.inf\
                             , args=(parameters[0], parameters[1], parameters[2], parameters[3]))
        peak_integral = peak_integral * maximum_intensity
        return True, parameters, peak_borders, r_squared, peak_integral
    else:
        return False, parameters, None, r_squared, None

def model_peak(x, amplitude, mean, stddev, alpha):
    """
    For now use skewed gaussian
    :param x:
    :param amplitude:
    :param mean:
    :param stddev:
    :param alpha: Parameter to describe skewness
    :return:
    """
    return amplitude * skewnorm.pdf(x, alpha, loc=mean, scale=stddev)

def model_peak_lorentz(x, amplitude, x0, gamma, alpha):  # Currently not used
    """
    Skewed lorentzian as model peak.
    :param x:
    :param amplitude:
    :param x0:
    :param gamma:
    :param alpha: Skewness, >0 -> goes to the right
    :return:
    """
    lorentzian = amplitude * (gamma**2/(x-x0)**2 + gamma**2)
    skew_factor = 1+ alpha*(x-x0)
    return lorentzian * skew_factor

class NNMFPeak:
    def __init__(self, fit_parameters, peak_borders, r_squared, peak_integral, height, masses):
        self.fit_parameters = fit_parameters
        self.peak_borders = peak_borders
        self.r_squared = r_squared
        self.peak_integral = peak_integral
        self.height = height
        self.masses = [float(mass) for mass in masses]

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))