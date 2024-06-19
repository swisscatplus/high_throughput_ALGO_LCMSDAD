import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt
from scipy.optimize import curve_fit, root_scalar
from scipy.stats import skewnorm
from scipy.integrate import quad
from . import ms_spectra_comparison
import pywt
import os
from netCDF4 import Dataset
from . import ms_superimposed_peaks as ms_spr



"""
Functions to pick peaks in MS chromatogram and consequently perform deconvolution
"""

def ms_create_peaks(full_analysis, ms_chr, settings, plot=False):

    background_masses_list = load_background_masses_list(full_analysis.info["LC Method"], settings)
    print(background_masses_list)

    entropy = np.zeros(shape=len(full_analysis.ms_data3d["MS spectra"]))
    ms_intensity = np.zeros(shape=len(full_analysis.ms_data3d["time"]))
    ms_value = np.zeros(shape=len(full_analysis.ms_data3d["time"]))
    if plot:
        columns_names = [str(i) for i in range(0, 5)]
        df_top_masses = pd.DataFrame(index = full_analysis.ms_data3d["time"], columns=columns_names)
        ms_single_values = pd.DataFrame(index = range(len(full_analysis.ms_data3d["time"])))
        for i in np.arange(225, 226, 0.1):
            i = round(i, 1)
            ms_single_values[i] = ms_value
    df_heatmap = pd.DataFrame(index = np.arange(0.1, 1000.0, 0.1), columns = full_analysis.ms_data3d["time"])
    # prefilling the df with zeros is much slower when replacing them with intensity values -> fillna at the end

    index = 0
    window = np.ones(10)/10
    for i in full_analysis.ms_data3d["time"]:
        ms_spectrum = ms_chr.extract_single_ms(i)
        # print(index)

        entropy[index] = ms_spectra_comparison.calculate_entropy(ms_spectrum)

        ms_intensity[index] = sum(ms_spectrum.data["Intensity"])
        sorted_index = np.argsort(-ms_spectrum.data["Intensity"])[:5]

        if plot:
            for value in np.arange(225, 226, 0.1):
                value = round(value, 1)
                for mass in ms_spectrum.data["m/z"]:
                    if mass == value:
                        value_index = ms_spectrum.data.index[ms_spectrum.data["m/z"] == value]
                        ms_single_values[value][index] = ms_spectrum.data["Intensity"][value_index]
        for mass, intensity in zip(ms_spectrum.data["m/z"], ms_spectrum.data["Intensity"]):
            if mass in df_heatmap.index and intensity > 2000:  # 10k Probably too high
                df_heatmap.at[mass, i] = intensity  # Tried set to 1 for clustering, didn't work well...
        index += 1
    df_heatmap.fillna(0, inplace=True)
    filtered_time = full_analysis.ms_data3d["time"]

    # To allow selection of one polarity when creating the spectrum

    if full_analysis.info["plus_minus_acq"]:
        if settings["ion_detection_mode"] == "positive":
            filtered_time = full_analysis.ms_data3d["time"][::2]
            new_entropy = entropy[::2]
            new_ms_intensity = ms_intensity[::2]
            # Now copy polarized into full_analysis to later extract spectra
            ms_data_polarized = np.empty((len(full_analysis.ms_data3d["MS spectra"][:, 0][::2]), 2), dtype=object)
            ms_data_polarized[:, 0] = full_analysis.ms_data3d["MS spectra"][:, 0][::2]
            ms_data_polarized[:, 1] = full_analysis.ms_data3d["MS spectra"][:, 1][::2]
            ms_data3d_polarized = {
                "time": full_analysis.ms_data3d["time"][::2],
                "total intensity": full_analysis.ms_data3d["total intensity"][::2],
                "MS spectra": np.array(ms_data_polarized),
            }
            full_analysis.ms_data3d_polarized = ms_data3d_polarized
        elif settings["ion_detection_mode"] == "negative":
            filtered_time = full_analysis.ms_data3d["time"][1::2]
            new_entropy = entropy[1::2]
            new_ms_intensity = ms_intensity[1::2]
            # Now copy polarized into full_analysis to later extract spectra
            ms_data_polarized = np.empty((len(full_analysis.ms_data3d["MS spectra"][:, 0][1::2]), 2), dtype=object)
            ms_data_polarized[:, 0] = full_analysis.ms_data3d["MS spectra"][:, 0][1::2]
            ms_data_polarized[:, 1] = full_analysis.ms_data3d["MS spectra"][:, 1][1::2]
            ms_data3d_polarized = {
                "time": full_analysis.ms_data3d["time"][1::2],
                "total intensity": full_analysis.ms_data3d["total intensity"][1::2],
                "MS spectra": np.array(ms_data_polarized),
            }
            full_analysis.ms_data3d_polarized = ms_data3d_polarized
        else:
            raise TypeError("Ion detection mode must be either positive or negative.")
    else:
        new_entropy = entropy
        new_ms_intensity = ms_intensity

    times_to_remove = [t for t in df_heatmap.columns if t not in filtered_time]
    reduced_heatmap = df_heatmap.drop(columns = times_to_remove)

    entropy_peaks = ms_entropy_peaks(new_entropy)

    if plot:
        plt.figure()
        plt.pcolormesh(reduced_heatmap, cmap="BuPu")
        plt.xlabel('Time [s]')
        plt.ylabel('m/z')
        plt.colorbar()
        plt.show()

        for i in np.arange(225, 226, 0.1):
            i = round(i, 1)
            plt.plot(full_analysis.ms_data3d["time"][::2], ms_single_values[i][::2], label=i)
        plt.legend()
        plt.show()

        #plt.plot(full_analysis.ms_data3d["time"], ms_intensity)
        plt.plot(filtered_time, new_ms_intensity)
        plt.xlabel("Time [s]")
        plt.ylabel("Total Ion Count")
        plt.show()

        #plt.plot(full_analysis.ms_data3d["time"], entropy)
        plt.plot(filtered_time, new_entropy)
        plt.xlabel("Time [s]")
        plt.ylabel("Entropy")
        plt.title("Spectral Entropy")
        plt.show()

    ms_peak_list = ms_summation(reduced_heatmap, entropy_peaks, background_masses_list, settings)
    return ms_peak_list

def determine_background_masses_list(full_analysis, ms_chr, settings):
    ms_intensity = np.zeros(shape=len(full_analysis.ms_data3d["time"]))

    df_heatmap = pd.DataFrame(index=np.arange(0.1, 1000.0, 0.1), columns=full_analysis.ms_data3d["time"])
    # prefilling the df with zeros is much slower when replacing them with intensity values -> fillna at the end

    index = 0
    for i in full_analysis.ms_data3d["time"]:
        ms_spectrum = ms_chr.extract_single_ms(i)
        # print(index)

        ms_intensity[index] = sum(ms_spectrum.data["Intensity"])

        for mass, intensity in zip(ms_spectrum.data["m/z"], ms_spectrum.data["Intensity"]):
            if mass in df_heatmap.index and intensity > 2000:
                df_heatmap.at[mass, i] = intensity
        index += 1
    df_heatmap.fillna(0, inplace=True)
    filtered_time = full_analysis.ms_data3d["time"]

    # To allow selection of one polarity when creating the spectrum

    if full_analysis.info["plus_minus_acq"]:
        if settings["ion_detection_mode"] == "positive":
            filtered_time = full_analysis.ms_data3d["time"][::2]
            new_ms_intensity = ms_intensity[::2]
            # Now copy polarized into full_analysis to later extract spectra
            ms_data_polarized = np.empty((len(full_analysis.ms_data3d["MS spectra"][:, 0][::2]), 2), dtype=object)
            ms_data_polarized[:, 0] = full_analysis.ms_data3d["MS spectra"][:, 0][::2]
            ms_data_polarized[:, 1] = full_analysis.ms_data3d["MS spectra"][:, 1][::2]
            ms_data3d_polarized = {
                "time": full_analysis.ms_data3d["time"][::2],
                "total intensity": full_analysis.ms_data3d["total intensity"][::2],
                "MS spectra": np.array(ms_data_polarized),
            }
            full_analysis.ms_data3d_polarized = ms_data3d_polarized
        elif settings["ion_detection_mode"] == "negative":
            filtered_time = full_analysis.ms_data3d["time"][1::2]
            new_ms_intensity = ms_intensity[1::2]
            # Now copy polarized into full_analysis to later extract spectra
            ms_data_polarized = np.empty((len(full_analysis.ms_data3d["MS spectra"][:, 0][1::2]), 2), dtype=object)
            ms_data_polarized[:, 0] = full_analysis.ms_data3d["MS spectra"][:, 0][1::2]
            ms_data_polarized[:, 1] = full_analysis.ms_data3d["MS spectra"][:, 1][1::2]
            ms_data3d_polarized = {
                "time": full_analysis.ms_data3d["time"][1::2],
                "total intensity": full_analysis.ms_data3d["total intensity"][1::2],
                "MS spectra": np.array(ms_data_polarized),
            }
            full_analysis.ms_data3d_polarized = ms_data3d_polarized
        else:
            raise TypeError("Ion detection mode must be either positive or negative.")
    else:
        new_ms_intensity = ms_intensity

    times_to_remove = [t for t in df_heatmap.columns if t not in filtered_time]
    reduced_heatmap = df_heatmap.drop(columns=times_to_remove)

    sum_labels = reduced_heatmap.index * 10 // 5
    data_sum = reduced_heatmap.groupby(sum_labels).sum()
    new_indices = data_sum.index * 0.5 + 0.25  # addition to correct that m/z gives the avg of the m/z kernel
    data_sum.set_index(new_indices, inplace=True)

    if settings["ion_detection_mode"] == "negative":
        background_threshold = 2000000
    else:
        background_threshold = 5000000

    list_background_masses = []
    for i in data_sum.index:
        intensity = np.array(data_sum.loc[i])
        total_intensity = np.sum(intensity)
        if total_intensity > background_threshold:
            list_background_masses.append(i)
    return list_background_masses

def ms_entropy_peaks(filtered_entropy, plot = False):
    """
    Take in filtered entropy entries, return list of entropy-based peaks (local minima)
    :return:
    """
    window = np.ones(10)/10

    # For local minima
    baseline = -als_baseline(-filtered_entropy)  # Negative sign bcs we want to emphatize minima
    corrected_entropy = -(filtered_entropy - baseline)  # Not sure if this is good, creates artefacts at beginning and end!
    averaged_entropy = np.convolve(corrected_entropy, window, "same")

    # Wavelet denoising to remove high-frequency noise
    wvt_coef_entropy = pywt.wavedec(averaged_entropy, "db1", mode="symmetric")
    threshold_wvt = .05
    wvt_coef_ent_thr = [pywt.threshold(c, value=threshold_wvt, mode="soft") for c in wvt_coef_entropy]

    denoised_entropy = pywt.waverec(wvt_coef_ent_thr, "db1", mode="symmetric")
    smoothed_entropy = savgol_filter(denoised_entropy, window_length=4, polyorder=2)

    peaks = find_peaks_cwt(smoothed_entropy, 5, min_snr=1.)  # No negative sign bcs entropy was inversed beforehand
    filtered_peaks = [peak for peak in peaks if smoothed_entropy[peak] > .15]  # Soft criterium, not to loose highly fragmenting compounds

    # print("Entropy-Peaks: " + str(filtered_peaks))
    if plot:
        plt.plot(range(len(smoothed_entropy)), smoothed_entropy)
        plt.plot(range(len(averaged_entropy)), averaged_entropy)
        plt.plot(range(len(baseline)), baseline)
        plt.plot(range(len(filtered_entropy)), filtered_entropy)

        plt.show()

    # For local maxima
    baseline2 = als_baseline(filtered_entropy)
    corrected_entropy2 = (filtered_entropy - baseline2)
    averaged_entropy2 = np.convolve(corrected_entropy2, window, "same")


    wvt_coef_entropy2 = pywt.wavedec(averaged_entropy2, "db1", mode="symmetric")
    threshold_wvt = .05
    wvt_coef_ent_thr2 = [pywt.threshold(c, value=threshold_wvt, mode="soft") for c in wvt_coef_entropy2]

    denoised_entropy2 = pywt.waverec(wvt_coef_ent_thr2, "db1", mode="symmetric")
    smoothed_entropy2 = savgol_filter(denoised_entropy2, window_length=4, polyorder=2)

    peaks2 = find_peaks_cwt(smoothed_entropy2, 5, min_snr=1.)
    filtered_peaks2 = [peak for peak in peaks2 if smoothed_entropy2[peak] > .15]  # Soft criterium, not to loose highly fragmenting compounds

    return filtered_peaks + filtered_peaks2  # adding local minima and local maxima


def ms_summation(data, entropy_peaks, background_masses_list, settings):
    """
    Fct to determine peaks by summation beforehand.
    :param data:
    :return:
    """
    max_peak_width = 20

    sum_labels = data.index*10 // 5
    data_sum = data.groupby(sum_labels).sum()
    new_indices = data_sum.index *0.5 + 0.25  # addition to correct that m/z gives the avg of the m/z kernel
    data_sum.set_index(new_indices, inplace=True)

    peak_dict = ms_peak_picking(data_sum, background_masses_list, settings)

    peak_clusters, inverse_peaklist = determine_peak_clusters(peak_dict, max_peak_width)
    ms_peak_list = process_ms_peaks(peak_clusters, inverse_peaklist, data_sum, entropy_peaks, max_peak_width)
    return ms_peak_list

def determine_peak_clusters(peak_dict, max_peak_width, print_=False):
    """
    Fct to take in the list of peaks and return the clustered m/z belonging together. And the times corresponding
    to the groups.
    :return:
    """
    inverse_peaklist = defaultdict(list)
    peak_clusters = defaultdict(set)
    mass_to_cluster = defaultdict(set)

    for mass, peaks in peak_dict.items():
        for peak in peaks:
            cluster_ids = find_cluster(peak, peak_clusters, max_peak_width)
            if not cluster_ids:
                cluster_id = len(peak_clusters)
                peak_clusters[cluster_id].add(peak)
                mass_to_cluster[mass].add(cluster_id)
            else:
                for cluster_id in cluster_ids:
                    peak_clusters[cluster_id].add(peak)
                    mass_to_cluster[mass].add(cluster_id)

            if mass not in mass_to_cluster.values():
                mass_to_cluster[mass].add(cluster_id)
    peak_clusters = {m: list(p) for m, p in peak_clusters.items()}

    for mass, cluster_ids in mass_to_cluster.items():
        for cluster_id in cluster_ids:
            inverse_peaklist[cluster_id].append(mass)

    inverse_peaklist = {m: list(set(p)) for m, p in inverse_peaklist.items()}

    # Print analytics
    if print_:
        print(peak_dict)
        print(peak_clusters)
        print(mass_to_cluster)
        print(inverse_peaklist)
    return peak_clusters, inverse_peaklist

def find_cluster(peak, clusters, max_peak_width):
    cluster_ids = []
    for cluster_id, peaks in clusters.items():
        if all(abs(peak - other) <= max_peak_width for other in peaks):
            cluster_ids.append(cluster_id)
    return cluster_ids

def time_in_ent_peaks(times, entropy_peaks, max_peak_width):
    """
    Fct to check if an entropy peak is present at any given time of a custer peak.
    :param times:
    :param entropy_peaks:
    :return: Boolean
    """
    """for time in times:
        if any(abs(time - entropy_peak) <= max_peak_width for entropy_peak in entropy_peaks):
            return True  # Maybe change to more strict criteria (e.g. peak_avg instead of any()
    return False"""
    peak_avg, peak_left, peak_right = peak_average(times)
    for entropy_peak in entropy_peaks:
        if abs(peak_avg - entropy_peak) <= max_peak_width:  # More harsh criteria, kept for now
            # print(peak_avg, entropy_peak)
            return True
    return True  # Inversed for now, don't want entropy criteria CHANGE LATER

def peak_average(peak_list):
    peak_sum = sum(peak_list)
    peak_avg = peak_sum/len(peak_list)
    peak_left = int(peak_avg - 50)
    peak_right = int(peak_avg + 100)  # Asymmetric bcs real peaks tend to tail
    if peak_left < 0:
        peak_left = 0
    return peak_avg, peak_left, peak_right

def als_baseline(intensities, lam=1e5, p=0.01, niter=10):
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    L = len(intensities)
    D = diags([1, -2, 1], [0, -1, -2], shape=(L, L), format='csr')
    w = np.ones(L)
    for i in range(niter):
        W = diags([w], [0], shape=(L, L), format='csr')
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * intensities)
        w = p * (intensities > z) + (1-p) * (intensities < z)
    return z

def process_ms_peaks(peak_clusters, inverse_peaklist, data_sum, entropy_peaks, max_peak_width, plot = False):
    """
    Compare summed peaks with entropy peaks. Fit those who exist in both with gaussian fct.
    :return:
    """

    entropy_confirmed_peaks = []
    ms_peak_list = []
    for name, times in peak_clusters.items():
        if time_in_ent_peaks(times, entropy_peaks, max_peak_width):
            entropy_confirmed_peaks.append(name)
        """for time in times:
            print(data_sum.columns[time])
        print("--------")
    print(entropy_confirmed_peaks)"""
    for name, masses in inverse_peaklist.items():
        if name in entropy_confirmed_peaks:
            total_peak_intensity = np.zeros(shape=len(data_sum.columns))
            peak_avg, peak_left, peak_right = peak_average(peak_clusters[name])
            if peak_right >= len(total_peak_intensity):
                peak_right = len(total_peak_intensity) - 1
            for mass in masses:
                for sec_mass in data_sum.index:  # Not sure why, but necessary to avoid key_error
                    if mass == sec_mass:
                        total_peak_intensity[peak_left:peak_right] += data_sum.loc[mass].iloc[peak_left:peak_right]
            # print(name, masses)
            time_index = np.round(np.array(data_sum.columns), 3)
            """test_time = range(len(total_peak_intensity))
            limited_time_index = time_index[peak_left:peak_right]"""
            limited_total_peak_intensity = total_peak_intensity[peak_left:peak_right]
            time_fits = range(len(limited_total_peak_intensity))
            nnmf_peaks = []
            try:
                successful_fit, fit_parameters, peak_borders, r_squared, peak_integral = ms_spr.fit_custom_peak_fct(name, limited_total_peak_intensity, time_fits)
                if len(masses) > 1 and r_squared > .5 and not successful_fit:
                    nnmf_peaks = ms_spr.superimposed_peak_deconvolution(time_fits, data_sum, masses, peak_left, peak_right)
                    successful_fit = False
            except RuntimeError:
                print("Peak: " + str(name) + " did not converge.")
                successful_fit = False

            if plot: #and successful_fit:
                plt.plot(time_index, total_peak_intensity)
                plt.show()

            if successful_fit:
                height = fit_parameters[0] * np.max(limited_total_peak_intensity)
                mean = time_index[round(fit_parameters[1] + peak_left)]  # testing here
                stwd = fit_parameters[2]
                skewness = fit_parameters[3]

                if round(peak_borders[0] + peak_left) > 0:
                    peak_left_fit = time_index[round(peak_borders[0] + peak_left)]
                else:
                    peak_left_fit = time_index[0]
                if round(peak_borders[1] + peak_left) < time_index.size:
                    peak_right_fit = time_index[round(peak_borders[1] + peak_left)]
                else:
                    peak_right_fit = time_index[-1]

                new_ms_peak = ms_peak(height, mean, stwd, skewness, masses, r_squared,
                                      peak_left_fit, peak_right_fit, peak_integral)
                ms_peak_list.append(new_ms_peak)
            elif len(nnmf_peaks)>0:
                for nnmf_peak in nnmf_peaks:
                    # height needs to be parsed already as different total peak intensity
                    mean = time_index[round(nnmf_peak.fit_parameters[1] + peak_left)]  # testing here
                    stwd = nnmf_peak.fit_parameters[2]
                    skewness = nnmf_peak.fit_parameters[3]
                    if round(nnmf_peak.peak_borders[0] + peak_left) > 0:
                        peak_left_fit = time_index[
                            round(nnmf_peak.peak_borders[0] + peak_left)]  # Replace peak borders by gaussian hitting baseline?
                    else:
                        peak_left_fit = time_index[0]
                    if round(nnmf_peak.peak_borders[1] + peak_left) < time_index.size:
                        peak_right_fit = time_index[round(nnmf_peak.peak_borders[1] + peak_left)]
                    else:
                        peak_right_fit = time_index[-1]
                    new_ms_peak = ms_peak(nnmf_peak.height, mean, stwd, skewness, nnmf_peak.masses, nnmf_peak.r_squared, peak_left_fit,
                                          peak_right_fit, nnmf_peak.peak_integral, nnmf_decon=True)
                    ms_peak_list.append(new_ms_peak)
    return ms_peak_list

def ms_peak_picking(data_sum, background_masses_list, settings, plot = False, print_ = False):
    """
    :return:
    """
    window = np.ones(5) / 5
    peak_dict = {}
    for i in data_sum.index:  # [:] select certain mass values
        intensity = np.array(data_sum.loc[i])
        if i in background_masses_list:
            intensity[:] = 0
        baseline = als_baseline(intensity)
        corrected_intensity = intensity - baseline
        avg_intensity = np.convolve(corrected_intensity, window, "same")
        smoothed_intensity = savgol_filter(avg_intensity, window_length=4, polyorder=3)
        # baseline = savgol_filter(smoothed_intensity, len(smoothed_intensity), 5)  # doesn't effect peak picking
        """baseline = als_baseline(smoothed_intensity)
        corrected_intensity = smoothed_intensity - baseline"""
        peaks = find_peaks_cwt(smoothed_intensity, 10)  # try increasing width?
            # Maybe exchange this for individual background subtraction and checking if peaks exceeds threshold
        if settings["ion_detection_mode"] == "negative":
            peak_threshold = 3000  # 15 000 before
        else:
            peak_threshold = 10000
        filtered_peaks = [peak for peak in peaks if smoothed_intensity[peak] > peak_threshold]  # Maybe change threshold?
        # Changed to smoothed_intensity, so that peaks directly between the window frames, are not ignored due to random 0 at peak.
        # -> readjust thresholds? Could be necessary to lower them. TEST! -> Yes, it is necessary
        index_test = 0

        if print_:
            for x in intensity:
                print("--------")
                print("intensity")
                print(x)
                print(data_sum.columns[index_test])
                print(index_test)
                index_test += 1
            index_test = 0
            for x in smoothed_intensity:
                print("--------")
                print(x)
                print(data_sum.columns[index_test])
                print(index_test)
                index_test +=1
            print("mass: " + str(i))
            print(peaks)
            print(filtered_peaks)
        peak_dict[i] = filtered_peaks
        if plot:
            time_index = np.round(np.array(data_sum.columns), 3)
            for t in peaks:
                print(time_index[t])
            plt.plot(time_index, intensity)
            plt.title("m/z value: " + str(i))
            plt.xlabel("Time [s]")
            plt.ylabel("Counts")
            plt.show()
            plt.plot(time_index, smoothed_intensity)
            plt.title("m/z value: " + str(i))
            plt.xlabel("Time [s]")
            plt.ylabel("Counts")
            plt.show()
    return peak_dict

def load_background_masses_list(method, settings):
    """
    Searches for the background file corresponding to a method, returns list of background masses.
    :param:
    :return:
    """
    background_folder = os.path.join(settings["directory_project"], "Data_examples", "background_spectra")
    background_list = os.listdir(background_folder)
    background_masses_list = None

    for background_file_name in background_list:
        background_file_path = os.path.join(background_folder, background_file_name)
        if os.path.isfile(background_file_path) and background_file_path.lower().endswith(".cdf"):
            if background_file_name[:-4] == method:  # -4 to remove .cdf ending!
                background_masses_list = extract_background_masses_list(background_file_path)
        else:
            print("!!! File " + background_file_name + " is not a supported background file type!")

    if background_masses_list.all() == None:
        print("No backgroundfile exists for the method " + method)
        return
    return background_masses_list

def extract_background_masses_list(background_file_path):
    background_file = Dataset(background_file_path, "r")

    background_masses_list = background_file.groups["MS Data"].variables["Background masses list"][:].compressed(),
    background_file.close()

    return background_masses_list[0]

class ms_peak:
    """
    Class defined to hold the necessary information on a ms_peak
    """
    def __init__(self, height, mean, stdw, skewness, mass_values, r_squared, left, right, integral, nnmf_decon = False):
        self.height = height
        self.max = mean
        self.stdw = stdw
        self.skewness = skewness
        self.mass_values = mass_values
        self.left = left
        self.right = right
        self.integral = integral
        if r_squared > 0.95 and not nnmf_decon:  # Later this needs to be evaluated by a different criteria
            self.pure = True
        else:
            self.pure = False
        self.nnmf = nnmf_decon

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))