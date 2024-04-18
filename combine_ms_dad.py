import data_processing_dad as dpr_dad
import data_processing as dpr
import ms_deconvolution as msd
import numpy as np
from collections import Counter

"""
Deal with MS and DAD peak picking and alignment of the corresponding signals.
"""

def assign_peaks(full_analysis, method_name, settings):
    """
    Main function, to controll handling of the others.
    :param full_analysis:
    :param settings:
    :return:
    """
    ms_peaks, dad_peaks = create_individual_peaks(full_analysis, settings)
    print(ms_peaks)
    print(dad_peaks)

    times_ms = np.array([peak.max for peak in ms_peaks])
    times_dad = np.array([peak.max for peak in dad_peaks])
    offset = estimate_offset(times_ms, times_dad)
    print(offset)
    ms_peaks = correct_times_offset(ms_peaks, offset)
    times_ms = times_ms - offset

    integrals_ms = np.array([peak.integral for peak in ms_peaks])
    integrals_dad = np.array([peak.integral for peak in dad_peaks])

    matches = match_times(times_ms, times_dad, integrals_ms, integrals_dad)
    print(matches)

    # Now create signal files from the matches peak pairs
    signal_nr = 1
    for match in matches:
        create_signal_from_match(full_analysis, dad_peaks[match[0]], ms_peaks[match[1]], signal_nr, method_name, settings)
        signal_nr += 1

    return

def create_signal_from_match(full_analysis, peak_dad, peak_ms, signal_nr, method_name, settings):
    """
    Function to extract information from the peaks and run through dpr fct which processes the individual spectrum and creates the signal file.
    :param peak_dad:
    :param peak_ms:
    :param settings:
    :return:
    """
    if full_analysis.info["plus_minus_acq"]:
        ms_chr = dpr.MS_full_chr(full_analysis.ms_data3d_polarized, full_analysis.info)
    else:
        ms_chr = dpr.MS_full_chr(full_analysis.ms_data3d, full_analysis.info)

    peak_ms.data = ms_chr.extract_ms_timespan(peak_ms.left, peak_ms.right)  # For now extract full spectrum
    # Later change so that only selected m/z values are put into the spectrum.

    # Rewrite information that is important to signal file
    peak_ms.data.info["Signal Nr"] = str(signal_nr)
    peak_ms.data.info["Run Nr"] = str(full_analysis.info["Run Nr"])
    peak_ms.data.info["Time"] = peak_ms.max
    peak_ms.data.info["LC Method"] = method_name
    peak_dad.data.info["Run Name"] = full_analysis.info["Name"]
    peak_dad.data.info["Pure"] = peak_dad.pure
    peak_dad.data.info["Relative Purity"] = peak_dad.relative_purity
    peak_dad.data.info["Integral"] = peak_dad.integral

    # Make sure that the .info part contains the correct information. It goes into the signal file from there.
    dpr.process_deconvoluted_peak(peak_ms.data, peak_dad.data, full_analysis.info["LC Method"], settings)
    return

def create_individual_peaks(full_analysis, settings):
    """
    Takes in analysis run and puts out list of ms peaks and dad peaks
    :param full_analysis:
    :param settings:
    :return:
    """
    ms_chr = dpr.MS_full_chr(full_analysis.ms_data3d, full_analysis.info)
    return msd.ms_create_peaks(full_analysis, ms_chr, settings), dpr_dad.run_mocca(full_analysis, settings)

def estimate_offset(times_ms, times_dad):
    """
    Function to estimate the offset of the ms peaks towards the dad peaks based on the cross-correlation matrix
    between the time values.
    The offset is constraint to realistic, experimentally observed values.
    :param times_ms:
    :param times_dad:
    :return: Offset value
    """
    # Preset offset max and minimum:
    min_offset = 2
    max_offset = 20

    time_differences = np.abs(times_ms[:, None] - times_dad)
    flat_differences = time_differences.flatten()
    filtered_differences = flat_differences[(flat_differences >= min_offset) & (flat_differences <= max_offset)]
    if len(filtered_differences) > 0:
        offset_estimate = np.median(filtered_differences)  # Median to avoid strong outlier influence (alias peaks that only occur in one spectrum)
    else:
        offset_estimate = max_offset
    return offset_estimate

def correct_times_offset(peak_list, offset):
    """
    Correct the peak maxima accord to the offset. Left and right are not corrected.
    :param peak_list:
    :param offset:
    :return:
    """
    for peak in peak_list:
        peak.max = peak.max - offset  # Only max value
        # left and right need to stay from original, since its needed for spectrum extraction
    return peak_list

def match_times(times_ms, times_dad, integrals_ms, integrals_dad):
    """
    Fct to match dad and ms times. Return indexes of both lists corresponding to the peaks.
    :param times_ms:
    :param times_dad:
    :param offset:
    :return:
    """
    # Preset time tolerance of 8 secounds
    tolerance = 8.

    matches = []

    for i, time_dad in enumerate(times_dad):
        time_differences = np.abs(times_ms - time_dad)
        close_indices_ms = np.where(time_differences <= tolerance)[0]

        if len(close_indices_ms) == 0:  # If nothing in tolerance there are no matches
            # print("No fitting MS")
            continue

        close_indices_dad = find_close_indices(indices_known=close_indices_ms, times_indices_known=times_ms,
                                               times_to_check=times_dad, tolerance=tolerance)

        if len(close_indices_ms) == 1 and len(close_indices_dad) == 1:  # Only 1 DAD and MS peak in tolerance window
            matches.append((i ,close_indices_ms[0]))
            # print("1 and 1")
        elif len(close_indices_ms) > 1 and len(close_indices_dad) == 1:  # Several MS peaks for 1 DAD peak
            relevant_integrals = []
            for index in close_indices_ms:
                relevant_integrals.append(integrals_ms[index])
            largest_integral_index = relevant_integrals.index((max(relevant_integrals)))  # Largest ingtegral should be correct MS peak
            matches.append((i, close_indices_ms[largest_integral_index]))
            # print(" 1 DAD and more MS")
        elif len(close_indices_ms) == 1 and len(close_indices_dad) > 1:  # Several DAD peaks for 1 MS peak
            relevant_integrals = []
            for index in close_indices_dad:
                relevant_integrals.append(integrals_dad[index])
            largest_integral_index = relevant_integrals.index((max(relevant_integrals)))  #  Largest integral as correct DAD peak if only one present
            matches.append((close_indices_dad[largest_integral_index], close_indices_ms[0]))
            # print("1 MS and more DAD")
        elif len(close_indices_ms) > 1 and len(close_indices_dad) > 1:  # Several DAD and MS peaks close by

            # First, iterate over present dad and ms indices and add new indices until no more are being added
            # This is to completely capture "chained" peaks
            new_ms_values = True
            new_dad_values = True
            iteration_counter = 0
            while new_ms_values or new_dad_values:
                new_close_indices_ms = find_close_indices(close_indices_dad, times_dad, times_ms, tolerance)
                if Counter(new_close_indices_ms) == Counter(close_indices_ms):
                    new_ms_values = False
                else:
                    close_indices_ms = new_close_indices_ms
                new_close_indices_dad = find_close_indices(close_indices_ms, times_ms, times_dad, tolerance)
                if Counter(new_close_indices_dad) == Counter(close_indices_dad):
                    new_dad_values = False
                else:
                    close_indices_dad = new_close_indices_dad
                iteration_counter += 1
                if iteration_counter > 1000:  # More than 1000 chained peaks is unrealistic
                    break  # Just in case sth weird happens this acts as an emergency break

            # Now, we calculate the corresponding integral values and order time- and integral- values
            relevant_integrals_ms = []
            relevant_integrals_dad = []
            relevant_times_ms = []
            relevant_times_dad = []
            for index in close_indices_dad:
                relevant_integrals_dad.append(integrals_dad[index])
                relevant_times_dad.append(times_dad[index])
            for index in close_indices_ms:
                relevant_integrals_ms.append(integrals_ms[index])
                relevant_times_ms.append(times_ms[index])
            pairs_ordered_dad = sorted(zip(relevant_times_dad, relevant_integrals_dad))
            pairs_ordered_ms = sorted(zip(relevant_times_ms, relevant_integrals_ms))
            times_ordered_dad, integrals_ordered_dad = zip(*pairs_ordered_dad)
            times_ordered_ms, integrals_ordered_ms = zip(*pairs_ordered_ms)
            total_integral_dad = sum(integrals_ordered_dad)
            total_integral_ms = sum(integrals_ordered_ms)
            rela_integrals_dad = integrals_ordered_dad/total_integral_dad
            rela_integrals_ms = integrals_ordered_ms/total_integral_ms

            min_index_ms = -1
            for index in range(len(times_ordered_dad)):
                rel_integral_dad = rela_integrals_dad[index]
                # Find index of closest relative integral in MS:
                closest_index_ms = min(enumerate(rela_integrals_ms), key = lambda x: abs(x[1]-rel_integral_dad))[0]

                # Find index of closest relative integral in DAD for the MS index found before
                rel_integral_ms = rela_integrals_ms[closest_index_ms]
                closest_index_dad = min(enumerate(rela_integrals_dad), key = lambda  x: abs(x[1]-rel_integral_ms))[0]

                # min_index_ms ensures that peak order stays correct.
                if index == closest_index_dad and closest_index_ms > min_index_ms:
                    # Find the corresponding index values in the full peak lists
                    true_index_dad = np.where(times_dad == times_ordered_dad[index])[0]
                    true_index_ms = np.where(times_ms == times_ordered_ms[closest_index_ms])[0]
                    matches.append((true_index_dad[0], true_index_ms[0]))
                    min_index_ms = closest_index_ms

        else:
            print("Warning: Scenario with neither MS nor DAD Peak shouldn't occur.")
    matches = list(dict.fromkeys(matches))  # Removes double entries
    return matches

def find_close_indices(indices_known, times_indices_known, times_to_check, tolerance):
    """
    Fct to find dad peaks within tolerance.
    :return:
    """
    close_indices_new = np.array(None)
    for i_known in indices_known:
        time_differences_new = np.abs(times_to_check - times_indices_known[i_known])
        indices_in_tolerance = np.where(time_differences_new <= tolerance)[0]
        if close_indices_new == None:
            close_indices_new = [indices_in_tolerance]
        else:
            close_indices_new += [indices_in_tolerance]
    close_indices_new = [tuple(element) if isinstance(element, np.ndarray) else element for element in
                         close_indices_new]  # Extract indices from arrays

    close_indices_new_flatten = []
    for list_ in close_indices_new:
        for index in list_:
            close_indices_new_flatten.append(index)
    close_indices_new = close_indices_new_flatten

    close_indices_new = list(
        dict.fromkeys(close_indices_new))  # Removes doubles values (due to several ms in one dad peak tolerance)
    return close_indices_new
