import os
import csv
import json
import data_processing_dad as dpr_dad
import data_processing as dpr
import ms_deconvolution as msd
import signal_peak_handling as sp_handling
from databank_handling import load_dtb_entry
import ms_spectra_comparison as ms_comp
import dad_spectra_comparison as dad_comp

def comparison_dtb_named_files(settings):

    analyte_folders = ["peaks_dtb_comp"]
    comparison_types_ms = ["dot_product", "weighted_dot_product", "bhat1", "entropy_sim", "all"]
    comparison_types_dad = ["pearson", "dot_product", "ft", "derivative", "all"]
    weighting_fct_types = ["None", "exponential", "exponential2", "logarithmic", "sin", "sin2"]
    file_headings = ["Comparison_type", "Similarity_true", "Similarity_false", "Diff_distribution",
                        "lowest_diff", "avg_diff"]
    lower_limit_ft_range = range(0, 30, 5)
    upper_limit_ft_range = range(200, 701, 50)

    for folder in analyte_folders:
        csv_filename = "ms_optimization_" + folder + ".csv"
        csv_filepath = os.path.join(settings["directory_project"], "Data_examples", csv_filename)
        with open(csv_filepath, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=file_headings)
            writer.writeheader()
        for weighting_fct in weighting_fct_types:
            settings["ms_weighting_fct"] = weighting_fct
            for comparison_type in comparison_types_ms:
                settings["ms_algorithm"] = comparison_type
                scores = compare_analyte_dtb_ms(folder, settings)
                with open(csv_filepath, "a", newline="") as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=file_headings)
                    writer.writerow({
                        file_headings[0]: comparison_type + "_" + weighting_fct,
                        file_headings[1]: json.dumps(scores[0]),
                        file_headings[2]: json.dumps(scores[1]),
                        file_headings[3]: json.dumps(scores[2]),
                        file_headings[4]: json.dumps(scores[3]),
                        file_headings[5]: json.dumps(scores[4])
                    })
    print("Done with MS")
    for folder in analyte_folders:
        csv_filename = "dad_optimization_" + folder + ".csv"
        csv_filepath = os.path.join(settings["directory_project"], "Data_examples", csv_filename)
        with open(csv_filepath, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=file_headings)
            writer.writeheader()
        for comparison_type in comparison_types_dad:
            settings["dad_algorithm"] = comparison_type
            if comparison_type == "ft":
                for lower_limit in lower_limit_ft_range:
                    for upper_limit in upper_limit_ft_range:
                        settings["FT_low_border"] = lower_limit
                        settings["FT_upper_border"] = upper_limit
                        scores = compare_analyte_dtb_dad(folder, settings)
                        with open(csv_filepath, "a", newline="") as csvfile:
                            writer = csv.DictWriter(csvfile, fieldnames=file_headings)
                            writer.writerow({
                                file_headings[0]: comparison_type + "_low_lim_" + str(lower_limit) + "_up_lim_" + str(upper_limit),
                                file_headings[1]: json.dumps(scores[0]),
                                file_headings[2]: json.dumps(scores[1]),
                                file_headings[3]: json.dumps(scores[2]),
                                file_headings[4]: json.dumps(scores[3]),
                                file_headings[5]: json.dumps(scores[4])
                            })
            elif comparison_type == "derivative":
                for nth in range(1, 5):
                    settings["nth_derivative"] = nth
                    scores = compare_analyte_dtb_dad(folder, settings)
                    with open(csv_filepath, "a", newline="") as csvfile:
                        writer = csv.DictWriter(csvfile, fieldnames=file_headings)
                        writer.writerow({
                            file_headings[0]: comparison_type + "_n_" + str(nth),
                            file_headings[1]: json.dumps(scores[0]),
                            file_headings[2]: json.dumps(scores[1]),
                            file_headings[3]: json.dumps(scores[2]),
                            file_headings[4]: json.dumps(scores[3]),
                            file_headings[5]: json.dumps(scores[4])
                        })
            elif comparison_type == "all":
                for nth in range(1, 5):
                    settings["nth_derivative"] = nth
                    for lower_limit in lower_limit_ft_range:
                        for upper_limit in upper_limit_ft_range:
                            settings["FT_low_border"] = lower_limit
                            settings["FT_upper_border"] = upper_limit
                            scores = compare_analyte_dtb_dad(folder, settings)
                            with open(csv_filepath, "a", newline="") as csvfile:
                                writer = csv.DictWriter(csvfile, fieldnames=file_headings)
                                writer.writerow({
                                    file_headings[0]: comparison_type + "_low_lim_" + str(
                                        lower_limit) + "_up_lim_" + str(upper_limit) + "_n_" + str(nth),
                                    file_headings[1]: json.dumps(scores[0]),
                                    file_headings[2]: json.dumps(scores[1]),
                                    file_headings[3]: json.dumps(scores[2]),
                                    file_headings[4]: json.dumps(scores[3]),
                                    file_headings[5]: json.dumps(scores[4])
                                })
            else:
                scores = compare_analyte_dtb_dad(folder, settings)
                with open(csv_filepath, "a", newline="") as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=file_headings)
                    writer.writerow({
                        file_headings[0]: comparison_type,
                        file_headings[1]: json.dumps(scores[0]),
                        file_headings[2]: json.dumps(scores[1]),
                        file_headings[3]: json.dumps(scores[2]),
                        file_headings[4]: json.dumps(scores[3]),
                        file_headings[5]: json.dumps(scores[4])
                    })
    return

def compare_analyte_dtb_ms(analyte_folder, settings, processed = True):
    """
    Fct to get metadata from opt with certain comparison fct and weighting fct (MS).
    """
    directory_database = os.path.join(settings["directory_project"], "Data_examples", "database")
    database_spectra = os.listdir(directory_database)  # List of all files in the database
    new_peaks_directory = os.path.join(settings["directory_project"], "Data_examples", analyte_folder)
    new_peaks_list = os.listdir(new_peaks_directory)

    false_similarity_scores = []
    true_similarity_score = []
    lowest_false_similarity_score = []  # index and value
    all_diff_similarity_scores = []
    avg_diff_similarity_scores = []

    for peak_name in new_peaks_list:
        analyte_path = os.path.join(new_peaks_directory, peak_name)
        if os.path.isfile(analyte_path) and analyte_path.lower().endswith(".cdf"):

            # analyte_signal = sp_handling.load_signal_file(analyte_path, processed)  # for later when its just reloaded peaks..
            analyte_signal = load_dtb_entry(analyte_path, processed)
            analyte_spectrum_ms = analyte_signal.ms_spectrum

            analyte_spectrum_ms.weighting_fct(settings)
            new_false_similarity_scores = []
            new_diff_similarity_scores = []

            for database_entry in database_spectra:
                single_entry = os.path.join(directory_database, database_entry)
                if os.path.isfile(single_entry) and single_entry.lower().endswith(".cdf"):
                    dtb_spectrum = load_dtb_entry(single_entry, processed)  # extract file
                    dtb_spectrum_ms = dtb_spectrum.ms_spectrum  # Obtains MS spectrum from signal object

                    dtb_spectrum_ms.weighting_fct(settings)

                    if settings["ms_algorithm"] == "dot_product":
                        similarity_score = ms_comp.dot_product_similarity(analyte_spectrum_ms, dtb_spectrum_ms)
                    elif settings["ms_algorithm"] == "weighted_dot_product":
                        similarity_score = ms_comp.weighted_dp_similarity(analyte_spectrum_ms, dtb_spectrum_ms)
                    elif settings["ms_algorithm"] == "bhat1":
                        similarity_score = ms_comp.bhattacharya1_distance(analyte_spectrum_ms, dtb_spectrum_ms)
                    elif settings["ms_algorithm"] == "entropy_sim":
                        similarity_score = ms_comp.entropy_similarity(analyte_spectrum_ms, dtb_spectrum_ms)
                    elif settings["ms_algorithm"] == "all":
                        sum_score = ms_comp.dot_product_similarity(analyte_spectrum_ms, dtb_spectrum_ms)
                        sum_score += ms_comp.weighted_dp_similarity(analyte_spectrum_ms, dtb_spectrum_ms)
                        sum_score += ms_comp.bhattacharya1_distance(analyte_spectrum_ms, dtb_spectrum_ms)
                        sum_score += ms_comp.entropy_similarity(analyte_spectrum_ms, dtb_spectrum_ms)
                        similarity_score = sum_score/4
                    else:
                        raise SyntaxError("Error: Unknown comparison algorithm MS")

                    if database_entry == peak_name:
                        true_similarity_score.append(similarity_score)
                    else:
                        new_false_similarity_scores.append(similarity_score)

                else:  # Prints error in case there is a wrong file
                    print("!!! File " + single_entry + " is not a supported database file type!")
            false_similarity_scores.extend(new_false_similarity_scores)
            for score in new_false_similarity_scores:
                diff_scores = score - true_similarity_score[-1]
                new_diff_similarity_scores.append(diff_scores)
            avg_diff_similarity_scores.append((sum(new_diff_similarity_scores)/len(new_diff_similarity_scores)))
            lowest_diff = min(new_diff_similarity_scores)
            lowest_diff_index = new_diff_similarity_scores.index(lowest_diff)
            lowest_false_similarity_score.append(lowest_diff_index)
            lowest_false_similarity_score.append(lowest_diff)
            all_diff_similarity_scores.extend(new_diff_similarity_scores)
        else:
            print("Error: " + peak_name + " is not a valid file.")
    return true_similarity_score, false_similarity_scores, all_diff_similarity_scores, lowest_false_similarity_score, avg_diff_similarity_scores

def compare_analyte_dtb_dad(analyte_folder, settings, processed = True):
    """
    Fct to get metadata from opt with certain comparison fct and weighting fct (MS).
    """
    directory_database = os.path.join(settings["directory_project"], "Data_examples", "database")
    database_spectra = os.listdir(directory_database)  # List of all files in the database
    new_peaks_directory = os.path.join(settings["directory_project"], "Data_examples", analyte_folder)
    new_peaks_list = os.listdir(new_peaks_directory)

    false_similarity_scores = []
    true_similarity_score = []
    lowest_false_similarity_score = []  # index and value
    all_diff_similarity_scores = []
    avg_diff_similarity_scores = []

    for peak_name in new_peaks_list:
        analyte_path = os.path.join(new_peaks_directory, peak_name)
        if os.path.isfile(analyte_path) and analyte_path.lower().endswith(".cdf"):

            # analyte_signal = sp_handling.load_signal_file(analyte_path, processed)  # for later when its just reloaded peaks..
            analyte_signal = load_dtb_entry(analyte_path, processed)
            analyte_spectrum_dad = analyte_signal.dad_spectrum

            new_false_similarity_scores = []
            new_diff_similarity_scores = []

            for database_entry in database_spectra:
                single_entry = os.path.join(directory_database, database_entry)
                if os.path.isfile(single_entry) and single_entry.lower().endswith(".cdf"):
                    dtb_spectrum = load_dtb_entry(single_entry, processed)  # extract file
                    dtb_spectrum_dad = dtb_spectrum.dad_spectrum  # Obtains MS spectrum from signal object

                    if settings["dad_algorithm"] == "pearson":
                        similarity_score = dad_comp.pearson_similarity(analyte_spectrum_dad, dtb_spectrum_dad)
                    elif settings["dad_algorithm"] == "dot_product":
                        similarity_score = dad_comp.dot_product_similarity(analyte_spectrum_dad, dtb_spectrum_dad)
                    elif settings["dad_algorithm"] == "ft":
                        similarity_score = dad_comp.ft_similarity(analyte_spectrum_dad, dtb_spectrum_dad, settings)
                    elif settings["dad_algorithm"] == "derivative":
                        similarity_score = dad_comp.derivative_similarity(analyte_spectrum_dad, dtb_spectrum_dad, settings)
                    elif settings["ms_algorithm"] == "all":
                        sum_score = dad_comp.pearson_similarity(analyte_spectrum_dad, dtb_spectrum_dad)
                        sum_score += dad_comp.dot_product_similarity(analyte_spectrum_dad, dtb_spectrum_dad)
                        sum_score += dad_comp.ft_similarity(analyte_spectrum_dad, dtb_spectrum_dad, settings)
                        sum_score += dad_comp.derivative_similarity(analyte_spectrum_dad, dtb_spectrum_dad, settings)
                        similarity_score = sum_score/4
                    else:
                        raise SyntaxError("Error: Unknown comparison algorithm MS")

                    if database_entry == peak_name:
                        true_similarity_score.append(similarity_score)
                    else:
                        new_false_similarity_scores.append(similarity_score)

                else:  # Prints error in case there is a wrong file
                    print("!!! File " + single_entry + " is not a supported database file type!")
            false_similarity_scores.extend(new_false_similarity_scores)
            for score in new_false_similarity_scores:
                diff_scores = score - true_similarity_score[-1]
                new_diff_similarity_scores.append(diff_scores)
            avg_diff_similarity_scores.append((sum(new_diff_similarity_scores)/len(new_diff_similarity_scores)))
            lowest_diff = min(new_diff_similarity_scores)
            lowest_diff_index = new_diff_similarity_scores.index(lowest_diff)
            lowest_false_similarity_score.append(lowest_diff_index)
            lowest_false_similarity_score.append(lowest_diff)
            all_diff_similarity_scores.extend(new_diff_similarity_scores)
        else:
            print("Error: " + peak_name + " is not a valid file.")
    return true_similarity_score, false_similarity_scores, all_diff_similarity_scores, lowest_false_similarity_score, avg_diff_similarity_scores