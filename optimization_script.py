import os
import csv
import json
import pandas as pd
import numpy as np
from databank_handling import load_dtb_entry
import ms_spectra_comparison as ms_comp
import dad_spectra_comparison as dad_comp
from matplotlib import pyplot as plt
import seaborn as sns


def plot_optimization_dad(csv_file_name, settings):
    directory_project = settings["directory_project"]
    csv_file_path = os.path.join(directory_project, "Data_examples", csv_file_name)

    data = pd.read_csv(csv_file_path)
    for column in data.columns:
        if not column == "Comparison_type":
            data[column] = data[column].apply(json.loads)

    plot_extract_histogram(data, "total", "avg_diff")

    # Plot heatmap average difference
    df_ft_avg = extract_ft_avg(data, "lowest_diff")
    plot_heatmap_all_alg(df_ft_avg, "Fourier Transform Average Difference")
    df_ft_min = extract_ft_min(data, "lowest_diff")
    plot_heatmap_all_alg(df_ft_min, "Fourier Transform Min Average Difference")

    df_others_min = extract_min_column_dad(data, "lowest_diff")
    plot_heatmap_all_alg(df_others_min, "DAD Algorithms Average Difference")

    df_others_avg = extract_avg_column_dad(data, "lowest_diff")
    plot_heatmap_all_alg(df_others_avg, "DAD Algorithms Average Difference")

    return

def plot_optimization_ms(csv_file_name, settings):
    directory_project = settings["directory_project"]
    csv_file_path = os.path.join(directory_project, "Data_examples", csv_file_name)

    data = pd.read_csv(csv_file_path)
    for column in data.columns:
        if not column == "Comparison_type":
            data[column] = data[column].apply(json.loads)

    plot_extract_histogram(data, "dot_product_sin2", "Similarity_false")

    # Plot heatmap average difference
    df_avg_diff = extract_avg_column_ms(data, "avg_diff")
    plot_heatmap_all_alg(df_avg_diff, "Averages Difference")
    df_min_avg_diff = extract_min_column_ms(data, "avg_diff")
    plot_heatmap_all_alg(df_min_avg_diff, "Minimum avg Difference")

    # Plot heatmap lowest difference
    df_low_diff = extract_avg_column_ms(data, "lowest_diff")
    plot_heatmap_all_alg(df_low_diff, "Lowest Difference Average")
    df_min_low_diff = extract_min_column_ms(data, "lowest_diff")
    plot_heatmap_all_alg(df_min_low_diff, "Minimum Lowest Difference")
    return

def plot_extract_histogram(data, algorithm_name, type_name):
    index = 0
    histogram_data = []
    for values in data[type_name]:
        if algorithm_name == "total":  # total for addition of all rows
            histogram_data += [value for value in values if type(value) == float]
        else:
            if data["Comparison_type"][index] == algorithm_name:
                histogram_data = [value for value in values if type(value) == float]
        index += 1
    if len(histogram_data) == 0:
        raise ValueError("Wrong algorithm name or type name")
    sns.histplot(histogram_data, kde=True)
    plt.title(f"Histogram of {type_name} for {algorithm_name}")
    plt.xlabel("Score")
    plt.ylabel("Counts")
    plt.show()
    return

def extract_min_column_dad(data, column_name):
    indices = ["pearson", "dot_product", "ft_low_lim_0_up_lim_200", "derivative_n_1", "derivative_n_2", "derivative_n_3", "derivative_n_4"]
    df = pd.DataFrame({
        "Algorithm": indices,
    })
    df.set_index("Algorithm", inplace=True)
    index = 0
    for i in indices:
        df[i] = np.nan
    indices_index = 0
    for values in data[column_name]:
        values = [value for value in values if type(value) == float]
        if data["Comparison_type"][index].startswith(indices[indices_index]):
            df.loc[indices[indices_index]] = np.min(values)
            indices_index += 1
        if indices_index > 6:
            indices_index = 0
        index += 1
    return df

def extract_avg_column_dad(data, column_name):
    indices = ["pearson", "dot_product", "ft_low_lim_0_up_lim_200", "derivative_n_1", "derivative_n_2", "derivative_n_3", "derivative_n_4"]
    df = pd.DataFrame({
        "Algorithm": indices,
    })
    df.set_index("Algorithm", inplace=True)
    index = 0
    for i in indices:
        df[i] = np.nan
    indices_index = 0
    for values in data[column_name]:
        values = [value for value in values if type(value) == float]
        if data["Comparison_type"][index].startswith(indices[indices_index]):
            df.loc[indices[indices_index]] = np.average(values)
            indices_index += 1
        if indices_index > 6:
            indices_index = 0
        index += 1
    return df

def extract_ft_avg(data, column_name):
    indices = [str(low_lim) for low_lim in range(0, 30, 5)]
    columns = [str(up_lim) for up_lim in range(25, 500, 25)]
    df = pd.DataFrame({
        "Lower Limit": indices,
    })
    df.set_index("Lower Limit", inplace=True)
    index = 0
    for up_lim in columns:
        df[up_lim] = np.nan

    indices_index = 0
    columns_index = 0
    for values in data[column_name]:
        values = [value for value in values if type(value) == float]
        if data["Comparison_type"][index].startswith(f"ft_low_lim_{indices[indices_index]}"):
            if data["Comparison_type"][index].endswith(f"up_lim_{columns[columns_index]}"):
                df[columns[columns_index]].loc[indices[indices_index]] = np.average(values)
                columns_index += 1
        if columns_index > 18:
            columns_index = 0
            indices_index += 1
        if indices_index > 5:
            indices_index = 0
        index += 1
    return df

def extract_ft_min(data, column_name):
    indices = [str(low_lim) for low_lim in range(0, 30, 5)]
    columns = [str(up_lim) for up_lim in range(25, 500, 25)]
    df = pd.DataFrame({
        "Lower Limit": indices,
    })
    df.set_index("Lower Limit", inplace=True)
    index = 0
    for up_lim in columns:
        df[up_lim] = np.nan

    indices_index = 0
    columns_index = 0
    for values in data[column_name]:
        values = [value for value in values if type(value) == float]
        if data["Comparison_type"][index].startswith(f"ft_low_lim_{indices[indices_index]}"):
            if data["Comparison_type"][index].endswith(f"up_lim_{columns[columns_index]}"):
                df[columns[columns_index]].loc[indices[indices_index]] = np.min(values)
                columns_index += 1
        if columns_index > 18:
            columns_index = 0
            indices_index += 1
        if indices_index > 5:
            indices_index = 0
        index += 1
    return df
def extract_min_column_ms(data, column_name):
    indices = ["Dot Product", "Weighted Dot Product", "Bhat 1", "Entropy Similarity", "All"]
    df = pd.DataFrame({
        "Algorithm": indices,
    })
    df.set_index("Algorithm", inplace=True)
    index = 0
    df["None"] = np.nan
    df["Exponential"] = np.nan
    df["Exponential 2"] = np.nan
    df["Logarithmic"] = np.nan
    df["sin"] = np.nan
    df["sin^2"] = np.nan
    algorithm_index = 0
    for values in data[column_name]:
        values = [value for value in values if type(value) == float]
        if data["Comparison_type"][index].endswith("_None"):
            df["None"].loc[indices[algorithm_index]] = np.min(values)
        elif data["Comparison_type"][index].endswith("_exponential"):
            df["Exponential"].loc[indices[algorithm_index]] = np.min(values)
        elif data["Comparison_type"][index].endswith("_exponential2"):
            df["Exponential 2"].loc[indices[algorithm_index]] = np.min(values)
        elif data["Comparison_type"][index].endswith("_logarithmic"):
            df["Logarithmic"].loc[indices[algorithm_index]] = np.min(values)
        elif data["Comparison_type"][index].endswith("_sin"):
            df["sin"].loc[indices[algorithm_index]] = np.min(values)
        elif data["Comparison_type"][index].endswith("_sin2"):
            df["sin^2"].loc[indices[algorithm_index]] = np.min(values)
        index += 1
        algorithm_index += 1
        if algorithm_index > 4:
            algorithm_index = 0
    return df
def extract_avg_column_ms(data, column_name):
    indices = ["Dot Product", "Weighted Dot Product", "Bhat 1", "Entropy Similarity", "All"]
    df = pd.DataFrame({
        "Algorithm": indices,
    })
    df.set_index("Algorithm", inplace=True)
    index = 0
    df["None"] = np.nan
    df["Exponential"] = np.nan
    df["Exponential 2"] = np.nan
    df["Logarithmic"] = np.nan
    df["sin"] = np.nan
    df["sin^2"] = np.nan
    algorithm_index = 0
    for values in data[column_name]:
        values = [value for value in values if type(value) == float]
        if data["Comparison_type"][index].endswith("_None"):
            df["None"].loc[indices[algorithm_index]] = np.average(values)
        elif data["Comparison_type"][index].endswith("_exponential"):
            df["Exponential"].loc[indices[algorithm_index]] = np.average(values)
        elif data["Comparison_type"][index].endswith("_exponential2"):
            df["Exponential 2"].loc[indices[algorithm_index]] = np.average(values)
        elif data["Comparison_type"][index].endswith("_logarithmic"):
            df["Logarithmic"].loc[indices[algorithm_index]] = np.average(values)
        elif data["Comparison_type"][index].endswith("_sin"):
            df["sin"].loc[indices[algorithm_index]] = np.average(values)
        elif data["Comparison_type"][index].endswith("_sin2"):
            df["sin^2"].loc[indices[algorithm_index]] = np.average(values)
        index += 1
        algorithm_index += 1
        if algorithm_index > 4:
            algorithm_index = 0
    return df

def plot_heatmap_all_alg(df, name):
    sns.heatmap(df, annot=True, cmap="coolwarm")
    plt.title(name)
    plt.show()
    return

def comparison_dtb_named_files(settings, processed = True):

    analyte_folders = ["peaks_dtb_comp"]
    comparison_types_ms = ["dot_product", "weighted_dot_product", "bhat1", "entropy_sim", "all"]
    comparison_types_dad = ["pearson", "dot_product", "ft", "derivative", "all"]
    weighting_fct_types = ["None", "exponential", "exponential2", "logarithmic", "sin", "sin2"]
    file_headings = ["Comparison_type", "Similarity_true", "Similarity_false", "Diff_distribution",
                        "lowest_diff", "avg_diff"]
    lower_limit_ft_range = range(0, 30, 5)
    upper_limit_ft_range = range(25, 500, 25)

    for folder in analyte_folders:
        csv_filename = "ms_optimization_" + folder + str(processed) + ".csv"
        csv_filepath = os.path.join(settings["directory_project"], "Data_examples", csv_filename)
        with open(csv_filepath, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=file_headings)
            writer.writeheader()
        for weighting_fct in weighting_fct_types:
            settings["ms_weighting_fct"] = weighting_fct
            for comparison_type in comparison_types_ms:
                settings["ms_algorithm"] = comparison_type
                scores = compare_analyte_dtb_ms(folder, settings, processed)
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
        csv_filename = "dad_optimization_" + folder + str(processed) + ".csv"
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
                        scores = compare_analyte_dtb_dad(folder, settings, processed)
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
                    scores = compare_analyte_dtb_dad(folder, settings, processed)
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
                            scores = compare_analyte_dtb_dad(folder, settings, processed)
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
                scores = compare_analyte_dtb_dad(folder, settings, processed)
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