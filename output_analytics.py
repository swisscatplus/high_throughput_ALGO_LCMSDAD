import os
from netCDF4 import Dataset
import pandas as pd
import data_processing as dpr
import data_processing_dad as dpr_dad
import initialize as init
import output_files as out_files
import numpy as np

def determine_best_runs(peak_path):
    all_runs = out_files.all_runs_details(peak_path)
    number_analysis_runs = 3  # Sth like 3 should be reasonable. Can be adjusted at will.
    analysis_runs_index = []
    analysis_runs_nr = []

    runs_properties = list(zip(all_runs["Pure"], all_runs["Relative Purity"], all_runs["Integral"]))

    if len(runs_properties) <= number_analysis_runs:
        for i in range(len(runs_properties)):
            analysis_runs_index.append(i)
        for i in analysis_runs_index:
            analysis_runs_nr.append(int(all_runs["Run Nr."][i]))
        return analysis_runs_nr  # Less runs than we want anyway, so just analyse all

    best_runs_index = pareto_optimality(runs_properties, 0)

    if len(best_runs_index) >= number_analysis_runs:
        for i in range(number_analysis_runs):
            analysis_runs_index.append(best_runs_index[i])  # Simply taking the first ones (as many as wanted). They should all be equal anyway
    else:  # This should be the most common scenario.
        max_superior_count = 1
        while len(best_runs_index) < number_analysis_runs:
            add_best_run_index = pareto_optimality(runs_properties, max_superior_count)
            best_runs_index.extend(add_best_run_index)
            max_superior_count += 1
        for i in range(number_analysis_runs):
            analysis_runs_index.append(best_runs_index[i])  # Just in case we now have more analysis runs than the max given

    for i in analysis_runs_index:
        analysis_runs_nr.append(int(all_runs["Run Nr."][i]))

    return analysis_runs_nr

def pareto_optimality(runs_properties, max_superior_count):
    superior_count = [0] * len(runs_properties)
    for i, run_prop in enumerate(runs_properties):  # Pareto optimality function
        for j, sec_run_prop in enumerate(runs_properties):
            if i != j and compare_runs(run_prop, sec_run_prop):
                superior_count[i] += 1
    best_run_index = []
    index = 0
    for i in superior_count:
        if i == max_superior_count:
            best_run_index.append(index)  # Adds to analysis runs only one time, index gives index in run_properties (equal to index in superior_count)
        index += 1
    return best_run_index

def compare_runs(run1_prop, run2_prop):
    """
    Returns True of the second run is valued better than the first one.
    """
    ratio_rel_purity = run1_prop[1]/run2_prop[1]
    ratio_integral = run1_prop[2]/run2_prop[2]
    if run1_prop[0] and not run2_prop[0]:
        return False
    elif run2_prop[0] and not run1_prop[0]:
        return True
    elif run1_prop[0] and run2_prop[0]:
        if run1_prop[1] >= run2_prop[1] and run1_prop[2] >= run2_prop[2]:
            return False
        elif run1_prop[1] < run2_prop[1] and run1_prop[2] < run2_prop[2]:
            return True
        elif run1_prop[1] >= run2_prop[1] and run1_prop[2] < run2_prop[2]:
            if ratio_rel_purity >= (1/ratio_integral):  # Ratio for comparison alias no bias between the two
                return False
            else:
                return True
        elif run1_prop[1] < run2_prop[1] and run1_prop[2] >= run2_prop[2]:
            if ratio_integral >= (1/ratio_rel_purity):
                return False
            else:
                return True
        else:
            return SyntaxError("Wrong syntax of run properties.")
    else:
        return SyntaxError("Wrong syntax of run properties.")

def fill_all_runs(all_runs_details, number_of_runs):
    total_runs_list = [number+1. for number in range(number_of_runs)]

    for number in total_runs_list:
        if number not in all_runs_details["Run Nr."].values:
            new_row = {"Run Nr.": number,
                       "Retention Time": np.NaN,
                       "Pure": False,
                       "Relative Purity": np.NaN,
                       "Integral": np.NaN,
                       }
            all_runs_details = pd.concat([all_runs_details, pd.DataFrame([new_row])], ignore_index=True)
    all_runs_details_sorted = all_runs_details.sort_values(by="Run Nr.")
    all_runs_details_sorted.reset_index(drop=True, inplace=True)
    return all_runs_details_sorted
