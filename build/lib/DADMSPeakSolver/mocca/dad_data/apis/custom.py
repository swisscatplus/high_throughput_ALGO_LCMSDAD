# flake8: noqa
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 15:17:53 2021

@author: haascp
"""
import numpy as np
import pandas as pd
from mocca.dad_data.utils import df_to_array



def read_custom_data(experiment):  # modified function to deal with our data
    """
    Returns the given custom data without any preprocessing
    """
    if experiment.custom_data is None:
        raise AttributeError("Custom data has to be given if data should be "
                             "processed with custom hplc_system_tag.")
    custom_data = experiment.custom_data
    return read_dad_object(custom_data)

def read_dad_object(dad_object):  # fct to read our data within mocca

    # absorbance = np.swapaxes(dad_object["DAD spectra"][:, 1], 0, 1)
    transformed_absorbance = np.array(dad_object["DAD spectra"][:,1].tolist())
    time = list(dad_object["time"])
    wvl_array = dad_object["DAD spectra"][0][0]
    wavelength = wvl_array.tolist()
    # print(transformed_absorbance.reshape(len(transformed_absorbance), len(wavelength)).shape)

    df = pd.DataFrame(transformed_absorbance.reshape(len(transformed_absorbance), len(wavelength)), columns=wavelength)
    df.insert(0, "time", time)  # ignore warning that time is a list, seems to work in mocca
    df = preprocess_df(df)

    df = pd.melt(df, id_vars='time', value_vars=df.columns[1:],
                 var_name='wavelength', value_name='absorbance')
    df['wavelength'] = df['wavelength'].astype(float)
    data, time, wavelength = df_to_array(df)
    return data, time, wavelength

def preprocess_df(df):
    """
    Preprocesses the df time column to be in line with the Chemstation API.
    """
    acq_time = df.time.max() / len(df)

    # generate new time column
    time_series = pd.Series(range(1, (len(df) + 1))).astype(float) * acq_time / 60
    df['time'] = time_series
    return df
