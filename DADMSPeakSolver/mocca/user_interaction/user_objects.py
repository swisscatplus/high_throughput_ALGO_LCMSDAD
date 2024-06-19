#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 10:46:55 2022

@author: haascp
"""

import os
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List

from mocca.dad_data.models import GradientData



@dataclass()
class Gradient():
    """
    Data container to store user input regarding gradients.
    """
    path : str
    dataset : GradientData = field(init=False)

    def __post_init__(self):
        if not os.path.exists(self.path):
            raise ValueError(f"Given gradient path {self.path} does not exist.")

    def __repr__(self):
        return f"Gradient({self.path})"


@dataclass()
class Compound():
    """
    Data container to store user input regarding added compounds.
    """
    key: str
    conc: Optional[float] = None
    # following are only for ordering batch of runs
    is_solvent: bool = False
    is_istd: bool = False


@dataclass()
class InternalStandard():
    """
    Data container to store user input regarding added internal standards.
    """
    key: str
    conc: Optional[float] = None


@dataclass()
class CustomData():
    """
    Data container to store custom data like, e.g., from HPLC chromatogram
    simulations.
    """
    data: np.ndarray
    time: list
    wavelength: list

    def __post_init__(self):
        self._check_custom_data()

    def _check_custom_data(self):
        if (self.data.shape[0] != len(self.wavelength) or
                self.data.shape[1] != len(self.time)):
            raise ValueError("Data must be given as a two-dimensional numpy "
                             "ndarray with the shape (len(wavelenght), "
                             "len(time))")


@dataclass()
class HplcInput():
    """
    Data container to store user input.
    """
    path : str
    gradient : Optional[Gradient]
    compound: Optional[Compound] = None
    istd: Optional[List[InternalStandard]] = None
    processed: bool = False
    custom_data: CustomData = None

    def __post_init__(self):
        """if self.custom_data is None and not os.path.exists(self.path):
            raise ValueError(f"Given path {self.path} does not exist.")"""  # This is not necessary, we work with custom_data
        if self.istd is not None and type(self.istd) != list:
            self.istd = [self.istd]
        if self.compound and self.istd:
            if self.compound.is_solvent:
                raise ValueError("Solvent run has an internal standard added. Use "
                                 "solvent == True only for pure solvent runs. These "
                                 "runs will be analyzed first and should cover the "
                                 "case that all samples are recorded in a UV-Vis "
                                 "active solvent. Solvents can also be added as "
                                 "compounds later on with solvent == False.")

            if self.compound.is_istd and self.compound.key == self.istd.key:
                raise ValueError("Internal standard cannot be analyzed relative "
                                 "to itself. If the internal standard should be "
                                 "added as compound, do not give internal "
                                 "standard parameter. If a run containing "
                                 "internal standard should be analyzed do not "
                                 "give internal standard as a compound.")
