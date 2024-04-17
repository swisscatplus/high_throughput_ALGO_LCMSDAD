# flake8: noqa
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 13:22:39 2021

@author: haascp
"""

from mocca.chromatogram.utils import check_same_dataset

class Chromatogram():
    """
    Class containing data regarding chromatograms.
    """
    def __init__(self, experiment, dataset):
        self.experiment = experiment
        self.dataset = dataset
        self.peaks = []
        self.warnings = []
        self.bad_data = False
        self.parafac_models = []

    def __iter__(self):
        """
        Yields all items inside the database.
        Allows statements such as "for component in component_database:"
        """
        for peak in self.peaks:
            yield peak

    def __contains__(self, compound_id):
        """
        Allows statements such as "'Component 1' in component_database"
        to see if that Component is inside the database.
        """
        return compound_id in [peak.compound_id for peak in self.peaks]

    def __getitem__(self, compound_id):
        """
        Allows statements such as "component_database['Component 1']" to access
        that Component.
        """
        if compound_id in self:
            return next(peak for peak in self.peaks if
                        peak.compound_id == compound_id)
        else:
            raise AttributeError("{} not found in Database!"
                                 "".format(compound_id))
    
    def __eq__(self, other):
        """
        Checks if two data sets are the same based on its underlying dataset.
        """
        if not isinstance(other, Chromatogram):
            # don't attempt to compare against unrelated types
            raise ValueError("Both chromatograms must be of the mocca "
                             "Chromatogram type!")
        return self.dataset == other.dataset

    def insert_peak(self, peak):
        """
        Adds a peak to the chromatogram.
        """
        if self.dataset:
            check_same_dataset(peak, self)
        self.peaks.append(peak)
