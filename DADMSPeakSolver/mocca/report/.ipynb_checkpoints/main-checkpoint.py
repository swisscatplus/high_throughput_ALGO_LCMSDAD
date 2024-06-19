# -*- coding: utf-8 -*-
"""
Created on Wed May 25 08:39:04 2022

@author: HaasCP
"""

from ...report.hplc_input import report_hplc_input
from ...report.gradient import report_gradients
from ...report.chroms import report_chroms
from ...report.bad_chroms import report_bad_chroms
from ...report.results import report_runs
from ...report.peaks import report_peaks
from ...report.quali_comps import report_quali_comps
from ...report.quant_comps import report_quant_comps
from ...report.parafac import report_parafac


def report(camp, export_path=''):
    """
    Consolidated report function.
    """
    report_hplc_input(camp.hplc_inputs, export_path)
    report_gradients(camp.hplc_inputs, export_path)
    report_peaks(camp.peak_db, export_path)
    report_chroms(camp.chroms, camp.settings, export_path)
    report_bad_chroms(camp.chroms, camp.settings, export_path)
    report_runs(camp.chroms, camp.quali_comp_db, camp.quant_comp_db, export_path)
    report_parafac(camp.chroms, export_path)
    report_quali_comps(camp.quali_comp_db, export_path)
    report_quant_comps(camp.quant_comp_db, export_path)
