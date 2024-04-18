import os
import output_files as out_files
import output_visualization as out_vis
import datapane as dp
import altair as alt


def dtb_molecule_list(settings):
    directory_project = settings["directory_project"]
    directory_database = os.path.join(directory_project, "Data_examples", "database")
    database_list = os.listdir(directory_database)

    dtb_nr = 1
    for entry in database_list:
        entry_path = os.path.join(directory_database, entry)
        if os.path.isfile(entry_path) and entry_path.lower().endswith(".cdf"):
            molecule_name = out_files.dtb_molecule_name(entry_path)
            methods_list = out_files.dtb_molecule_methods(entry_path)
            print("----------")
            print("File: " + entry)
            print("Entry Nr.: " + str(dtb_nr))
            print("Molecule name: " + molecule_name)
            print("Available methods: " + str(methods_list))
            print("----------")
            dtb_nr += 1
        else:
            print("The file " + entry + " is not a valid database entry.")
    return

def dtb_molecule_full_data(file_name, settings, processed = True):
    directory_project = settings["directory_project"]
    file_path = os.path.join(directory_project, "Data_examples", "database", file_name)

    ms_spectrum, dad_spectrum = out_files.dtb_molecule_spectra(file_path, processed)
    methods_list = out_files.dtb_molecule_methods(file_path)

    print("-----------")
    print("File Name: " + file_name)
    print("Available methods: " + str(methods_list))
    print("Processed spectra: " + str(processed))
    out_vis.plot_ms_spectrum(ms_spectrum)
    out_vis.plot_dad_spectrum(dad_spectrum)
    return

def create_analysis_report(settings, peak_folder = None, report_name = None):
    """
    Create output HTML file including visual output.
    :param settings:
    :param peak_folder: Needs to be set as a string.
    :param report_name: Needs to be set as a string.
    :return:
    """
    if peak_folder == None:
        peak_folder = settings["peak_folder_time"]
    if report_name == None:
        report_name = settings["peak_folder_time"]

    peaks_directory = os.path.join(settings["directory_project"], "Data_examples", "Peak_files", peak_folder)
    peak_list = os.listdir(peaks_directory)
    reports_directory = os.path.join(settings["directory_project"], "Data_examples", "Reports")
    report_path = os.path.join(reports_directory, report_name) + ".html"

    for peak_name in peak_list:
        peak_path = os.path.join(peaks_directory, peak_name)
        ms_spectrum, dad_spectrum = out_files.peak_spectra(peak_path, False)
        ms_spec_proc, dad_spec_proc = out_files.peak_spectra(peak_path, True)

        ms_chart = out_vis.altair_plot_ms(ms_spectrum, False)
        dad_chart = out_vis.altair_plot_dad(dad_spectrum, False)
        ms_chart_proc = out_vis.altair_plot_ms(ms_spec_proc, True)
        dad_chart_proc = out_vis.altair_plot_dad(dad_spec_proc, True)
        # spectra_unproc = alt.vconcat(ms_chart, dad_chart).resolve_scale(color="independent")
        # spectra_proc = alt.vconcat(ms_chart_proc, dad_chart_proc).resolve_scale(color="independent")

        molecule_name = out_files.peak_molecule_name(peak_path)

    select_spectra = dp.Select(blocks=[
        dp.Plot(ms_chart, label="Unprocessed MS Spectrum"),
        dp.Plot(ms_chart_proc, label="Processed MS Spectrum"),
        dp.Plot(dad_chart, label="Unprocessed DAD Spectrum"),
        dp.Plot(dad_chart_proc, label="Processed DAD Spectrum")
    ])

    select_table = dp.Select(blocks=[
        dp.DataTable(ms_spectrum.data, label="Unprocessed MS Spectrum"),
        dp.DataTable(ms_spec_proc.data, label="Processed MS Spectrum"),
        dp.DataTable(dad_spectrum.data, label="Unprocessed DAD Spectrum"),
        dp.DataTable(dad_spec_proc.data, label="Processed DAD Spectrum")
    ])


    report = dp.Report(

            dp.Text("File: " + peak_name),
            dp.Text("Molecule: " + molecule_name),
            select_spectra,
            select_table,

    )
    report.save(path=report_path, open=True)
    return