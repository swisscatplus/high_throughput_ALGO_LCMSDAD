import os
import output_files as out_files
import output_visualization as out_vis
import output_analytics as out_an
import datapane as dp
import altair as alt
alt.data_transformers.enable('vegafusion')

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

def create_analysis_report(settings, run_folder, peak_folder = None, report_name = None):
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

    if len(peak_list) == 0:
        print("No peaks to analyse.")
        return  # Later change so we go to run analysis directly

    reports_directory = os.path.join(settings["directory_project"], "Data_examples", "Reports")
    report_path = os.path.join(reports_directory, report_name) + ".html"

    runs_directory = os.path.join(settings["directory_project"], "Data_examples", "testfiles", run_folder)
    run_list = os.listdir(runs_directory)

    if len(run_list) == 0:
        print("No runs to analyse.")
        return

    """
    Handling of peaks
    """

    peak_report_list = []
    peak_analytics_list = []
    for peak_name in peak_list:
        peak = create_peak_report(peak_name, peaks_directory)
        peak_report_list.append(peak)
        peak_analytics = create_peak_heatmap(peak_name, peaks_directory, len(run_list), settings)
        peak_analytics_list.append(peak_analytics)

    peaks_groups = [dp.Group(*peak[0], label=peak[1]) for peak in peak_report_list]
    while len(peaks_groups) <= 1:
        peaks_groups.append(dp.Group(
            dp.Text("#")
        ))

    peaks_analytic_groups = [dp.Group(*peak[0], label=peak[1]) for peak in peak_analytics_list]
    while len(peaks_analytic_groups) <= 1:
        peaks_analytic_groups.append(dp.Group(
            dp.Text("#")
        ))

    peak_rep_select = dp.Select(
        blocks=peaks_groups,
        type=dp.SelectType.DROPDOWN
    )

    analytics_select = dp.Select(
        blocks=peaks_analytic_groups,
        type=dp.SelectType.DROPDOWN
    )

    """
    Handling of runs
    """

    run_report_list = []
    run_nr = 1
    for run_name in run_list:
        run = create_run_report(run_name, runs_directory, run_nr, peak_folder, settings)
        run_report_list.append(run)
        run_nr += 1

    runs_groups = [dp.Group(*run[0], label=("Run Nr.: " + str(run[1]))) for run in run_report_list]
    while len(runs_groups) <= 1:
        runs_groups.append(dp.Group(
            dp.Text("#")
        ))

    run_rep_select = dp.Select(
        blocks=runs_groups,
        type=dp.SelectType.DROPDOWN
    )
    """
    Handling groups to build in.
    """

    all_peaks_group = dp.Group(
        peak_rep_select,
        dp.Text("#"),
        label="Peaks"
    )
    all_runs_group = dp.Group(
        run_rep_select,
        dp.Text("#"),
        label="Runs"
    )
    all_analysis_group = dp.Group(
        analytics_select,
        dp.Text("#"),
        label="Analysis"
    )

    select_report = dp.Select(
        blocks=[
        all_peaks_group,
        all_runs_group,
        all_analysis_group
        ]
    )

    master_report = dp.Report(
        select_report
    )

    master_report.save(path=report_path, open=True)
    return

def create_peak_report(peak_name, peaks_directory):
    peak_path = os.path.join(peaks_directory, peak_name)
    ms_spectrum, dad_spectrum = out_files.peak_spectra(peak_path, False)
    ms_spec_proc, dad_spec_proc = out_files.peak_spectra(peak_path, True)

    ms_chart = out_vis.altair_plot_ms(ms_spectrum, False)
    dad_chart = out_vis.altair_plot_dad(dad_spectrum, False)
    ms_chart_proc = out_vis.altair_plot_ms(ms_spec_proc, True)
    dad_chart_proc = out_vis.altair_plot_dad(dad_spec_proc, True)
    # spectra_unproc = alt.vconcat(ms_chart, dad_chart).resolve_scale(color="independent")
    # spectra_proc = alt.vconcat(ms_chart_proc, dad_chart_proc).resolve_scale(color="independent")
    retention_time = round(out_files.peak_retention_time(peak_path), 3)

    all_runs = out_files.all_runs_details(peak_path)

    analysis_runs = out_an.determine_best_runs(peak_path)

    run_table = dp.Group(
        dp.Text("Best runs for analysis (Run Nr.): " + str(analysis_runs)),
        dp.Text("**All runs with this peak:**"),
        dp.DataTable(all_runs, label="All runs with this peak")
    )

    molecule_name = out_files.peak_molecule_name(peak_path)
    inchi = out_files.peak_inchi(peak_path)
    if not inchi == None:
        molecule_image_html = out_vis.draw_molecule(inchi)
    else:
        molecule_image_html = None

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

    if not molecule_image_html == None:
        molecule_display = dp.Group(
            dp.Group(
                dp.Text("File: " + peak_name),
                dp.Text("Molecule: " + molecule_name),
                dp.Text("Retention time [s]: " + str(retention_time)),
            ),
            dp.HTML(molecule_image_html),
            columns=2
        )
    else:
        molecule_display = dp.Group(
            dp.Text("File: " + peak_name),
            dp.Text("Molecule: " + molecule_name),
            dp.Text("Retention time [s]: " + str(retention_time)),
        )

    report = [
        molecule_display,
        select_spectra,
        run_table,
        select_table,
    ]

    return report, peak_name

def create_run_report(run_name, run_directory, run_nr, peak_folder, settings):
    run_path = os.path.join(run_directory, run_name)
    chrom_ms, chrom_dad = out_files.run_chromatograms(run_path, settings)
    ms_chart = out_vis.altair_plot_chrom_ms(chrom_ms)
    dad_chart = out_vis.altair_plot_chrom_dad(chrom_dad)
    chrom_combined = alt.vconcat(ms_chart, dad_chart).resolve_scale(color="independent")

    all_peaks = out_files.all_peaks_details(run_nr, peak_folder, settings)

    peaks_table = dp.Group(
        dp.Text("**All peaks within this run:**"),
        dp.DataTable(all_peaks, label="All peaks within this run")
    )

    chrom_group = dp.Group(
        dp.Plot(chrom_combined, label="Chromatograms")
    )
    run_display = dp.Group(
        dp.Text("File: " + run_name),
        dp.Text("Method: " + settings["method_name"]),
        dp.Text("Ion detection mode: " + settings["ion_detection_mode"])
    )

    report = [
        run_display,
        chrom_group,
        peaks_table,
    ]

    return report, run_nr

def create_peak_heatmap(peak_name, peaks_directory, number_of_runs, settings):
    peak_path = os.path.join(peaks_directory, peak_name)

    all_runs_details = out_files.all_runs_details(peak_path)
    all_runs_details = out_an.fill_all_runs(all_runs_details, number_of_runs)
    heatmap_integral = out_vis.altair_heatmap_integral(all_runs_details, settings["Number Columns"])
    heatmap__relative_purity = out_vis.altair_heatmap_relative_purity(all_runs_details, settings["Number Columns"])

    molecule_name = out_files.peak_molecule_name(peak_path)
    inchi = out_files.peak_inchi(peak_path)
    if not inchi == None:
        molecule_image_html = out_vis.draw_molecule(inchi)
    else:
        molecule_image_html = None

    select_heatmap = dp.Select(blocks=[
        dp.Plot(heatmap_integral, label="Integral Values"),
        dp.Plot(heatmap__relative_purity, label="Relative Purity Values"),
    ])

    if not molecule_image_html == None:
        molecule_display = dp.Group(
            dp.Group(
                dp.Text("File: " + peak_name),
                dp.Text("Molecule: " + molecule_name),
            ),
            dp.HTML(molecule_image_html),
            columns=2
        )
    else:
        molecule_display = dp.Group(
            dp.Text("File: " + peak_name),
            dp.Text("Molecule: " + molecule_name),
        )

    report = [
        molecule_display,
        select_heatmap,
    ]

    return report, peak_name