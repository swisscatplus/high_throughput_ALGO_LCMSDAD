import os
import output_files as out_files

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

def plot_dtb_molecule_spectra(file_name, settings, processed = False):
    directory_project = settings["directory_project"]
    file_path = os.path.join(directory_project, "Data_examples", "database", file_name)

    ms_spectrum, dad_spectrum = out_files.dtb_molecule_spectra(file_path)

    return