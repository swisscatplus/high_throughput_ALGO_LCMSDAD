import os
import requests
from netCDF4 import Dataset
from rdkit import Chem
from rdkit.Chem import AllChem

"""
Creates 4 different types of files:
Database file - All different analytics to describe a registered molecule. Created from a "Peak file".
Method file - Contains the background spectra for a given method.
Peak file - A peak, found in one or multiple analysis runs, contains data and if applicable information on 
associated molecules.
Signal file - Any peak found in a analysis run, used for comparison with dtb and other runs.
Run file - A summary of the peaks of a given run and their association.
"""

def create_dtb_entry(directory_project, peak_file_path, molecule_name, molecule_inchi = None, molecule_smiles = None):
    """
    Fct to create a background file including the DAD and MS spectrum for a molecule, as well as the information
    on the method used and corresponding Rf values.
    :param peak_file_path:
    :param dtb_folder_path:
    :param name:
    :return:
    """
    directory_database = os.path.join(directory_project, "Data_examples", "database")  # Database path
    peak_file = Dataset(peak_file_path, "r")

    if molecule_inchi is not None:  # If both, inchi and smiles are entered, priority lies at InChI.
        try:
            inchi_key = AllChem.InchiToInchiKey(molecule_inchi) #  Problems before, need to check if this really works (name and smiles work!)
        except:
            print("ERROR: An invalid InChI was entered.")
            peak_file.close()
            return
        molecule = Chem.MolFromInchi(molecule_inchi)
        molecule_smiles = Chem.MolToSmiles(molecule)
    elif molecule_smiles is not None:
        try:
            molecule = Chem.MolFromSmiles(molecule_smiles)
            molecule_inchi = Chem.MolToInchi(molecule)
            inchi_key = AllChem.InchiToInchiKey(molecule_inchi)
        except:
            print("ERROR: An invalid SMILES was entered.")
            peak_file.close()
            return
    else:
        molecule_smiles = get_smiles_from_pubchem(molecule_name)
        try:
            molecule = Chem.MolFromSmiles(molecule_smiles)
            molecule_inchi = Chem.MolToInchi(molecule)
            inchi_key = AllChem.InchiToInchiKey(molecule_inchi)
        except:
            print("ERROR: No valid SMILES was obtained from PubChem. Try entering a SMILES or InChI directly.")
            peak_file.close()
            return

    dtb_file_name = inchi_key + ".cdf" # give inch key name either from name directly or from inchi
    dtb_file_path = os.path.join(directory_database, dtb_file_name)

    if dtb_file_name in os.listdir(directory_database): # List of all files in the database
        print("The molecule " + molecule_name + " already has an entry in the database: " + dtb_file_name)
        return
    else:
        dtb_file = Dataset(dtb_file_path, "w", format="NETCDF4")

    copy_file(peak_file, dtb_file)

    retention_time_dim = dtb_file.createDimension("Retention Time", 1)

    grp_methods = dtb_file.createGroup("Methods")

    first_method = grp_methods.createGroup(peak_file.Method)
    retention_time = first_method.createVariable("Retention Time", "f8", retention_time_dim)
    retention_time[:] = float(peak_file.time)

    dtb_file.delncattr("run_nr")
    dtb_file.delncattr("Method")
    dtb_file.delncattr("time")
    dtb_file.delncattr("peak_nr")

    dtb_file.molecule_name = molecule_name
    dtb_file.molecule_inchi = molecule_inchi
    dtb_file.molecule_inchi_key = inchi_key
    dtb_file.molecule_smiles = molecule_smiles

    peak_file.close()
    dtb_file.close()
    return

def add_method_to_dtb_entry(dtb_file_path, peak_file_path):
    """
    Fct which adds a new group to a dtb entry, containing the name of a new method and the corresponding
    retention time.
    :param dtb_file_path:
    :param peak_file_path:
    :return:
    """

    dtb_file = Dataset(dtb_file_path, "r+")
    peak_file = Dataset(peak_file_path, "r")

    grp_methods = dtb_file.groups["Methods"]

    if peak_file.Method in grp_methods.groups:
        print("The Method " + peak_file.Method + " already exists in the database file.")
    else:
        new_method = grp_methods.createGroup(peak_file.Method)
        retention_time = new_method.createVariable("Retention Time", "f8", dtb_file.dimensions["Retention Time"])
        retention_time[:] = float(peak_file.time)

    dtb_file.close()
    peak_file.close()
    return

def create_peak_file(signal_file_path, peak_folder_path, peak_nr):
    """
    Fct. to create a peak file from a signal file
    :param signal_file_path:
    :param peak_nr:
    :return:
    """
    signal_file = Dataset(signal_file_path, "r")

    # Create new file
    peak_file_path = os.path.join(peak_folder_path, ("peak_" + str(peak_nr) + ".cdf"))
    peak_file = Dataset(peak_file_path, "w", format="NETCDF4")

    copy_file(signal_file, peak_file)
    peak_file.peak_nr = str(peak_nr)  # Add Peak Nr as new attribute
    peak_file.delncattr("sgl_nr")  # Delete Signal Nr attribute
    peak_file.dtb_hits = ""
    peak_file.dtb_hits_inchi = ""
    peak_file.several_hits = "No"
    hits_different_method = peak_file.createGroup("dtb hits different method")
    hits_different_method.dtb_hits = ""
    hits_only_ms = peak_file.createGroup("dtb hits only MS")
    hits_only_ms.dtb_hits = ""
    hits_only_dad = peak_file.createGroup("dtb hits only DAD")
    hits_only_dad.dtb_hits = ""

    signal_file.close()
    peak_file.close()
    return

def create_signal_file(ms_spec, dad_spec, directory_project, plot = True):
    directory_sglfiles = os.path.join(directory_project, "Data_examples", "Signal_files", "Sgl_file")  # Database path
    file_name = directory_sglfiles + "_" + ms_spec.info["Run Nr"] + "_" + ms_spec.info["Signal Nr"] + ".cdf"

    if plot:
        print("----------")
        ms_spec.plot_spectra()
        dad_spec.plot_dad()
        print("Retention time: " + str(dad_spec.info["Time"]))

    sglfile = Dataset(file_name, 'w', format='NETCDF4')

    sglfile.Method = ms_spec.info["LC Method"]
    sglfile.run_nr = ms_spec.info["Run Nr"]
    sglfile.run_name = dad_spec.info["Run Name"]
    sglfile.sgl_nr = ms_spec.info["Signal Nr"]
    sglfile.time = dad_spec.info["Time"]
    sglfile.pure = str(dad_spec.info["Pure"])
    sglfile.relative_purity = dad_spec.info["Relative Purity"]
    sglfile.integral = dad_spec.info["Integral"]

    sglfile.corresponding_peakfile = ""

    intensity = sglfile.createDimension("intensity", None)

    grp_ms = sglfile.createGroup("MS Data")
    grp_dad = sglfile.createGroup("DAD Data")

    mass = grp_ms.createDimension("mass", None)
    grp_raw_ms = grp_ms.createGroup("Raw MS Data")
    raw_mass = grp_raw_ms.createVariable("mass", "f8", mass)
    raw_ms_intensity = grp_raw_ms.createVariable("Intensity", "f8", intensity)
    raw_mass[:] = ms_spec.raw_data.loc[:, "m/z"]
    raw_ms_intensity[:] = ms_spec.raw_data.loc[:, "Intensity"]

    grp_proc_ms = grp_ms.createGroup("Processed MS Data")
    proc_mass = grp_proc_ms.createVariable("mass", "f8", mass)
    proc_ms_intensity = grp_proc_ms.createVariable("Intensity", "f8", intensity)
    proc_mass[:] = ms_spec.data.loc[:, "m/z"]
    proc_ms_intensity[:] = ms_spec.data.loc[:, "Intensity"]

    wavelength = grp_dad.createDimension("Wavelength", None)
    ft_wvl = grp_dad.createDimension("ft_wvl", None)

    grp_raw_dad = grp_dad.createGroup("Raw DAD Data")
    raw_wvl = grp_raw_dad.createVariable("Wavelength", "f8", wavelength)
    raw_dad_intensity = grp_raw_dad.createVariable("Intensity", "f8", intensity)
    raw_wvl[:] = dad_spec.raw_data.loc[:, "Wavelength"]
    raw_dad_intensity[:] = dad_spec.raw_data.loc[:, "Intensity"]

    grp_proc_dad = grp_dad.createGroup("Processed DAD Data")
    proc_wvl = grp_proc_dad.createVariable("Wavelength", "f8", wavelength)
    proc_dad_intensity = grp_proc_dad.createVariable("Intensity", "f8", intensity)
    proc_ft_wvl = grp_proc_dad.createVariable("FT Wavelength", "f8", ft_wvl)  # delete if ft comparison not working well
    proc_wvl[:] = dad_spec.data.loc[:, "Wavelength"]
    proc_dad_intensity[:] = dad_spec.data.loc[:, "Intensity"]
    # proc_ft_wvl[:] = dad_spec.data.loc[:, "FT Wavelength"]  # Uncomment once FT has been constructed.

    sglfile.close()
    return

def copy_file(file_old, file_new):
    """
    Function to copy one netCDF file completely to another netCDF file without changing the initial file.
    :param file_old: Old file (open netCDF file)
    :param file_new: New file (open netCDF fle)
    :return:
    """

    for name, dimension in file_old.dimensions.items():
        file_new.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in file_old.variables.items():
        new_variable = file_new.createVariable(name, variable.datatype, variable.dimensions)
        new_variable.setncatts({attr: variable.getncattr(attr) for attr in variable.ncattrs()})
        new_variable[:] = variable[:]

    # Copy all attributes
    file_new.setncatts({attr: file_old.getncattr(attr) for attr in file_old.ncattrs()})

    for name, group in file_old.groups.items():
        new_group = file_new.createGroup(name)
        copy_file(group, new_group)  # Recursion step to copy all groups including their content

    return


def get_smiles_from_pubchem(compound_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200: # 200 status means request was successful.
        data = response.json() # parse response to python dictionary.
        smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES'] # extracts smiles from data
        return smiles
    else:
        return None
