from matplotlib import pyplot as plt

def plot_ms_spectrum(ms_spectrum):
    plt.vlines(ms_spectrum.data["m/z"], ymin=0, ymax=ms_spectrum.data["Intensity"])
    plt.title(ms_spectrum.info["Molecule"])
    plt.show()
    return

def plot_dad_spectrum(dad_spectrum):
    plt.plot(dad_spectrum.data["Wavelength"], dad_spectrum.data["Intensity"])
    plt.title(dad_spectrum.info["Molecule"])
    plt.show()
    return