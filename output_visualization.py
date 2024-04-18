from matplotlib import pyplot as plt
import altair as alt
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io
import base64


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

def altair_plot_ms(ms_spectrum, processed):
    brush = alt.selection_interval(encodings=["x", "y"])
    if processed:
        base_chart = alt.Chart(ms_spectrum.data).mark_bar(size=2).encode(  # size to control bar line width
            x= alt.X("m/z", title="m/z"),
            y= alt.Y("Intensity", title="Intensity")
        ).properties(
            title = alt.TitleParams(
                text='Processed MS Spectrum',
                fontSize=20,
                color='black',
                align='center'
            )
        )
    else:
        base_chart = alt.Chart(ms_spectrum.data).mark_bar(size=2).encode(
            x=alt.X("m/z", title="m/z"),
            y=alt.Y("Intensity", title="Intensity")
        ).properties(
            title=alt.TitleParams(
                text='Unprocessed MS Spectrum',
                fontSize=20,
                color='black',
                align='center'
            )
        )
    left = base_chart.add_params(brush).properties(
        width=700
    )

    right = base_chart.transform_filter(
        brush
    ).properties(
        width=350,
        title=alt.TitleParams(
            text='Zoom',
            fontSize=20,
            color='black',
            align='center'
        )
    )

    return alt.hconcat(left, right)

def altair_plot_dad(dad_spectrum, processed):
    brush = alt.selection_interval(encodings=["x", "y"])

    if processed:
        base_chart = alt.Chart(dad_spectrum.data).mark_line().encode(
            x= alt.X("Wavelength", title = "Wavelength [nm]"),
            y= alt.Y("Intensity", title = "Intensity")
        ).properties(
            title=alt.TitleParams(
                text='Processed DAD Spectrum',
                fontSize=20,
                color='black',
                align='center'
            )
        )
    else:
        base_chart = alt.Chart(dad_spectrum.data).mark_line().encode(
            x=alt.X("Wavelength", title="Wavelength [nm]"),
            y=alt.Y("Intensity", title="Intensity")
        ).properties(
            title=alt.TitleParams(
                text='Unprocessed DAD Spectrum',
                fontSize=20,
                color='black',
                align='center'
            )
        )
    left = base_chart.add_params(brush).properties(
        width=700
    )

    right = base_chart.transform_filter(
        brush
    ).properties(
        width=350,
        title=alt.TitleParams(
            text='Zoom',
            fontSize=20,
            color='black',
            align='center'
        )
    )
    return alt.hconcat(left, right)

def draw_molecule(inchi):
    molecule = Chem.MolFromInchi(inchi)
    image = Draw.MolToImage(molecule)

    img_byte_arr = io.BytesIO()
    image.save(img_byte_arr, format='PNG')
    img_byte_arr = img_byte_arr.getvalue()
    b64_string = base64.b64encode(img_byte_arr).decode()
    html_image = f'<img src="data:image/png;base64,{b64_string}" alt="Molecule Image" style="width: 30%; height: auto;">'

    return html_image

def altair_plot_chrom_ms(chrom_ms):
    brush = alt.selection_interval(encodings=["x", "y"])

    base_chart = alt.Chart(chrom_ms).mark_line().encode(
        x=alt.X("time", title="Time [s]"),
        y=alt.Y("total intensity", title="Intensity [Count]")
    ).properties(
        title=alt.TitleParams(
            text='MS Chromatogram',
            fontSize=20,
            color='black',
            align='center'
        )
    )

    left = base_chart.add_params(brush).properties(
        width=700
    )

    right = base_chart.transform_filter(
        brush
    ).properties(
        width=350,
        title=alt.TitleParams(
            text='Zoom',
            fontSize=20,
            color='black',
            align='center'
        )
    )
    return alt.hconcat(left, right)

def altair_plot_chrom_dad(chrom_dad):
    brush = alt.selection_interval(encodings=["x", "y"])

    base_chart = alt.Chart(chrom_dad).mark_line().encode(
        x=alt.X("time", title="Time [s]"),
        y=alt.Y("total intensity", title="Intensity [mAU]")
    ).properties(
        title=alt.TitleParams(
            text='DAD Chromatogram',
            fontSize=20,
            color='black',
            align='center'
        )
    )

    left = base_chart.add_params(brush).properties(
        width=700
    )

    right = base_chart.transform_filter(
        brush
    ).properties(
        width=350,
        title=alt.TitleParams(
            text='Zoom',
            fontSize=20,
            color='black',
            align='center'
        )
    )
    return alt.hconcat(left, right)
