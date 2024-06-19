from setuptools import setup, find_packages

setup(
    name='DADMSPeakSolver',
    version='0.1.3',
    author='Leander Choudhury',
    author_email='leander.choudhury@epfl.ch',
    description='Program for the automated analysis and sample selection of high-throughout LC-DAD-MS data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/swisscatplus/high_throughput_ALGO_LCMSDAD',
    packages=find_packages(),  # Automatically find packages in the directory
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',  # Specify the required Python version
    install_requires=[
        'numpy==1.26.4',
        'pandas==1.5.3',
        'scipy==1.13.1',
        'matplotlib==3.8.4',
        'scikit_learn==1.4.2',
        'requests==2.31.0',
        'netCDF4==1.6.5',
        'rdkit==2023.9.5',
        'PyWavelets==1.6.0',
        'tensorly==0.8.1',
        'seaborn==0.13.2',
        'altair==5.3.0',
        'datapane==0.17.0',
        'dill==0.3.8',
        'h5py==3.11.0',
        'importlib_metadata==7.1.0',
        'rdflib==7.0.0',
        'Pillow==10.3.0',
    ],
    include_package_data=False,
)
