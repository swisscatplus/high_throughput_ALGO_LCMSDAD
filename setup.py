from setuptools import setup, find_packages

setup(
    name='mypackage',
    version='0.1.0',
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
        'required_package1',
        'required_package2',
    ],
    include_package_data=False,  # Include files from MANIFEST.in
)
