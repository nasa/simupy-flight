from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# get the version
exec(open(path.join(here, "simupy_flight", "version.py")).read())

# Get the long description from the README file
with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="simupy_flight",
    version=__version__,
    description="A framework for modeling and simulating dynamical systems.",
    long_description=long_description,
    packages=find_packages(),
    author="Benjamin W. L. Margolis",
    author_email="benjamin.margolis@nasa.gov",
    license="NASA Open Source Agreement Version 1.3",
    python_requires=">=3",
    install_requires=[
        "numpy>=1.11.3",
        "scipy>=0.18.1",
        "simupy>=1.1.2",
        "pyerfa",
        "fluids",
    ],
    extras_require={
        "test": ["pandas", "ndsplines", "matplotlib"],
        "derivation": ["sympy==1.4", "jupyterlab"],
        "docs": ["sphinx", "sphinx-gallery", "sphinx-rtd-theme"],
        "daveml": ["lxml", "ndsplines", "sympy==1.4"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
)
