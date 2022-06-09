#!/usr/bin/env python3
from setuptools import setup, find_packages
import pathlib

setup(
    name="beta_mining",
    version="1.0.0",
    author="Paula Berry",
    author_email="paula.l.berry@gmail.com, pberry@stowers.org",
    scripts=["bin/beta_mining"],
    url="https://github.com/paulaberry/beta_mining",
    license="GPL v3",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Academic Researchers",
        "License :: GPL v3",
        "Programming Language :: Python :: 3",

    ],
    packages=find_packages(include=["bin/beta_mining", "beta_mining/beta_mining_algorithm.py", "beta_mining/beta_mining_functions.py"]),
    python_requires=">=3",
    install_requires=[
        "biopython>=1.79",
        "biopandas>=0.4.1",
        "json5>=0.9.6",
        "numpy>=1.21",
        "pandas>=1.4.1",
        "ProDy>=2.0.2",
        "python-dateutil>=2.8.2",
        "PyYAML>=5.3.1",
        "scikit-image>=0.19.2",
        "scikit-learn>=1.0.2",
        "scipy>=1.8.0",
        "simplejson>=3.16.0",
    ],
    data_files=[("", ["beta_mining/default_config.yaml", "beta_mining/structure_dictionary.json"])],

    entry_points={
    "console_scripts": [
        "beta_mining=beta_mining:main",
    ],
    },
    )
