#!/usr/bin/env python3
"""Structural variants simulation, detection and ressembly using Hi-C contacts.
"""

import re
import codecs
from setuptools import setup, find_packages

with open("hiscram/__init__.py", "r") as init:
    init_conts = init.read()
    vers_regex = r"^__version__ = ['\"]([^'\"]*)['\"]"
    VERSION = re.search(vers_regex, init_conts, re.MULTILINE)
    if not VERSION:
        raise RuntimeError(f"Cannot find version string in __init__.py")

DESCRIPTION = __doc__.strip("\n")

with codecs.open("README.md", encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

with open("requirements.txt", "r") as f:
    REQUIREMENTS = f.read().splitlines()


CLASSIFIERS = [
    "Development Status :: 2 - Pre-Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Operating System :: OS Independent",
]


setup(
    name="HiC_reassembler",
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS,
    packages=find_packages(),
    entry_points={"console_scripts": ["hiscram=hiscram.cli:entry_point"]},
    install_requires=REQUIREMENTS,
)
