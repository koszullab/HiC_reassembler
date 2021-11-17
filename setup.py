import re
from setuptools import setup, find_packages

with open("hiscram/__init__.py", "r") as init:
    init_conts = init.read()
    vers_regex = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version = re.search(vers_regex, init_conts, re.MULTILINE)
    if not version:
        raise RuntimeError(f"Cannot find version string in __init__.py")

setup(
    name="HiC_reassembler",
    version=version,
    packages=find_packages(),
    entry_points={"console_scripts": ["hiscram=hiscram.cli:entry_point"]},
)
