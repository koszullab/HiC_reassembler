import pathlib

__version__ = "0.0.1"
__author__ = "Cyril Matthey-Doret, Jules Olayé, Amaury Bignaud"
__copyright__ = "Copyright © 2017-2018, Institut Pasteur, Paris, France"
__LICENSE__ = "MIT"
SCRAMBLE_CONFIG_PATH = str(
    pathlib.Path(__file__).parent / "config" / "sv_profiles.json"
)
