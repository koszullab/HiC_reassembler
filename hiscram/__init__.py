import pathlib

__version__ = "0.0.1"
SCRAMBLE_CONFIG_PATH = str(
    pathlib.Path(__file__).parent / "scramble" / "config" / "template.json"
)
