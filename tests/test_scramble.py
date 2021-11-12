import tempfile
import pytest
from hiscram.scramble import sv


@pytest.fixture
def mixer():
    fname = tempfile.NamedTemporaryFile(delete=False).name
    with open(fname, "w") as fh:
        fh.write(">chr1\nAACCCAAACC\n")
        fh.write(">chr2\nGTGTGT\n")
    yield sv.Mixer(fname, config_profile="Debug")


def test_mixer(mixer):
    assert len(mixer.genome.chroms) == 2
