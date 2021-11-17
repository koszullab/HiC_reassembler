import os
import random
import tempfile
import pytest
import numpy as np
from hiscram.scramble import sv


@pytest.fixture
def tmp_file_rw():
    tmp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    yield tmp_file
    os.unlink(tmp_file.name)


@pytest.fixture
def mixer():
    np.random.seed(10)
    random.seed(10)
    fname = tempfile.NamedTemporaryFile(delete=False).name
    with open(fname, "w") as fh:
        fh.write(">chr1\nAACCCAAACC\n")
        fh.write(">chr2\nGTGTGT\n")
    yield sv.Mixer(fname, config_profile="Debug")


def test_get_random_region(mixer):
    frag = sv.get_random_region(mixer.genome, region_size=3)
    assert len(frag) == 3
    with pytest.raises(ValueError):
        sv.get_random_region(mixer.genome, region_size=300)


def test_mixer(mixer):
    assert len(mixer.genome.chroms) == 2


def test_mixer_load_profile(mixer):
    """validate the debug profile (correct entries, probabilities sum to 1 and
    expected structure)."""
    assert list(mixer.config.keys()) == ["SV_freq", "SV_types"]
    assert isinstance(mixer.config["SV_freq"], float)
    tot_probs = 0.0
    for sv in mixer.config["SV_types"].keys():
        sv_cfg = mixer.config["SV_types"][sv]
        assert all(
            [
                val in sv_cfg.keys()
                for val in ["prop", "mean_size", "sd_size", "code"]
            ]
        )
        tot_probs += sv_cfg["prop"]
    assert tot_probs == 1.0


def test_mixer_generate_sv(mixer):
    n_sv = mixer.generate_sv()
    assert sum(n_sv.values()) == 2
    assert len(mixer.genome.get_breakpoints()) == 2


def test_mixer_write_genome(mixer, tmp_file_rw):
    """Check if the genome produced by mixer matches the SV introduced."""
    mixer.generate_sv()
    mixer.write_genome(tmp_file_rw)
    tmp_file_rw.close()
    expected_genome = ">chr1\nATTTGCC\n>chr2\nGTGTGT\n"
    with open(tmp_file_rw.name, "r") as fh:
        assert fh.read() == expected_genome


def test_mixer_write_breakpoints(mixer, tmp_file_rw):
    """Validate bedpe files produced by mixer.write_breakpoints."""
    mixer.generate_sv()
    mixer.write_breakpoints(tmp_file_rw)
    tmp_file_rw.close()
    expected_bedpe = "chr1\t1\t1\tchr1\t4\t4\t.\t.\t+\t+\nchr1\t8\t8\tchr1\t8\t8\t.\t.\t-\t-\n"
    with open(tmp_file_rw.name, "r") as fh:
        assert fh.read() == expected_bedpe