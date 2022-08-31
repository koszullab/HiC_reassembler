# End-to-end tests of scrambling, detection and reassembly
import tempfile
import pytest
from click.testing import CliRunner
from hiscram.cli import entry_point

GENOME = "test_data/seq.fa"
R1 = "test_data/sample.reads_for.fastq.gz"
R2 = "test_data/sample.reads_rev.fastq.gz"


@pytest.fixture
def tmp_dir():
    tmp_dir = tempfile.TemporaryDirectory()
    yield tmp_dir
    # tmp_dir.cleanup()


def test_scramble(tmp_dir):
    """Test the end-to-end scrambling process."""
    runner = CliRunner()
    result = runner.invoke(
        entry_point, ["scramble", "-1", R1, "-2", R2, GENOME, tmp_dir.name]
    )
    assert result.exit_code == 0
