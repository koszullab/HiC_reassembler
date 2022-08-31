import click

from hiscram import __version__
import hiscram.cli as scramble_cli
import hiscram.cli.common as common_cli


@click.group()
@click.version_option(version=__version__)
def entry_point():
    pass


entry_point.add_command(scramble_cli.run_scrambles, name="scramble")
entry_point.add_command(common_cli.list_profiles, name="profiles")
