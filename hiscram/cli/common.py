import json
import click
import hiscram


@click.command()
def list_profiles():
    """Display the names of available SV profiles for scrambling."""
    format_settings = lambda stg: (
        f'SV freq={stg["SV_freq"]},\tSV types=['
        f'{", ".join(list(stg["SV_types"].keys()))}]'
    )
    with open(hiscram.SCRAMBLE_CONFIG_PATH, "r") as cfg:
        profiles = [
            f"{name}:\t{format_settings(stg)}"
            for name, stg in json.load(cfg).items()
        ]

    profiles = "\n".join(profiles)
    click.echo(profiles)
