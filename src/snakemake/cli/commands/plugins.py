"""Click subcommand for listing installed Snakemake plugins."""

import click

from snakemake.cli.plugin_adapter import (
    PLUGIN_REGISTRIES,
    get_install_path,
    get_registry,
)


def _print_plugin(plugin, module_prefix, verbose=False):
    """Print info for a single plugin."""
    install_path = get_install_path(plugin.name, module_prefix)
    click.echo(f"  {click.style(plugin.name, bold=True)}")

    if install_path == "built-in":
        click.echo(f"    {click.style('(built-in)', dim=True)}")
    else:
        click.echo(f"    Package:  {module_prefix}{plugin.name.replace('-', '_')}")
        click.echo(f"    Path:     {install_path}")

    if not verbose:
        return

    settings_info = plugin.get_settings_info()
    if not settings_info:
        click.echo("    Options:  (none)")
        return

    click.echo("    Options:")
    for info in settings_info:
        line = f"      {info['cliarg']}"
        parts = []
        if info.get("type"):
            type_val = info["type"]
            parts.append(
                type_val.__name__ if hasattr(type_val, "__name__") else str(type_val)
            )
        if info.get("choices"):
            parts.append(f"choices: {','.join(str(c) for c in info['choices'])}")
        if info.get("required"):
            parts.append("required")
        elif info.get("default") is not None:
            parts.append(f"default: {info['default']}")
        if info.get("env_var"):
            parts.append(f"env: {info['env_var']}")
        if parts:
            line += f"  ({'; '.join(parts)})"
        click.echo(line)
        if info.get("help"):
            click.echo(f"        {info['help']}")


@click.group()
def plugins():
    """Manage and inspect Snakemake plugins."""


@plugins.command(name="list")
@click.option(
    "-t",
    "--type",
    "plugin_types",
    multiple=True,
    type=click.Choice(list(PLUGIN_REGISTRIES.keys()), case_sensitive=False),
    help="Only show plugins of this type. Can be specified multiple times.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Show plugin CLI options and settings.",
)
def list_plugins(plugin_types, verbose):
    """List all detected Snakemake plugins."""
    types_to_show = plugin_types if plugin_types else PLUGIN_REGISTRIES.keys()
    found_any = False

    for type_name in types_to_show:
        registry = get_registry(type_name)
        if registry is None:
            continue

        plugin_names = registry.get_registered_plugins()
        if not plugin_names:
            continue

        found_any = True
        click.echo(
            f"\n{click.style(f'{type_name} plugins', fg='green', bold=True)} "
            f"({len(plugin_names)} installed)"
        )
        click.echo(f"{'─' * 40}")

        for name in sorted(plugin_names):
            plugin = registry.get_plugin(name)
            _print_plugin(plugin, registry.module_prefix, verbose=verbose)
            click.echo()

    if not found_any:
        click.echo("No plugins detected.")
