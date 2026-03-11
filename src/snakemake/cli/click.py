"""Click-based CLI entry point for Snakemake.

This module provides a Click group that serves as the main entry point.
When no subcommand is given, it falls through to the existing argparse-based
main() for full backward compatibility.
"""

import sys

import click
from snakemake.cli.legacy import main as legacy_main
from snakemake.cli.commands import info, clean, dagviz, unlock, lint, utils, software, new


def _find_workflow_profile(snakefile):
    """Discover workflow-specific profile relative to snakefile."""
    from pathlib import Path

    from snakemake.api import resolve_snakefile

    resolved = resolve_snakefile(snakefile, allow_missing=True)
    if resolved is None:
        return None

    default_path = Path("profiles/default")
    candidates = [default_path, Path(resolved).parent / default_path]
    for profile in candidates:
        if profile.exists():
            return str(profile)
    return None


def _load_profile_defaults(profile, workflow_profile):
    """Load profile YAML files and return merged defaults dict."""
    from snakemake.cli.legacy import get_profile_dir
    from snakemake.profiles import ProfileConfigFileParser

    defaults = {}
    for p in [profile, workflow_profile]:
        if p is None:
            continue
        entry = get_profile_dir(p)
        if entry is not None:
            _, config_file = entry
            with open(config_file) as f:
                parsed = ProfileConfigFileParser().parse(f)
                defaults.update(parsed)
    return defaults


def _extract_flag(args, *flags):
    """Pre-scan raw argv for a flag value without modifying args."""
    for flag in flags:
        try:
            idx = args.index(flag)
            if idx + 1 < len(args):
                return args[idx + 1]
        except ValueError:
            continue
    return None


class ProfileAwareGroup(click.Group):
    """Click group that loads snakemake profiles as default_map.

    When no registered subcommand is found, all args are passed through
    to the legacy argparse-based main() for full backward compatibility.
    """

    def make_context(self, info_name, args, parent=None, **extra):
        profile = _extract_flag(args, "--profile")
        snakefile = _extract_flag(args, "--snakefile", "-s")
        workflow_profile = _find_workflow_profile(snakefile)

        if profile or workflow_profile:
            defaults = _load_profile_defaults(profile, workflow_profile)
            defaults = {k.replace("-", "_"): v for k, v in defaults.items()}
            extra.setdefault("default_map", {}).update(defaults)

        return super().make_context(info_name, args, parent=parent, **extra)

    def parse_args(self, ctx, args):
        # Let Click parse the group's own options.
        rest = super().parse_args(ctx, args)

        # If the first remaining arg is not a registered subcommand,
        # move everything into ctx.args so that invoke_without_command
        # fires and we fall through to legacy_main.
        if ctx.protected_args:
            cmd_name = ctx.protected_args[0]
            if cmd_name not in self.commands:
                ctx.args = ctx.protected_args + ctx.args
                ctx.protected_args = []

        return rest


@click.group(
    cls=ProfileAwareGroup,
    invoke_without_command=True,
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    ),
    add_help_option=False,
)
@click.pass_context
def cli(ctx):
    """Snakemake workflow engine."""
    if ctx.invoked_subcommand is None:
        legacy_main(sys.argv[1:])

@cli.command()
@click.pass_context
def help(ctx):
    """Show available subcommands."""
    click.echo(ctx.parent.get_help())


cli.add_command(lint)
cli.add_command(unlock)
cli.add_command(clean)
cli.add_command(dagviz)
cli.add_command(utils)
cli.add_command(info)
cli.add_command(software)
cli.add_command(new)
