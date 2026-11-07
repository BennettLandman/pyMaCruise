#!/usr/bin/env python3
"""Regenerate the complete standalone command reference from argparse."""
from __future__ import annotations
import argparse
from pathlib import Path
from pymacruise.cli import build_parser


def render():
    parser = build_parser()
    commands = next(action for action in parser._actions if isinstance(action, argparse._SubParsersAction)).choices
    rows = ["# Complete command reference", "", 
        "Generated from the current command parser with `python scripts/generate_cli_reference.py`. All commands below run independently; supply the required inputs shown in their usage. Use new output directories, retain provenance, and keep images on the physical grids required by [the detailed I/O guide](cli-help.md).",
        "", "For full PLUS testing and a stage-by-stage recipe, start with [the PLUS walkthrough](testing/full-plus-walkthrough.md) and [the individual tools guide](testing/individual-tools.md). Ordinary argparse defaults are listed below; the PLUS testing kit deliberately overrides initial atlas alignment to 50 iterations/four levels.", ""]
    for name, command in commands.items():
        rows += [f"## `{name}`", "", "```text", command.format_help().rstrip(), "```", ""]
        effective_defaults = {}
        if name == "base-run":
            from pymacruise.base_workflow import BaseRunOptions
            effective_defaults = dict(zip(
                ("tgdm_initial_iterations", "tgdm_inner_iterations", "tgdm_central_iterations"),
                BaseRunOptions().effective_tgdm_iterations(), strict=True,
            ))
        defaults = [(action.option_strings[0], effective_defaults.get(action.dest, action.default)) for action in command._actions
                    if action.option_strings and action.dest != 'help' and not action.required
                    and action.default is not argparse.SUPPRESS]
        if defaults:
            rows += ["Defaults: " + "; ".join(f"`{key}` = `{value}`" for key, value in defaults) + ".", ""]
    from pymacruise.lesion_output_cli import build_parser as output_parser
    rows += ["## `pymacruise-lesion-output`", "",
             "Also available as `python -m pymacruise.lesion_output_cli`. This writes already-computed arrays; it does not run segmentation.",
             "", "```text", output_parser().format_help().rstrip(), "```", ""]
    return "\n".join(rows)


if __name__ == '__main__':
    target = Path(__file__).resolve().parents[1] / 'docs/tools-reference.md'
    target.write_text(render())
    print(target)
