"""Simple command-line interface for computing TRCSS from CSV/BED input.

Usage:
    trcss-compute INPUT.csv --fd FD_COL --txdir TXDIR_COL [-o OUTPUT.csv]

If the input is a BED-like file (chrom, start, end, name, ...) with FD and
TxDir columns, this will append a ``trcss`` column and write the result.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .core import compute_trcss_dataframe


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="trcss-compute",
        description=(
            "Compute TRCSS (Transcription-Replication Context Score) for a "
            "CSV/BED table containing fork-directionality (FD) and "
            "transcription-direction (TxDir) columns."
        ),
    )
    p.add_argument(
        "input",
        type=Path,
        help="Path to input CSV/TSV/BED file. Columns are auto-detected.",
    )
    p.add_argument(
        "--fd",
        default="fd",
        help="Name of the fork-directionality column (default: fd).",
    )
    p.add_argument(
        "--txdir",
        default="txdir",
        help="Name of the transcription-direction column (default: txdir).",
    )
    p.add_argument(
        "--out-col",
        default="trcss",
        help="Name of the output TRCSS column (default: trcss).",
    )
    p.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output CSV path. Defaults to stdout if not specified.",
    )
    p.add_argument(
        "--sep",
        default="auto",
        help="CSV separator. Use 'auto' (default), ',', '\\t', or ' '.",
    )
    return p


def _detect_sep(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in (".tsv", ".bed"):
        return "\t"
    return ","


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    sep = _detect_sep(args.input) if args.sep == "auto" else args.sep
    try:
        df = pd.read_csv(args.input, sep=sep)
    except FileNotFoundError:
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001
        print(f"Error reading input: {exc}", file=sys.stderr)
        return 2

    for col in (args.fd, args.txdir):
        if col not in df.columns:
            print(
                f"Error: required column '{col}' not present. "
                f"Found: {list(df.columns)}",
                file=sys.stderr,
            )
            return 3

    compute_trcss_dataframe(df, fd_col=args.fd, txdir_col=args.txdir,
                            out_col=args.out_col, inplace=True)

    if args.output is None:
        df.to_csv(sys.stdout, index=False)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
