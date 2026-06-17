#!/usr/bin/env python3
"""Small Python wrapper around the wavecar_compressor executable."""

from __future__ import annotations

import argparse
import bz2
import csv
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_BINARY = Path(__file__).resolve().parent / "build" / "wavecar_compressor"


@dataclass(frozen=True)
class AnalysisPoint:
    cutoff_fraction: float
    encut: float
    kept_coefficients: int
    total_coefficients: int
    coefficient_fraction: float
    retained_norm_fraction: float
    lost_norm_fraction: float
    min_band_retained_norm_fraction: float
    estimated_record_length: int
    estimated_file_size: int
    estimated_size_ratio: float


def _binary_path(binary: str | os.PathLike[str] | None = None) -> Path:
    if binary is not None:
        return Path(binary)
    env_binary = os.environ.get("WAVECAR_COMPRESSOR_BIN")
    if env_binary:
        return Path(env_binary)
    return DEFAULT_BINARY


def _run(
    args: list[str],
    binary: str | os.PathLike[str] | None = None,
) -> subprocess.CompletedProcess[str]:
    command = [str(_binary_path(binary)), *args]
    return subprocess.run(command, check=True, capture_output=True, text=True)


def bzip2_compress(
    input_path: str | os.PathLike[str],
    output_path: str | os.PathLike[str] | None = None,
    compresslevel: int = 9,
    chunk_size: int = 1024 * 1024,
) -> Path:
    source = Path(input_path)
    destination = Path(output_path) if output_path is not None else source.with_name(
        source.name + ".bz2"
    )

    with source.open("rb") as input_file:
        with bz2.open(destination, "wb", compresslevel=compresslevel) as output_file:
            while True:
                chunk = input_file.read(chunk_size)
                if not chunk:
                    break
                output_file.write(chunk)

    return destination


def analyze(
    input_path: str | os.PathLike[str],
    fractions: Iterable[float] | None = None,
    binary: str | os.PathLike[str] | None = None,
) -> list[AnalysisPoint]:
    args = ["analyze", str(input_path)]
    if fractions is not None:
        args.extend(str(fraction) for fraction in fractions)

    completed = _run(args, binary)
    rows = csv.DictReader(completed.stdout.splitlines())
    return [
        AnalysisPoint(
            cutoff_fraction=float(row["cutoff_fraction"]),
            encut=float(row["encut"]),
            kept_coefficients=int(row["kept_coefficients"]),
            total_coefficients=int(row["total_coefficients"]),
            coefficient_fraction=float(row["coefficient_fraction"]),
            retained_norm_fraction=float(row["retained_norm_fraction"]),
            lost_norm_fraction=float(row["lost_norm_fraction"]),
            min_band_retained_norm_fraction=float(row["min_band_retained_norm_fraction"]),
            estimated_record_length=int(row["estimated_record_length"]),
            estimated_file_size=int(row["estimated_file_size"]),
            estimated_size_ratio=float(row["estimated_size_ratio"]),
        )
        for row in rows
    ]


def compress(
    input_path: str | os.PathLike[str],
    output_path: str | os.PathLike[str],
    cutoff_fraction: float,
    zero_small: bool = False,
    binary: str | os.PathLike[str] | None = None,
    *,
    create_bz2: bool = False,
    bz2_output: str | os.PathLike[str] | None = None,
) -> str:
    args = [
        "compress",
        str(input_path),
        str(output_path),
        str(cutoff_fraction),
        "true" if zero_small else "false",
    ]
    stdout = _run(args, binary).stdout
    if create_bz2:
        compressed_path = bzip2_compress(output_path, bz2_output)
        stdout += f"Bzip2 file written: {compressed_path}\n"
    return stdout


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary",
        default=None,
        help="Path to the wavecar_compressor executable",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze_parser = subparsers.add_parser("analyze")
    analyze_parser.add_argument("input")
    analyze_parser.add_argument("fractions", nargs="*", type=float)

    compress_parser = subparsers.add_parser("compress")
    compress_parser.add_argument("input")
    compress_parser.add_argument("output")
    compress_parser.add_argument("cutoff_fraction", type=float)
    compress_parser.add_argument("--zero-small", action="store_true")
    compress_parser.add_argument(
        "--bz2",
        action="store_true",
        help="Also write a bzip2-compressed copy of the output WAVECAR",
    )
    compress_parser.add_argument(
        "--bz2-output",
        default=None,
        help="Path for the bzip2 output; defaults to <output>.bz2",
    )

    return parser


def main() -> int:
    args = _build_parser().parse_args()
    if args.command == "analyze":
        fractions = args.fractions if args.fractions else None
        for row in analyze(args.input, fractions, args.binary):
            print(
                f"{row.cutoff_fraction:g}\t"
                f"norm={row.retained_norm_fraction:.8f}\t"
                f"min_band={row.min_band_retained_norm_fraction:.8f}\t"
                f"size={row.estimated_size_ratio:.6f}"
            )
        return 0

    if args.command == "compress":
        print(
            compress(
                args.input,
                args.output,
                args.cutoff_fraction,
                zero_small=args.zero_small,
                binary=args.binary,
                create_bz2=args.bz2,
                bz2_output=args.bz2_output,
            ),
            end="",
        )
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
