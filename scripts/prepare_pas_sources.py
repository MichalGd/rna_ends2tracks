#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[1] / "src"
if SOURCE_ROOT.is_dir():
    sys.path.insert(0, str(SOURCE_ROOT))

from rnaends2tracks.pas_sources import prepare_sources


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize official GENCODE and PolyA_DB inputs for PAS atlas v1"
    )
    parser.add_argument("--species", choices=("human", "mouse"), required=True)
    parser.add_argument("--assembly", choices=("GRCh38", "GRCm39"), required=True)
    parser.add_argument("--annotation-release", required=True)
    parser.add_argument("--gencode-polya-gtf", type=Path, required=True)
    parser.add_argument("--polyadb-zip", type=Path, required=True)
    parser.add_argument("--polyadb-assembly", required=True)
    parser.add_argument("--liftover-chain", type=Path)
    parser.add_argument("--lift-over", default="liftOver")
    parser.add_argument("--download-date", required=True, help="YYYY-MM-DD snapshot download date")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    provenance = prepare_sources(
        species=args.species,
        assembly=args.assembly,
        annotation_release=args.annotation_release,
        gencode_polya_gtf=args.gencode_polya_gtf,
        polyadb_zip=args.polyadb_zip,
        polyadb_assembly=args.polyadb_assembly,
        output=args.output,
        download_date=args.download_date,
        liftover_chain=args.liftover_chain,
        lift_over=args.lift_over,
    )
    print(json.dumps(provenance, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
