from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pysam

from .apa_a import annotate_site, load_genes
from .paths import workflow_asset
from .polyaseqtrap_adapter import (
    _load_json,
    environment_executable,
    extract_end_counts,
    r_subprocess_environment,
    reuse_exact_end_counts,
    verify_installation,
)


def _read(header: dict[str, object], name: str, start: int, flag: int) -> pysam.AlignedSegment:
    record = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
    record.query_name = name
    record.query_sequence = "CGTC" * 12 + "CG"
    record.flag = flag
    record.reference_id = 0
    record.reference_start = start
    record.mapping_quality = 60
    record.cigar = ((0, 50),)
    record.query_qualities = pysam.qualitystring_to_array("I" * 50)
    return record


def _write_coordinate_fixture(work: Path) -> tuple[Path, Path, Path]:
    fasta = work / "pilot.fa"
    fasta.write_text(">chrPilot\n" + ("CGTC" * 500) + "\n", encoding="ascii")
    pysam.faidx(str(fasta))
    gtf = work / "pilot.gtf"
    gtf.write_text(
        'chrPilot\tpilot\tgene\t101\t1000\t.\t+\t.\tgene_id "pilot_gene";\n'
        'chrPilot\tpilot\texon\t101\t200\t.\t+\t.\tgene_id "pilot_gene";\n'
        'chrPilot\tpilot\texon\t901\t1000\t.\t+\t.\tgene_id "pilot_gene";\n',
        encoding="utf-8",
    )
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrPilot", "LN": 2000}]}
    bam = work / "pilot.bam"
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        handle.write(_read(header, "plus_primary", 300, 16))
        handle.write(_read(header, "plus_duplicate_flagged", 300, 16 | 1024))
        handle.write(_read(header, "minus_primary", 700, 0))
    return fasta, gtf, bam


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _deepip_truth(deepip_script: Path, model: Path, work: Path) -> dict[str, object]:
    source = deepip_script.parent / "test_data" / "train_mini.fa"
    if not source.is_file():
        raise RuntimeError(f"Pinned DeepIP labeled truth sequences are unavailable: {source}")
    output = work / "deepip_truth_predictions.csv"
    log = work / "deepip_truth.log"
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [sys.executable, str(deepip_script),
             "-testSeq", str(source), "-trainedModel", str(model), "-outputFile", str(output)],
            stdout=handle, stderr=subprocess.STDOUT, check=True,
        )
    with output.open(encoding="utf-8-sig", newline="") as handle:
        predictions = list(csv.DictReader(handle))
    truth = [(int(row["true_label"]), int(row["predict_label"])) for row in predictions]
    correct_artifacts = sum(observed == 0 and predicted == 0 for observed, predicted in truth)
    correct_pas = sum(observed == 1 and predicted == 1 for observed, predicted in truth)
    accuracy = sum(observed == predicted for observed, predicted in truth) / max(1, len(truth))
    return {
        "deepip_truth_sequences": len(truth),
        "deepip_correct_artifacts": correct_artifacts,
        "deepip_correct_pas": correct_pas,
        "deepip_accuracy": accuracy,
        "deepip_artifact_rejected": correct_artifacts > 0 and accuracy >= 0.75,
        "deepip_true_pas_retained": correct_pas > 0 and accuracy >= 0.75,
    }


def execute(args: argparse.Namespace) -> int:
    output = Path(args.output).resolve()
    work = output.parent / (output.stem + ".work")
    work.mkdir(parents=True, exist_ok=True)
    installation = _load_json(Path(args.installation_manifest), "APA-B installation manifest")
    model, deepip_script, environment_hash = verify_installation(installation, "human")
    fasta, gtf, bam = _write_coordinate_fixture(work)
    ends = work / "exact_end_counts.tsv"
    extraction = extract_end_counts(bam, ends)
    observed = {(row["chrom"], int(row["position"]), row["strand"]): int(row["count"])
                for row in _rows(ends)}
    expected = {("chrPilot", 349, "+"): 2, ("chrPilot", 700, "-"): 1}

    # Prove that the optimized production path reconstructs precisely the same
    # raw endpoint universe from receipt-validated C1 + C1S, never from C2.
    reusable = work / "receipt_validated_exact_ends"
    reusable.mkdir(exist_ok=True)
    c1 = reusable / "C1_exact_ends.tsv.gz"
    c1s = reusable / "C1S_uncertain_ends.tsv.gz"
    with gzip.open(c1, "wt", encoding="utf-8") as handle:
        handle.write("chrom\tstart\tend\tstrand\tcount\nchrPilot\t349\t350\t+\t2\n")
    with gzip.open(c1s, "wt", encoding="utf-8") as handle:
        handle.write("chrom\tstart\tend\tstrand\tcount\nchrPilot\t700\t701\t-\t1\n")
    end_audit = reusable / "end_audit.json"
    end_audit.write_text(json.dumps({
        "sample_id": "pilot", "C0": 3, "C1": 2, "C1S": 1,
        "duplicate_flagged": 1, "duplicate_flagged_C0": 1,
    }), encoding="utf-8")
    receipt_dir = reusable / ".receipt"
    receipt_dir.mkdir(exist_ok=True)
    receipt_outputs = []
    for path in (c1, c1s, end_audit):
        receipt_outputs.append({
            "path": str(path.resolve()), "size": path.stat().st_size,
            "validation": "sha256", "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    receipt = receipt_dir / "run_receipt.json"
    receipt.write_text(json.dumps({
        "module": "exact_ends_sample", "exit_status": 0, "outputs": receipt_outputs,
    }), encoding="utf-8")
    reused_ends = work / "reused_exact_end_counts.tsv"
    reuse_audit = reuse_exact_end_counts({
        "sample_id": "pilot", "c1": str(c1), "c1s": str(c1s),
        "exact_end_audit": str(end_audit), "exact_end_receipt": str(receipt),
    }, reused_ends)
    reused = {(row["chrom"], int(row["position"]), row["strand"]): int(row["count"])
              for row in _rows(reused_ends)}

    clustered = work / "polyaseqtrap_candidates.tsv"
    cluster_audit = work / "polyaseqtrap_audit.tsv"
    r_script = workflow_asset("scripts/R/polyaseqtrap_quantseq_rev.R")
    with (work / "polyaseqtrap.log").open("w", encoding="utf-8") as handle:
        subprocess.run(
            [environment_executable("Rscript"), str(r_script), "--ends", str(ends), "--fasta", str(fasta),
             "--output", str(clustered), "--audit", str(cluster_audit), "--cluster-gap", "24"],
            env=r_subprocess_environment(), stdout=handle, stderr=subprocess.STDOUT, check=True,
        )
    clustered_rows = _rows(clustered)
    clustered_counts = sum(int(row["count"]) for row in clustered_rows)
    genes, bins = load_genes(str(gtf))
    gene_id, feature = annotate_site("chrPilot", "+", 349, genes, bins)
    evidence = {
        "coordinate_and_strand": observed == expected,
        "c1_c1s_reuse_equivalent": reused == observed == expected
        and reuse_audit["records_written"] == extraction["records_written"],
        "record_count_conserved": (
            extraction["eligible_records"] == extraction["records_written"] == clustered_counts == 3
        ),
        "duplicate_flagged_records_retained": extraction["duplicate_flagged_records_retained"] == 1,
        "intragenic_site_retained": gene_id == "pilot_gene" and feature == "intron",
        "polyaseqtrap_cluster_count": len(clustered_rows),
        "environment_sha256": environment_hash,
        **_deepip_truth(deepip_script, model, work),
    }
    required = (
        "coordinate_and_strand", "c1_c1s_reuse_equivalent", "record_count_conserved",
        "duplicate_flagged_records_retained",
        "deepip_artifact_rejected", "deepip_true_pas_retained", "intragenic_site_retained",
    )
    evidence["status"] = "PASS" if all(evidence[key] is True for key in required) else "FAIL"
    output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if evidence["status"] != "PASS":
        raise RuntimeError(f"APA-B synthetic pilot failed; inspect {output}")
    print(f"APA-B synthetic pilot: PASS ({output})")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the pinned APA-B synthetic coordinate and DeepIP pilot")
    parser.add_argument("--installation-manifest", required=True)
    parser.add_argument("--output", required=True)
    raise SystemExit(execute(parser.parse_args()))


if __name__ == "__main__":
    main()
