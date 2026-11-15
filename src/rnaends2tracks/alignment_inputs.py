"""Validate and record external C0 BAMs without recreating upstream outputs."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from .expression_methods import sample_bam
from .receipts import sha256, write_receipt


def validate_reused_alignments(project, samples, references, check_files=True):
    source = Path(project["alignment_source_dir"]).resolve()
    output = Path(project["output_dir"]).resolve()
    if output == source or output.is_relative_to(source) or source.is_relative_to(output):
        raise ValueError("BAM source and new output directory must be separate, non-nested directories")
    records = []
    for sample in samples:
        bam = sample_bam(project, output, sample["sample_id"])
        record = {"sample_id": sample["sample_id"], "genome": sample["genome"],
                  "bam": str(bam), "source_receipt": str(bam.parent / ".merge_receipt" / "run_receipt.json")}
        if check_files:
            import pysam
            if not bam.is_file():
                raise ValueError(f"Missing source BAM: {bam}; alignment will not be restarted")
            try:
                pysam.quickcheck(str(bam))
            except Exception as exc:
                raise ValueError(f"BAM quickcheck failed: {bam}: {exc}") from exc
            sizes = {line.split()[0]: int(line.split()[1]) for line in
                     Path(references[sample["genome"]]["chrom_sizes"]).read_text().splitlines()
                     if line.strip() and not line.startswith("#")}
            with pysam.AlignmentFile(str(bam), "rb", require_index=True) as handle:
                header = handle.header.to_dict()
                if header.get("HD", {}).get("SO") != "coordinate":
                    raise ValueError(f"Source BAM is not marked coordinate sorted: {bam}")
                if dict(zip(handle.references, handle.lengths)) != sizes:
                    raise ValueError(f"Source BAM contigs/lengths do not match reference: {bam}")
                names = {rg["SM"] for rg in header.get("RG", []) if rg.get("SM")}
                if names and names != {sample["sample_id"]}:
                    raise ValueError(f"Source BAM sample identity differs from samplesheet: {bam}: {names}")
                inspected = 0
                for read in handle.fetch(until_eof=True):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.has_tag("NH") or read.get_tag("NH") != 1:
                        raise ValueError(f"Source BAM violates mapped primary NH=1 C0 policy: {bam}")
                    inspected += 1
                    if inspected == 10000:
                        break
                record["sampled_c0_records"] = inspected
            index = next((path for path in (Path(str(bam)+".bai"), bam.with_suffix(".bai"), Path(str(bam)+".csi")) if path.is_file()), None)
            if index is None:
                raise ValueError(f"Missing BAM index: {bam}")
            record.update(index=str(index), size=bam.stat().st_size, mtime_ns=bam.stat().st_mtime_ns,
                          validation="quickcheck+indexed_header+sampled_C0_policy")
            receipt = Path(record["source_receipt"])
            record["source_receipt_status"] = "unavailable"
            if receipt.is_file():
                evidence = json.loads(receipt.read_text(encoding="utf-8"))
                outputs = [row for row in evidence.get("outputs", []) if Path(row.get("path", "")).name == bam.name]
                if evidence.get("exit_status") != 0 or len(outputs) != 1:
                    raise ValueError(f"Invalid source merge receipt: {receipt}")
                prior = outputs[0]
                if int(prior.get("size", -1)) != bam.stat().st_size:
                    raise ValueError(f"Source BAM size differs from merge receipt: {bam}")
                if prior.get("sha256") and sha256(bam) != prior["sha256"]:
                    raise ValueError(f"Source BAM checksum differs from merge receipt: {bam}")
                if prior.get("validation") == "size_mtime" and int(prior.get("mtime_ns", -1)) != bam.stat().st_mtime_ns:
                    raise ValueError(f"Source BAM timestamp differs from merge receipt: {bam}; review relocated provenance")
                record.update(source_receipt_status="matched", source_workflow_version=evidence.get("workflow_version", ""),
                              source_receipt_sha256=sha256(receipt))
        records.append(record)
    return records


def write_reuse_manifest(plan, results: Path):
    from .config import signature_for
    records = plan.project.get("reused_alignments", [])
    if not records:
        raise ValueError("No validated reused alignment records")
    directory = results / "provenance" / "alignment_inputs"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "reused_alignments.tsv"
    fields = list(dict.fromkeys(key for row in records for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    inputs = [Path(row[key]) for row in records for key in ("bam", "index") if row.get(key)]
    signature = signature_for(inputs, {"module": "alignment_inputs", "samples": plan.samples})
    # Only owned metadata is an output. External BAMs are never cleanup outputs.
    write_receipt("alignment_inputs", directory, signature, [path], ["rna-ends2tracks", "import-bams"])
    return path
