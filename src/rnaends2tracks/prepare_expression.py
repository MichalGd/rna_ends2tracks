"""Prepare a separate BAM-input expression project without starting analysis."""
from __future__ import annotations

import argparse
import csv
import json
import shlex
import shutil
from pathlib import Path

from .conf import read_conf, project_from_conf
from .config import build_conf_plan
from .receipts import sha256


def prepare_project(source_config: Path, directory: Path, project_id: str, jobs: int = 8) -> Path:
    if jobs < 1:
        raise ValueError("jobs must be positive")
    source_config = source_config.resolve()
    directory = directory.resolve()
    source_project, source_sheet = project_from_conf(source_config)
    source_results = Path(source_project["output_dir"]).resolve()
    bam_root = Path(source_project["alignment_source_dir"]) if source_project.get("input_mode") == "bam" else source_results / "02_alignment"
    if directory == source_results or directory.is_relative_to(source_results) or source_results.is_relative_to(directory):
        raise ValueError("Choose a new project outside the original result tree")
    if directory.exists() and any(directory.iterdir()):
        raise ValueError(f"New project directory must be empty: {directory}")
    values = {key: value for key, value in read_conf(source_config).items() if not key.startswith("_")}
    with Path(source_sheet).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        samples = list(reader)
    if not fields or not samples:
        raise ValueError("Original samplesheet is empty")
    for row in samples:
        for key in ("fastq_r1", "fastq_r2"):
            if row.get(key):
                path = Path(row[key])
                row[key] = str((path if path.is_absolute() else Path(source_sheet).parent / path).resolve())
    for prefix, genome in (("MM39", "GRCm39"), ("HG38", "GRCh38")):
        reference = source_project["references"][genome]
        for suffix, key in (("GTF","gtf"),("CHROM_SIZES","chrom_sizes"),("FASTA","fasta"),("STAR_INDEX","star_index"),("PAS_ATLAS","pas_atlas"),("RSEQC_BED","rseqc_bed")):
            values[f"{prefix}_{suffix}"] = reference.get(key, "")
    values.update(PROJECT_ID=project_id, SAMPLESHEET="samplesheet.csv", OUTPUT_DIR=str(directory / "results"),
        TMP_DIR=str(directory / "tmp"), INPUT_MODE="bam", ALIGNMENT_SOURCE_DIR=str(bam_root),
        DGE_METHODS="extended_exons_v1", DGE_PRIMARY_METHOD="extended_exons_v1", RUN_GENE_EXPRESSION="true",
        GEX_EXTENSION_MAX_NT="6000", GEX_TSS_BUFFER_NT="2000", GEX_FILTER_METHOD="filterByExpr", GEX_MIN_COUNT="5", GEX_MIN_TOTAL_COUNT="15",
        RUN_DGE_STRINGENT="true", RUN_DGE_RELAXED="true", DGE_STRINGENT_LFC_THRESHOLD="1.0", FDR="0.05",
        RUN_DGE_ENRICHMENT="true", ENRICHMENT_ORA="true", ENRICHMENT_GSEA="true", RUN_APA_ENRICHMENT="false",
        RUN_APA_A_MCELL2019="false", RUN_APA_A2="false", RUN_APA_B="false", RUN_RSEQC="false", RUN_TRACKS="false",
        CLEANUP_INTERMEDIATES="false", DGE_CONTRAST_PARALLEL_JOBS=str(jobs), ENRICHMENT_PARALLEL_JOBS=str(jobs),
        MAX_TOTAL_THREADS=str(max(int(values["MAX_TOTAL_THREADS"]),jobs)), PARALLEL_DOWNSTREAM_MODULES="false")
    def quoted(value):
        value = str(value)
        quote = '"' if '"' not in value else "'"
        if quote in value or "\n" in value or "\r" in value:
            raise ValueError("Configuration value cannot be represented safely")
        return quote+value+quote
    text = "# Prepared by alpha13 from an existing config; analysis has not started.\n" + "".join(f"{key}={quoted(value)}\n" for key,value in values.items())
    directory.mkdir(parents=True,exist_ok=True)
    sheet = directory / "samplesheet.csv"
    with sheet.open("x",encoding="utf-8",newline="") as handle:
        writer = csv.DictWriter(handle,fieldnames=fields)
        writer.writeheader(); writer.writerows(samples)
    config = directory / "config.conf"
    with config.open("x",encoding="utf-8") as handle:
        handle.write(text)
    plan = build_conf_plan(config,check_inputs=False)
    prior_contrasts = source_results / "00_metadata" / "contrasts.tsv"
    if prior_contrasts.is_file():
        with prior_contrasts.open(encoding="utf-8-sig",newline="") as handle:
            previous = list(csv.DictReader(handle,delimiter="\t"))
        keys = ("genome","contrast_id","numerator","denominator","resolved_design")
        old = {tuple(row.get(key,"") for key in keys) for row in previous}
        new = {tuple(str(row.get(key,"")) for key in keys) for row in plan.contrasts}
        if old != new:
            raise ValueError("Prepared contrasts/design differ from original validated contrasts; review the new project before launch")
    snapshot = directory / "source_metadata"
    snapshot.mkdir(exist_ok=True)
    sources = {"original_config.conf": source_config,"original_samplesheet.csv": Path(source_sheet)}
    for name in ("contrasts.tsv","validated_samples.tsv","resolved_config.json"):
        path = source_results / "00_metadata" / name
        if path.is_file(): sources[name] = path
    audit = []
    for name,path in sources.items():
        shutil.copyfile(path,snapshot / name)
        audit.append({"source":str(path),"snapshot":name,"sha256":sha256(path)})
    (snapshot / "source_manifest.json").write_text(json.dumps(audit,indent=2)+"\n",encoding="utf-8")
    (directory / "README.md").write_text(
        "# Prepared expression project\n\nNo analysis has started. Original alignments and outputs were not modified.\n\n"
        "The source design/pairing and condition order were retained. Validate actual BAMs in the alpha13 Linux environment first:\n\n"
        f"```bash\nrna-ends2tracks --config {shlex.quote(str(config))} --stop-after validate\n"
        f"rna-ends2tracks --config {shlex.quote(str(config))}\n```\n",encoding="utf-8")
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-config",type=Path,required=True)
    parser.add_argument("--new-project",type=Path,required=True)
    parser.add_argument("--project-id",required=True)
    parser.add_argument("--jobs",type=int,default=8)
    args = parser.parse_args()
    try:
        print(prepare_project(args.source_config,args.new_project,args.project_id,args.jobs))
    except (ValueError,OSError) as exc:
        parser.exit(2,f"ERROR: {exc}\nNo workflow was launched.\n")


if __name__ == "__main__":
    main()
