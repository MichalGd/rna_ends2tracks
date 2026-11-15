"""Portable, explicitly separated extended-exon DE/ORA and shared GSEA report."""
from __future__ import annotations

import csv
import html
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path

from .expression_methods import NEW_METHOD, expression_methods, expression_rows
from .receipts import sha256


def rows(path: Path):
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def enrichment_summary(index: dict, method: str, fdr: float, limit=10):
    grouped = defaultdict(list)
    for row in rows(Path(index["gsea_file" if method == "GSEA" else "ora_file"])):
        padj = number(row.get("padj"))
        if padj is None or padj > fdr:
            continue
        if method == "GSEA":
            nes = number(row.get("NES"))
            if nes is None or nes == 0:
                continue
            direction = "positive" if nes > 0 else "negative"
            condition = index.get("numerator", "numerator") if nes > 0 else index.get("denominator", "denominator")
            label = f"{direction} NES; favored: {condition}"
        else:
            direction = row.get("query", "")
            label = {"upregulated": "higher in " + index.get("numerator", "numerator"),
                     "downregulated": "higher in " + index.get("denominator", "denominator"),
                     "any_dge": "any DE (nondirectional)"}.get(direction, direction)
        grouped[(row.get("database", ""), direction)].append({
            "database": row.get("database", ""), "direction": label, "term_id": row.get("term_id", ""),
            "term_name": row.get("term_name", ""), "padj": row.get("padj", ""),
            **({"NES": row["NES"]} if method == "GSEA" else {"overlap_count": row.get("overlap_count", "")})})
    return [row for key in sorted(grouped) for row in sorted(grouped[key], key=lambda r: number(r["padj"]))[:limit]]


def render_expression_report(plan, results: Path, outdir: Path, image_assets: list, gallery, table):
    if NEW_METHOD not in expression_methods(plan.project) or not plan.project.get("modules", {}).get("gene_expression", True):
        return "", [], []
    indexes = [row for genome in sorted(plan.references) for row in expression_rows(plan.project, results, genome)
               if row["expression_method"] == NEW_METHOD]
    enrichment = [row for row in rows(outdir / "enrichment_summary" / "enrichment_index.tsv")
                  if row.get("expression_method") == NEW_METHOD]
    for index in indexes:
        data = rows(Path(index["result_file"]))
        valid = [row for row in data if row.get("valid_test", "").lower() == "true"]
        if int(index["tested_genes"]) != len(valid) or int(index["eligible_genes"]) != len(data):
            raise RuntimeError("Expression index universe disagrees with complete result table")
        if index["analysis_mode"] != "shared_gsea":
            significant = [row for row in valid if number(row.get("padj")) is not None and number(row["padj"]) <= float(index["fdr"])]
            if int(index["significant"]) != len(significant):
                raise RuntimeError("Expression index significance disagrees with complete result table")
            index["higher_numerator"] = sum(number(row["log2FoldChange"]) > 0 for row in significant)
            index["higher_denominator"] = sum(number(row["log2FoldChange"]) < 0 for row in significant)
    fragments = ["<h2 id='extended-expression'>Extended-exon gene expression</h2>"]
    markdown = ["", "## Extended-exon gene expression", ""]
    outputs, manifest = [], []

    def prose(text):
        fragments.append("<p>"+html.escape(text)+"</p>")
        markdown.extend([text, ""])

    def link(path, label):
        path = Path(path)
        if not path.is_file():
            raise RuntimeError(f"Missing expression report asset: {path}")
        relative = path.resolve().relative_to(results.resolve())
        target = outdir / "assets" / "expression" / relative
        if not target.parent.resolve().is_relative_to(outdir.resolve()):
            raise RuntimeError("Expression report asset directory escapes report")
        target.parent.mkdir(parents=True, exist_ok=True)
        digest = sha256(path)
        if target.is_symlink():
            target.unlink()
        if not target.is_file() or sha256(target) != digest:
            shutil.copyfile(path, target)
        local = target.relative_to(outdir).as_posix()
        outputs.append(target)
        manifest.append({"source": str(path), "asset": local, "sha256": digest})
        fragments.append(f"<p><a href='{html.escape(local, quote=True)}'>{html.escape(label)}</a></p>")
        markdown.append(f"[{label}]({local})\n")

    settings = plan.project["gene_expression"]
    prose(f"Counts use original exon unions plus bounded outer 3-prime extensions (maximum {settings['extension_max_nt']} nt; "
          f"same-strand TSS buffer {settings['tss_buffer_nt']} nt). They do not require PAS discovery or exact-end/internal-priming filtering. "
          "Extensions can capture genuine unannotated ends and internal-priming noise; recovered genes require biological review. "
          "Reads overlapping multiple genes on the compatible strand are excluded; exon overlaps within one gene count once.")
    prose(f"Per-contrast edgeR filterByExpr: min.count={settings['min_count']}, min.total.count={settings['min_total_count']}, "
          "using assigned library sizes and design leverage. Both DE tests share the eligible and valid tested universe. "
          "DESeq2 poscounts size factors are estimated per contrast; plotted fold changes are unshrunken MLEs. "
          "BH adjustment is separate for each hypothesis and contrast; no combined FDR across modes is claimed.")
    for genome in sorted(plan.references):
        root = results / "05_gene_expression" / genome / NEW_METHOD
        for name in ("extension_audit.tsv", "gene_metadata.tsv", "annotation_manifest.json"):
            link(root / "annotation" / name, f"{genome}: {name}")
        for name in ("raw_counts.tsv", "featureCounts.tsv.summary"):
            link(root / name, f"{genome}: {name}")
        assignment = rows(root / "featureCounts.tsv.summary")
        summary = []
        for column in (key for key in assignment[0] if key != "Status") if assignment else []:
            counts = {row["Status"]: int(row[column]) for row in assignment}
            total = sum(counts.values())
            summary.append({"sample_id": Path(column).stem, "input_units": total,
                "assigned": counts.get("Assigned", 0), "ambiguous": counts.get("Unassigned_Ambiguity", 0),
                "no_feature": counts.get("Unassigned_NoFeatures", 0),
                "assigned_pct": round(100*counts.get("Assigned", 0)/total, 2) if total else "NA"})
        if summary:
            fragments.append(table(summary,list(summary[0])))
        link(root / "qc" / "library_sizes.tsv", f"{genome}: library depths and size factors")
        fragments.append(gallery(outdir, [(f"{genome} extended-exon PCA", root / "qc" / "pca.png"),
                                          (f"{genome} extended-exon sample distances", root / "qc" / "sample_distances.png")], image_assets))
    for mode, title, method in (("stringent", "Stringent DE and ORA", "ORA"),
                                ("relaxed_fdr", "Relaxed DE and ORA", "ORA"),
                                ("shared_gsea", "Shared ranked GSEA", "GSEA")):
        selected = [row for row in indexes if row["analysis_mode"] == mode]
        fragments.append(f"<h3>{title}</h3>")
        markdown.extend([f"### {title}", ""])
        if not selected:
            prose("Disabled or no eligible contrasts.")
            continue
        if mode == "shared_gsea":
            prose("GSEA runs once on all common valid tested genes ranked by the ordinary Wald statistic. "
                  "It is not filtered to significant genes and is shared by both DE modes. Positive NES favors the numerator; negative NES favors the denominator.")
        else:
            prose((f"Formal minimum-effect test: abs(log2 fold change) > {settings['lfc_threshold']}; " if mode == "stringent" else "Ordinary zero-null test: log2 fold change != 0; ")+
                  f"adjusted p-value <= {plan.project['reporting']['fdr']}. ORA uses this mode's significant genes, separated by direction; any_dge is nondirectional.")
        fields = ["genome", "contrast_id", "numerator", "denominator", "eligible_genes", "tested_genes"]
        if mode != "shared_gsea":
            fields += ["significant", "higher_numerator", "higher_denominator", "hypothesis", "fdr"]
        fragments.append(table(selected, fields))
        markdown.extend(["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"]*len(fields)) + " |"])
        markdown.extend("| " + " | ".join(str(row.get(key,"")) for key in fields) + " |" for row in selected)
        for index in selected:
            ident = f"{index['genome']} {index['contrast_id']}"
            fragments.append(f"<h4>{html.escape(ident)}</h4>")
            markdown.extend(["", f"#### {ident}", ""])
            link(index["result_file"], f"{ident} {mode}: complete result table")
            common = results / "05_gene_expression" / index["genome"] / NEW_METHOD / "contrasts" / index["contrast_id"]
            if mode == "shared_gsea":
                for name in ("eligibility.tsv", "test_universe.tsv", "filter_parameters.tsv", "sessionInfo.txt"):
                    link(common / name, f"{ident}: {name}")
            else:
                fragments.append(gallery(outdir, [(f"{ident} {mode} {label}", Path(index[field]))
                    for field,label in (("ma_png","MA"),("volcano_png","volcano")) if index.get(field)], image_assets))
            jobs = [row for row in enrichment if row["genome"] == index["genome"] and row["contrast_id"] == index["contrast_id"] and row["analysis_mode"] == mode]
            if not jobs:
                prose(f"{method} disabled or unavailable for this mode.")
            for job in jobs:
                prose(f"{method} background: {job['background_genes']} valid tested genes before database mapping; foreground: {job['foreground_genes']}. "
                      f"Enrichment adjusted p <= {plan.project['enrichment']['padj']}. BH scope: " +
                      ("each query and database within this contrast/mode." if method == "ORA" else "each database within this contrast."))
                terms = enrichment_summary(job, method, float(plan.project["enrichment"]["padj"]))
                if terms:
                    fragments.append(table(terms,list(terms[0])))
                    fields = list(terms[0])
                    markdown += ["| "+" | ".join(fields)+" |", "| "+" | ".join(["---"]*len(fields))+" |"]
                    markdown += ["| "+" | ".join(str(row[key]).replace("|","/") for key in fields)+" |" for row in terms]
                prose("Selected summaries retain up to ten significant terms per database and direction. A missing direction means no term passed the threshold; consult complete tables for all tested terms.")
                for field in (("ora_file",) if method == "ORA" else ("gsea_file",)) + ("prepared_gene_table", "mapping_audit", "provenance_file"):
                    link(job[field], f"{ident} {mode}: {field}")
                mapping_genes = Path(job["mapping_audit"]).with_name("mapping_genes.tsv")
                if mapping_genes.is_file():
                    link(mapping_genes,f"{ident} {mode}: per-gene database coverage")
                plots = [(f"{ident} {mode} {plot['method']} {plot['query']} {plot['database']} {plot['plot_type']}",Path(plot["png"]))
                         for plot in rows(Path(job["plot_index"])) if plot.get("method") == method and plot.get("png")]
                fragments.append(gallery(outdir, plots, image_assets))
    target = outdir / "expression_report_assets.json"
    target.write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    return "\n".join(fragments), markdown, [*dict.fromkeys(outputs), target]
