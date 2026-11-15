"""Portable figure assets, also usable as a standalone existing-report repair tool."""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


FIGURE_ASSET_LAYOUT = 1
IMAGE_SUFFIXES = {".png", ".svg", ".jpg", ".jpeg", ".gif", ".webp", ".avif"}


def _atomic_write(path: Path, contents: bytes) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_bytes(contents)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def copy_report_image(outdir: Path, source: Path) -> dict[str, object]:
    """Copy one figure into the report; content-based names avoid collisions."""
    outdir = outdir.resolve()
    source = source.resolve()
    suffix = source.suffix.lower()
    if suffix not in IMAGE_SUFFIXES:
        raise RuntimeError(f"Unsupported report image type: {source}")
    contents = source.read_bytes()
    if not contents:
        raise RuntimeError(f"Report image is empty: {source}")
    digest = hashlib.sha256(contents).hexdigest()
    directory = outdir / "assets" / "figures"
    if not directory.resolve().is_relative_to(outdir):
        raise RuntimeError(f"Report assets directory points outside the report: {directory}")
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / (digest + suffix)
    # Replace symlinks rather than following a link outside this report.
    if target.is_symlink() or not target.is_file() or target.read_bytes() != contents:
        _atomic_write(target, contents)
    return {
        "source": str(source), "asset": target.relative_to(outdir).as_posix(),
        "sha256": digest, "size_bytes": len(contents),
    }


def write_image_manifest(path: Path, records: list[dict[str, object]]) -> None:
    payload = {
        "schema_version": FIGURE_ASSET_LAYOUT,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "images": records,
    }
    _atomic_write(path, (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8"))


class _ReportTags(HTMLParser):
    def __init__(self, contents: str):
        super().__init__(convert_charrefs=False)
        self.tags: list[tuple[str, list[tuple[str, str | None]], int, int, bool]] = []
        self.offsets = [0]
        for line in contents.splitlines(keepends=True):
            self.offsets.append(self.offsets[-1] + len(line))
        self.feed(contents)
        self.close()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag not in {"img", "a"}:
            return
        raw = self.get_starttag_text()
        line, column = self.getpos()
        start = self.offsets[line - 1] + column
        self.tags.append((tag, attrs, start, start + len(raw), raw.rstrip().endswith("/>")))

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.handle_starttag(tag, attrs)


def _local_image_path(report: Path, value: str) -> Path | None:
    url = urlsplit(value)
    if url.scheme not in {"", "file"} or url.netloc not in {"", "localhost"}:
        return None  # Embedded data and remote URLs are not fetched or modified.
    decoded = unquote(url.path)
    if os.name == "nt" and url.scheme == "file" and decoded.startswith("/"):
        decoded = decoded[1:]
    path = Path(decoded)
    return (path if path.is_absolute() else report.parent / path).resolve()


def repair_report_figures(report: Path, output: Path | None = None) -> dict[str, object]:
    """Write report.portable.html without rerunning statistics or changing receipts.

    Run where the referenced original figures still exist (normally the server).
    Local image references must resolve inside the report's results directory.
    All tables, scripts and non-figure links are retained byte-for-byte as text.
    """
    report = report.expanduser().resolve()
    if output:
        output = output.expanduser()
        output = (output if output.is_absolute() else report.parent / output).absolute()
    else:
        output = report.with_name("report.portable.html")
    if not report.is_file():
        raise RuntimeError(f"Report does not exist: {report}")
    if output.parent.resolve() != report.parent or output.resolve() == report:
        raise RuntimeError("Choose a separate HTML output file in the same report directory.")
    contents = report.read_bytes().decode("utf-8-sig")
    parser = _ReportTags(contents)
    records: list[dict[str, object]] = []
    references: dict[Path, str] = {}
    local_images: dict[Path, str] = {}
    remote_images = 0
    embedded_images = 0
    results = report.parent.parent
    # Validate every source before publishing a new HTML document.
    for tag, attrs, _start, _end, _closed in parser.tags:
        if tag != "img":
            continue
        attributes = dict(attrs)
        value = attributes.get("src")
        if not value:
            continue
        if attributes.get("srcset"):
            raise RuntimeError("Repair of srcset images is not supported; original HTML was left unchanged.")
        source = _local_image_path(report, value)
        if source is None:
            if value.lower().startswith("data:"):
                embedded_images += 1
            else:
                remote_images += 1
            continue
        if not source.is_relative_to(results):
            raise RuntimeError(f"Image points outside this results directory: {source}")
        if not source.is_file():
            raise RuntimeError(
                f"Missing figure: {source}\nRun this repair on the server where the complete results tree exists."
            )
        local_images[source] = attributes.get("alt") or source.name
    for source, label in local_images.items():
        record = copy_report_image(report.parent, source)
        record["label"] = label
        records.append(record)
        references[source] = str(record["asset"])

    replacements = []
    for tag, attrs, start, end, closed in parser.tags:
        field = "src" if tag == "img" else "href"
        value = dict(attrs).get(field)
        if not value:
            continue
        source = _local_image_path(report, value)
        if source not in references:
            continue
        new_value = references[source]
        fragment = urlsplit(value).fragment
        if fragment:
            new_value += "#" + fragment
        rewritten = "<" + tag + "".join(
            " " + name + ("" if attr_value is None else
                          '="' + html.escape(new_value if name == field else attr_value, quote=True) + '"')
            for name, attr_value in attrs
        ) + ("/>" if closed else ">")
        replacements.append((start, end, rewritten))
    for start, end, replacement in reversed(replacements):
        contents = contents[:start] + replacement + contents[end:]
    manifest = output.with_suffix(".assets.json")
    write_image_manifest(manifest, records)
    _atomic_write(output, contents.encode("utf-8"))
    return {
        "html": str(output), "manifest": str(manifest), "local_images": len(local_images),
        "unique_assets": len(set(references.values())), "embedded_images": embedded_images,
        "remote_images": remote_images,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect figures for an existing rna_ends2tracks HTML report.")
    parser.add_argument("--report", type=Path, required=True, help="Original 10_reports/report.html on the server")
    parser.add_argument("--output", type=Path, help="New HTML filename in the same 10_reports folder")
    args = parser.parse_args()
    try:
        summary = repair_report_figures(args.report, args.output)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"REPORT FIGURE REPAIR ERROR: {exc}\n")
    print(json.dumps(summary, indent=2))
    print("Download the complete 10_reports folder and open report.portable.html.")
    if summary["remote_images"]:
        print("Remote image URLs were retained and still require a network connection.")


if __name__ == "__main__":
    main()
