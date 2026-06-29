#!/usr/bin/env python3
"""Core tabular engine for the generic integrative RNA-seq + ChIP-seq pipeline."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import glob
import html
import math
import os
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


def env(name: str, default: str = "") -> str:
    return os.environ.get(name, default)


def outdir(*parts: str) -> Path:
    base = Path(env("INTEGRATION_OUTPUT_DIR", env("PROJECT_DIR", ".")))
    path = base.joinpath(*parts)
    path.mkdir(parents=True, exist_ok=True)
    return path


def norm_path(path: str) -> str:
    cleaned = path.replace("\\", "/")
    match = re.match(r"^([A-Za-z]):/(.*)$", cleaned)
    if match and os.name != "nt":
        drive, rest = match.groups()
        return f"/mnt/{drive.lower()}/{rest}"
    return cleaned


def path_exists(path: str) -> bool:
    return bool(path) and Path(norm_path(path)).is_file() and Path(norm_path(path)).stat().st_size > 0


def open_text(path: str):
    cleaned = norm_path(path)
    if cleaned.endswith(".gz"):
        return gzip.open(cleaned, "rt", encoding="utf-8-sig", newline="")
    return open(cleaned, "r", encoding="utf-8-sig", newline="")


def detect_delimiter(path: str) -> str:
    with open_text(path) as handle:
        first = handle.readline()
    return "\t" if first.count("\t") >= first.count(",") else ","


def read_table(path: str) -> tuple[list[str], list[dict[str, str]]]:
    if not path_exists(path):
        return [], []
    sep = detect_delimiter(path)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        header = list(reader.fieldnames or [])
        rows = [{k: (v if v is not None else "") for k, v in row.items()} for row in reader]
    return header, rows


def write_table(path: Path, rows: list[dict[str, object]], header: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def first_col(header: list[str], candidates: list[str]) -> str:
    low = {h.lower(): h for h in header}
    for name in candidates:
        if name.lower() in low:
            return low[name.lower()]
    return ""


def as_float(value: object, default: float = 0.0) -> float:
    try:
        if value is None or str(value).strip() in {"", "NA", "NaN", "nan"}:
            return default
        x = float(value)
        if math.isnan(x) or math.isinf(x):
            return default
        return x
    except (TypeError, ValueError):
        return default


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()) or "unknown"


def glob_existing(pattern: str) -> list[str]:
    if not pattern:
        return []
    return sorted(p for p in glob.glob(norm_path(pattern)) if path_exists(p))


def read_matrix_gene_ids(path: str) -> list[str]:
    header, rows = read_table(path)
    if not header:
        return []
    gene_col = first_col(header, [env("GENE_ID_COLUMN", "gene_id"), "gene", "id", "Geneid"])
    if not gene_col:
        gene_col = header[0]
    return [r.get(gene_col, "") for r in rows if r.get(gene_col, "")]


def read_set_file(path: str) -> set[str]:
    if not path_exists(path):
        return set()
    genes: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            value = line.strip().split("\t")[0].split(",")[0]
            if value and not value.lower().startswith("gene"):
                genes.add(value)
    return genes


def parse_gtf_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            bits = item.split(None, 1)
            if len(bits) != 2:
                continue
            key, value = bits
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def parse_annotation_genes(path: str) -> dict[str, dict[str, str]]:
    genes: dict[str, dict[str, str]] = {}
    if not path_exists(path):
        return genes
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _phase, raw_attrs = parts[:9]
            if feature.lower() not in {"gene", "mrna", "transcript"}:
                continue
            attrs = parse_gtf_attributes(raw_attrs)
            gene_id = attrs.get("gene_id") or attrs.get("ID") or attrs.get("Name")
            if not gene_id:
                continue
            gene_id = gene_id.replace("gene:", "")
            existing = genes.get(gene_id)
            row = {
                "gene_id": gene_id,
                "gene_name": attrs.get("gene_name") or attrs.get("Name") or gene_id,
                "chromosome": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "biotype": attrs.get("gene_biotype") or attrs.get("biotype") or attrs.get("gene_type") or "",
            }
            if existing:
                existing["start"] = str(min(as_int(existing.get("start")), as_int(start)))
                existing["end"] = str(max(as_int(existing.get("end")), as_int(end)))
                for key, value in row.items():
                    existing[key] = existing.get(key) or value
            else:
                genes[gene_id] = row
    return genes


def metadata_groups(path: str) -> dict[str, str]:
    header, rows = read_table(path)
    if not header:
        return {}
    sample_col = first_col(header, [env("SAMPLE_ID_COLUMN", "sample_id"), "sample", "file_prefix"])
    if not sample_col:
        return {}
    group_cols = [c.strip() for c in env("GROUP_COLUMNS", "").split(",") if c.strip()]
    groups = {}
    for row in rows:
        sid = row.get(sample_col, "")
        if not sid:
            continue
        values = [row.get(c, "") for c in group_cols if row.get(c, "")]
        groups[sid] = " | ".join(values) if values else sid
    return groups


def sample_group(sample: str, groups: dict[str, str]) -> str:
    if sample in groups:
        return groups[sample]
    cleaned = re.sub(r"[_-]R[0-9]+$", "", sample)
    cleaned = re.sub(r"[_-]rep[0-9]+$", "", cleaned, flags=re.IGNORECASE)
    return cleaned


def normalize_deg_rows() -> list[dict[str, str]]:
    path = env("RNA_DEG_RESULTS")
    header, rows = read_table(path)
    if not header:
        return []
    gene_col = first_col(header, [env("GENE_ID_COLUMN", "gene_id"), "gene", "id"])
    name_col = first_col(header, [env("GENE_NAME_COLUMN", "gene_name"), "gene", "symbol", "name"])
    contrast_col = first_col(header, [env("CONTRAST_ID_COLUMN", "contrast_id"), "contrast", "comparison"])
    lfc_col = first_col(header, ["log2FoldChange", "log2FC", "logFC", "lfc"])
    padj_col = first_col(header, ["padj", "FDR", "qvalue", "adj.P.Val"])
    pvalue_col = first_col(header, ["pvalue", "P.Value", "p_value"])
    out = []
    for row in rows:
        gene_id = row.get(gene_col, "") if gene_col else ""
        if not gene_id:
            continue
        lfc = as_float(row.get(lfc_col, "0"))
        padj = as_float(row.get(padj_col, "1"), 1.0)
        status = "not_significant"
        if padj <= as_float(env("DEG_PADJ_THRESHOLD", "0.05"), 0.05):
            if lfc >= as_float(env("DEG_LOG2FC_THRESHOLD", "1"), 1.0):
                status = "up"
            elif lfc <= -as_float(env("DEG_LOG2FC_THRESHOLD", "1"), 1.0):
                status = "down"
        item = dict(row)
        item.update(
            {
                "contrast_id": row.get(contrast_col, "") if contrast_col else "default",
                "gene_id": gene_id,
                "gene_name": row.get(name_col, gene_id) if name_col else gene_id,
                "log2FoldChange": f"{lfc:.8g}",
                "pvalue": row.get(pvalue_col, "") if pvalue_col else "",
                "padj": f"{padj:.8g}",
                "deg_status": status,
            }
        )
        out.append(item)
    return out


def load_epigenetic_catalog() -> list[dict[str, str]]:
    out = []
    seen = set()
    for source_label, path in [("primary", env("RNA_GENE_CATALOG")), ("supplemental", env("RNA_GENE_CATALOG_EXTRA"))]:
        header, rows = read_table(path)
        if not header:
            continue
        gene_col = first_col(header, ["matched_gene_id", "gene_id", "query", "id"])
        group_col = first_col(header, ["group", "category", "class"])
        name_col = first_col(header, ["gene_name", "symbol", "name"])
        query_col = first_col(header, ["query", "query_display"])
        for row in rows:
            gid = row.get(gene_col, "") if gene_col else ""
            if not gid:
                continue
            group = row.get(group_col, "epigenetic_machinery") if group_col else "epigenetic_machinery"
            key = (gid, group)
            if key in seen:
                continue
            seen.add(key)
            item = dict(row)
            item["gene_id"] = gid
            item["gene_name"] = row.get(name_col, gid) if name_col else gid
            item["machinery_group"] = group
            item["query_id"] = row.get(query_col, gid) if query_col else gid
            item["catalog_source"] = source_label
            item["catalog_path"] = path
            out.append(item)
    return out


def command_validate(_args: argparse.Namespace) -> None:
    rows = []

    def check_file(label: str, path: str, required: bool = False) -> None:
        ok = path_exists(path)
        status = "ok" if ok else ("error" if required else "warning")
        rows.append({"item": label, "path": path, "status": status, "message": "" if ok else "missing or empty"})

    check_file("RNA counts matrix", env("RNA_COUNTS_MATRIX"), True)
    check_file("RNA normalized matrix", env("RNA_NORMALIZED_MATRIX"), True)
    check_file("RNA DEG results", env("RNA_DEG_RESULTS"), True)
    check_file("RNA metadata", env("RNA_METADATA_FILE"), False)
    check_file("RNA epigenetic machinery catalog", env("RNA_GENE_CATALOG"), False)
    check_file("RNA supplemental epigenetic catalog", env("RNA_GENE_CATALOG_EXTRA"), False)
    check_file("RNA expression by context", env("RNA_EXPRESSION_CONTEXT"), False)
    check_file("RNA WGCNA hits", env("RNA_WGCNA_HITS"), False)
    check_file("RNA Mfuzz hits", env("RNA_MFUZZ_HITS"), False)
    check_file("RNA DTU hits", env("RNA_DTU_HITS"), False)
    check_file("RNA splicing hits", env("RNA_SPLICING_HITS"), False)
    check_file("ChIP metadata", env("CHIP_METADATA_FILE"), False)
    check_file("ChIP differential binding", env("CHIP_DIFF_BINDING_FILE"), False)
    check_file("Genome annotation", env("ANNOTATION_FILE"), False)
    check_file("Functional annotation", env("FUNCTIONAL_ANNOTATION"), False)

    for label, pattern in [
        ("ChIP annotated peaks", env("CHIP_ANNOTATED_PEAKS_GLOB")),
        ("ChIP peak BEDs", env("CHIP_PEAK_BED_GLOB")),
        ("ChIP peak counts", env("CHIP_PEAK_COUNT_GLOB")),
    ]:
        matches = glob_existing(pattern)
        rows.append(
            {
                "item": label,
                "path": pattern,
                "status": "ok" if matches else "warning",
                "message": f"{len(matches)} matching file(s)",
            }
        )

    counts_genes = set(read_matrix_gene_ids(env("RNA_COUNTS_MATRIX")))
    norm_genes = set(read_matrix_gene_ids(env("RNA_NORMALIZED_MATRIX")))
    deg_genes = {r["gene_id"] for r in normalize_deg_rows()}
    if counts_genes and norm_genes:
        overlap = len(counts_genes & norm_genes)
        rows.append({"item": "RNA matrix gene overlap", "path": "", "status": "ok" if overlap else "error", "message": f"{overlap} genes"})
    if deg_genes and norm_genes:
        overlap = len(deg_genes & norm_genes)
        rows.append({"item": "RNA DEG/expression gene overlap", "path": "", "status": "ok" if overlap else "warning", "message": f"{overlap} genes"})

    write_table(outdir("010-input-validation") / "validation_report.tsv", rows, ["item", "path", "status", "message"])
    errors = [r for r in rows if r["status"] == "error"]
    summary = ["# Input validation", "", f"Generated: {dt.datetime.now().isoformat(timespec='seconds')}", ""]
    summary.extend(f"- {r['status'].upper()}: {r['item']} {r['message']}".strip() for r in rows)
    write_text(outdir("010-input-validation") / "validation_report.md", "\n".join(summary) + "\n")
    if errors:
        print("ERROR: required inputs failed validation. See 010-input-validation/validation_report.tsv", file=sys.stderr)
        raise SystemExit(1)


def command_prepare(_args: argparse.Namespace) -> None:
    manifest = []
    for key in [
        "RNA_COUNTS_MATRIX",
        "RNA_NORMALIZED_MATRIX",
        "RNA_DEG_RESULTS",
        "RNA_METADATA_FILE",
        "CHIP_METADATA_FILE",
        "CHIP_DIFF_BINDING_FILE",
        "ANNOTATION_FILE",
        "FUNCTIONAL_ANNOTATION",
        "GENES_OF_INTEREST_FILE",
    ]:
        manifest.append({"input_key": key, "path": env(key), "exists": str(path_exists(env(key)))})
    for key in ["CHIP_ANNOTATED_PEAKS_GLOB", "CHIP_PEAK_BED_GLOB", "CHIP_PEAK_COUNT_GLOB"]:
        matches = glob_existing(env(key))
        manifest.append({"input_key": key, "path": env(key), "exists": str(bool(matches)), "matches": ",".join(matches)})
    write_table(outdir("020-prepared-inputs") / "input_manifest.tsv", manifest, ["input_key", "path", "exists", "matches"])

    deg = normalize_deg_rows()
    deg_header = ["contrast_id", "gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "deg_status"]
    extra = sorted({k for r in deg for k in r.keys()} - set(deg_header))
    write_table(outdir("020-prepared-inputs") / "rnaseq_deg_normalized.tsv", deg, deg_header + extra)

    meta_rows = []
    for source, path in [("rna", env("RNA_METADATA_FILE")), ("chip", env("CHIP_METADATA_FILE"))]:
        header, rows = read_table(path)
        for row in rows:
            item = {"assay": source}
            item.update(row)
            meta_rows.append(item)
    header = sorted({k for r in meta_rows for k in r.keys()})
    if "assay" in header:
        header.remove("assay")
        header = ["assay"] + header
    write_table(outdir("020-prepared-inputs") / "metadata_combined.tsv", meta_rows, header or ["assay"])


def command_harmonize(_args: argparse.Namespace) -> None:
    genes = parse_annotation_genes(env("ANNOTATION_FILE"))
    gene_interest = read_set_file(env("GENES_OF_INTEREST_FILE"))
    epigenetic_catalog = load_epigenetic_catalog()
    epigenetic_gene_ids = {r["gene_id"] for r in epigenetic_catalog}
    gene_interest.update(epigenetic_gene_ids)
    deg = normalize_deg_rows()
    deg_by_gene = {r["gene_id"]: r for r in deg}
    functional: dict[str, dict[str, str]] = {}
    f_header, f_rows = read_table(env("FUNCTIONAL_ANNOTATION"))
    f_gene_col = first_col(f_header, ["gene_id", "gene", "id"]) if f_header else ""
    for row in f_rows:
        if f_gene_col and row.get(f_gene_col):
            functional[row[f_gene_col]] = row
    all_gene_ids = set(genes)
    all_gene_ids.update(read_matrix_gene_ids(env("RNA_COUNTS_MATRIX")))
    all_gene_ids.update(read_matrix_gene_ids(env("RNA_NORMALIZED_MATRIX")))
    all_gene_ids.update(deg_by_gene)
    all_gene_ids.update(gene_interest)
    rows = []
    lost = []
    for gid in sorted(all_gene_ids):
        base = genes.get(gid, {"gene_id": gid, "gene_name": gid, "chromosome": "", "start": "", "end": "", "strand": "", "biotype": ""})
        deg_row = deg_by_gene.get(gid, {})
        fn_row = functional.get(gid, {})
        gene_name = base.get("gene_name") or deg_row.get("gene_name") or gid
        biotype = base.get("biotype") or deg_row.get("biotype") or fn_row.get("biotype", "")
        fn_text = fn_row.get("functional_annotation") or fn_row.get("description") or fn_row.get("term") or ""
        row = {
            "gene_id": gid,
            "gene_name": gene_name,
            "chromosome": base.get("chromosome", ""),
            "start": base.get("start", ""),
            "end": base.get("end", ""),
            "strand": base.get("strand", ""),
            "biotype": biotype,
            "functional_annotation": fn_text,
            "is_gene_of_interest": str(gid in gene_interest).lower(),
            "is_epigenetic_machinery": str(gid in epigenetic_gene_ids).lower(),
        }
        if not row["chromosome"]:
            lost.append({"gene_id": gid, "reason": "missing_genomic_annotation"})
        rows.append(row)
    header = [
        "gene_id",
        "gene_name",
        "chromosome",
        "start",
        "end",
        "strand",
        "biotype",
        "functional_annotation",
        "is_gene_of_interest",
        "is_epigenetic_machinery",
    ]
    write_table(outdir("030-id-harmonization") / "gene_master_table.tsv", rows, header)
    write_table(outdir("030-id-harmonization") / "unmapped_genes.tsv", lost, ["gene_id", "reason"])
    catalog_header = ["gene_id", "gene_name", "machinery_group", "query_id"] + [
        h for h in sorted({k for r in epigenetic_catalog for k in r}) if h not in {"gene_id", "gene_name", "machinery_group", "query_id"}
    ]
    write_table(outdir("030-id-harmonization") / "epigenetic_machinery_catalog.tsv", epigenetic_catalog, catalog_header or ["gene_id"])


def load_gene_master() -> dict[str, dict[str, str]]:
    _header, rows = read_table(str(outdir("030-id-harmonization") / "gene_master_table.tsv"))
    return {r["gene_id"]: r for r in rows if r.get("gene_id")}


def annotated_peak_rows() -> list[dict[str, str]]:
    rows = []
    for path in glob_existing(env("CHIP_ANNOTATED_PEAKS_GLOB")):
        header, data = read_table(path)
        if not header:
            continue
        for idx, row in enumerate(data, start=1):
            gene_col = first_col(header, ["nearest_gene_id", "gene_id", "associated_gene", "target_gene"])
            name_col = first_col(header, ["nearest_gene_name", "gene_name", "gene"])
            class_col = first_col(header, ["genomic_annotation", "annotation", "class", "peak_class"])
            mark = row.get(first_col(header, [env("MARK_COLUMN", "mark_or_factor"), "mark", "factor"]), "")
            if not mark:
                mark = Path(path).stem.replace(".annotated", "")
            peak_id = row.get(first_col(header, ["peak_id", "id", "name"]), "") or f"{Path(path).stem}_{idx}"
            rows.append(
                {
                    "peak_id": peak_id,
                    "chrom": row.get(first_col(header, ["chrom", "chr", "seqnames"]), ""),
                    "start": row.get(first_col(header, ["start", "chromStart"]), ""),
                    "end": row.get(first_col(header, ["end", "chromEnd"]), ""),
                    "mark_or_factor": mark,
                    "condition": row.get(first_col(header, [env("CONDITION_COLUMN", "condition"), "stage", "treatment"]), ""),
                    "associated_gene_id": row.get(gene_col, "") if gene_col else "",
                    "associated_gene_name": row.get(name_col, "") if name_col else "",
                    "distance_to_tss": row.get(first_col(header, ["distance_to_tss", "distanceToTSS", "distance"]), ""),
                    "genomic_annotation": row.get(class_col, "") if class_col else "",
                    "link_mode": "annotated_peak_table",
                    "source_file": norm_path(path),
                }
            )
    return rows


def bed_peak_rows_from_master(gene_master: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    gene_by_chrom: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in gene_master.values():
        if row.get("chromosome") and row.get("start") and row.get("end"):
            gene_by_chrom[row["chromosome"]].append(row)
    rows = []
    window = as_int(env("PEAK_GENE_WINDOW_BP", "5000"), 5000)
    for path in glob_existing(env("CHIP_PEAK_BED_GLOB")):
        with open_text(path) as handle:
            for idx, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom, start, end = parts[0], parts[1], parts[2]
                midpoint = (as_int(start) + as_int(end)) // 2
                best = None
                best_dist = None
                for gene in gene_by_chrom.get(chrom, []):
                    gs, ge = as_int(gene.get("start")), as_int(gene.get("end"))
                    dist = 0 if gs <= midpoint <= ge else min(abs(midpoint - gs), abs(midpoint - ge))
                    if dist <= window and (best_dist is None or dist < best_dist):
                        best, best_dist = gene, dist
                if best:
                    source_name = Path(path).stem.replace(".consensus", "")
                    rows.append(
                        {
                            "peak_id": parts[3] if len(parts) > 3 and parts[3] else f"{source_name}_{idx}",
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "mark_or_factor": source_name,
                            "condition": "",
                            "associated_gene_id": best["gene_id"],
                            "associated_gene_name": best.get("gene_name", ""),
                            "distance_to_tss": str(best_dist or 0),
                            "genomic_annotation": "window",
                            "link_mode": "window_based",
                            "source_file": norm_path(path),
                        }
                    )
    return rows


def command_map_peaks(_args: argparse.Namespace) -> None:
    gene_master = load_gene_master()
    rows = annotated_peak_rows()
    if not rows:
        rows = bed_peak_rows_from_master(gene_master)
    for row in rows:
        ann = row.get("genomic_annotation", "").lower()
        row["promoter_flag"] = str("promoter" in ann).lower()
    header = [
        "peak_id",
        "chrom",
        "start",
        "end",
        "mark_or_factor",
        "condition",
        "associated_gene_id",
        "associated_gene_name",
        "distance_to_tss",
        "genomic_annotation",
        "promoter_flag",
        "link_mode",
        "source_file",
    ]
    write_table(outdir("040-peak-gene-mapping") / "peak_to_gene.tsv", rows, header)
    promoter = [r for r in rows if r.get("promoter_flag") == "true"]
    distal = [r for r in rows if r.get("promoter_flag") != "true"]
    write_table(outdir("040-peak-gene-mapping") / "promoter_peak_gene_links.tsv", promoter, header)
    write_table(outdir("040-peak-gene-mapping") / "distal_peak_gene_links.tsv", distal, header)
    summary = []
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("associated_gene_id"):
            by_gene[row["associated_gene_id"]].append(row)
    for gid, links in sorted(by_gene.items()):
        anns = [x.get("genomic_annotation", "").lower() for x in links]
        summary.append(
            {
                "gene_id": gid,
                "total_associated_peaks": len(links),
                "promoter_peaks": sum("promoter" in a for a in anns),
                "gene_body_peaks": sum(a in {"gene", "exon", "intron", "gene_body"} for a in anns),
                "distal_peaks": sum("promoter" not in a and a not in {"gene", "exon", "intron", "gene_body"} for a in anns),
                "marks_or_factors": ";".join(sorted({x.get("mark_or_factor", "") for x in links if x.get("mark_or_factor")})),
                "conditions": ";".join(sorted({x.get("condition", "") for x in links if x.get("condition")})),
            }
        )
    write_table(
        outdir("040-peak-gene-mapping") / "gene_to_peak_summary.tsv",
        summary,
        ["gene_id", "total_associated_peaks", "promoter_peaks", "gene_body_peaks", "distal_peaks", "marks_or_factors", "conditions"],
    )


def command_summarize_rna(_args: argparse.Namespace) -> None:
    header, rows = read_table(env("RNA_NORMALIZED_MATRIX"))
    groups = metadata_groups(env("RNA_METADATA_FILE"))
    gene_col = first_col(header, [env("GENE_ID_COLUMN", "gene_id"), "gene", "id"]) if header else ""
    if not gene_col and header:
        gene_col = header[0]
    sample_cols = [h for h in header if h != gene_col]
    deg = normalize_deg_rows()
    by_gene_deg: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in deg:
        by_gene_deg[row["gene_id"]].append(row)
    summaries = []
    for row in rows:
        gid = row.get(gene_col, "")
        if not gid:
            continue
        values_by_group: dict[str, list[float]] = defaultdict(list)
        all_values = []
        for sample in sample_cols:
            value = as_float(row.get(sample, "0"))
            values_by_group[sample_group(sample, groups)].append(value)
            all_values.append(value)
        group_means = {g: (sum(v) / len(v) if v else 0.0) for g, v in values_by_group.items()}
        top_group = max(group_means, key=group_means.get) if group_means else ""
        hits = by_gene_deg.get(gid, [])
        statuses = [h.get("deg_status", "not_significant") for h in hits]
        padjs = [as_float(h.get("padj"), 1.0) for h in hits if h.get("padj")]
        lfcs = [abs(as_float(h.get("log2FoldChange"))) for h in hits]
        max_expr = max(all_values) if all_values else 0.0
        total_expr = sum(all_values) if all_values else 0.0
        summaries.append(
            {
                "gene_id": gid,
                "mean_expression": f"{statistics.mean(all_values):.8g}" if all_values else "0",
                "median_expression": f"{statistics.median(all_values):.8g}" if all_values else "0",
                "max_expression": f"{max_expr:.8g}",
                "min_expression": f"{min(all_values):.8g}" if all_values else "0",
                "context_with_highest_expression": top_group,
                "expression_specificity": f"{(max_expr / total_expr):.8g}" if total_expr else "0",
                "n_deg_records": len(hits),
                "n_significant_contrasts": sum(s in {"up", "down"} for s in statuses),
                "n_up": statuses.count("up"),
                "n_down": statuses.count("down"),
                "max_abs_log2FC": f"{max(lfcs):.8g}" if lfcs else "0",
                "min_padj": f"{min(padjs):.8g}" if padjs else "1",
                "transcriptional_dynamism_score": f"{(math.log2(max_expr + 1) * (max(lfcs) if lfcs else 0)):.8g}",
            }
        )
    sum_header = [
        "gene_id",
        "mean_expression",
        "median_expression",
        "max_expression",
        "min_expression",
        "context_with_highest_expression",
        "expression_specificity",
        "n_deg_records",
        "n_significant_contrasts",
        "n_up",
        "n_down",
        "max_abs_log2FC",
        "min_padj",
        "transcriptional_dynamism_score",
    ]
    write_table(outdir("050-rnaseq-summary") / "rna_gene_summary.tsv", summaries, sum_header)
    deg_header = ["contrast_id", "gene_id", "gene_name", "log2FoldChange", "pvalue", "padj", "deg_status"]
    write_table(outdir("050-rnaseq-summary") / "rna_deg_long.tsv", deg, deg_header + sorted({k for r in deg for k in r} - set(deg_header)))


def load_diff_binding() -> list[dict[str, str]]:
    header, rows = read_table(env("CHIP_DIFF_BINDING_FILE"))
    if not header:
        return []
    gene_col = first_col(header, ["gene_id", "nearest_gene_id", "associated_gene_id"])
    peak_col = first_col(header, ["peak_id", "id", "region", "feature_id"])
    mark_col = first_col(header, [env("MARK_COLUMN", "mark_or_factor"), "mark", "factor", "peak_set"])
    lfc_col = first_col(header, ["log2FoldChange", "log2FC", "logFC"])
    padj_col = first_col(header, ["padj", "FDR", "qvalue", "adj.P.Val"])
    out = []
    for row in rows:
        item = dict(row)
        item["gene_id"] = row.get(gene_col, "") if gene_col else ""
        item["peak_id"] = row.get(peak_col, "") if peak_col else ""
        item["mark_or_factor"] = row.get(mark_col, "") if mark_col else ""
        item["chip_log2FC"] = row.get(lfc_col, "0") if lfc_col else "0"
        item["chip_padj"] = row.get(padj_col, "1") if padj_col else "1"
        lfc = as_float(item["chip_log2FC"])
        padj = as_float(item["chip_padj"], 1.0)
        status = "not_significant"
        if padj <= as_float(env("DIFF_BINDING_PADJ_THRESHOLD", "0.05"), 0.05):
            if lfc >= as_float(env("DIFF_BINDING_LOG2FC_THRESHOLD", "1"), 1.0):
                status = "gained"
            elif lfc <= -as_float(env("DIFF_BINDING_LOG2FC_THRESHOLD", "1"), 1.0):
                status = "lost"
        item["binding_status"] = status
        out.append(item)
    return out


def command_summarize_chip(_args: argparse.Namespace) -> None:
    _header, links = read_table(str(outdir("040-peak-gene-mapping") / "peak_to_gene.tsv"))
    diff = load_diff_binding()
    diff_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in diff:
        if row.get("gene_id"):
            diff_by_gene[row["gene_id"]].append(row)
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in links:
        if row.get("associated_gene_id"):
            by_gene[row["associated_gene_id"]].append(row)
    rows = []
    for gid in sorted(set(by_gene) | set(diff_by_gene)):
        links_g = by_gene.get(gid, [])
        diff_g = diff_by_gene.get(gid, [])
        anns = [x.get("genomic_annotation", "").lower() for x in links_g]
        statuses = [x.get("binding_status", "not_tested") for x in diff_g]
        rows.append(
            {
                "gene_id": gid,
                "total_associated_peaks": len(links_g),
                "promoter_peaks": sum("promoter" in a for a in anns),
                "gene_body_peaks": sum(a in {"gene", "exon", "intron", "gene_body"} for a in anns),
                "distal_peaks": sum("promoter" not in a and a not in {"gene", "exon", "intron", "gene_body"} for a in anns),
                "marks_or_factors": ";".join(sorted({x.get("mark_or_factor", "") for x in links_g + diff_g if x.get("mark_or_factor")})),
                "conditions_with_peak": ";".join(sorted({x.get("condition", "") for x in links_g if x.get("condition")})),
                "n_differential_peaks": sum(s in {"gained", "lost"} for s in statuses),
                "max_abs_chip_log2FC": f"{max([abs(as_float(x.get('chip_log2FC'))) for x in diff_g] or [0]):.8g}",
                "min_chip_padj": f"{min([as_float(x.get('chip_padj'), 1.0) for x in diff_g] or [1]):.8g}",
                "binding_statuses": ";".join(sorted(set(statuses))) if statuses else "not_tested",
            }
        )
    header = [
        "gene_id",
        "total_associated_peaks",
        "promoter_peaks",
        "gene_body_peaks",
        "distal_peaks",
        "marks_or_factors",
        "conditions_with_peak",
        "n_differential_peaks",
        "max_abs_chip_log2FC",
        "min_chip_padj",
        "binding_statuses",
    ]
    write_table(outdir("060-chipseq-summary") / "chip_gene_summary.tsv", rows, header)
    write_table(outdir("060-chipseq-summary") / "chip_differential_long.tsv", diff, sorted({k for r in diff for k in r}) or ["gene_id"])
    mark_stage = chip_metadata_mark_stage_rows()
    write_table(
        outdir("060-chipseq-summary") / "chip_mark_stage_metadata.tsv",
        mark_stage,
        ["mark_or_factor", "stage_or_condition", "n_samples", "replicates", "batches", "sample_ids"],
    )


def load_mark_config() -> dict[str, dict[str, str]]:
    path = Path(env("PROJECT_DIR", ".")) / "config" / "chip_marks_config.tsv"
    header, rows = read_table(str(path))
    return {r.get("mark_or_factor", "unknown"): r for r in rows if r.get("mark_or_factor")}


def chip_metadata_mark_stage_rows() -> list[dict[str, object]]:
    header, rows = read_table(env("CHIP_METADATA_FILE"))
    if not header:
        return []
    mark_col = first_col(header, [env("MARK_COLUMN", "mark_or_factor"), "mark", "factor"])
    condition_col = first_col(header, [env("CONDITION_COLUMN", "condition"), "stage", "life_stage"])
    sample_col = first_col(header, [env("SAMPLE_ID_COLUMN", "sample_id"), "sample"])
    rep_col = first_col(header, ["replicate", "rep", "biological_replicate"])
    batch_col = first_col(header, ["batch"])
    counts: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        mark = row.get(mark_col, "") if mark_col else "unknown"
        stage = row.get(condition_col, "") if condition_col else "unknown"
        key = (mark or "unknown", stage or "unknown")
        item = counts.setdefault(
            key,
            {"mark_or_factor": key[0], "stage_or_condition": key[1], "n_samples": 0, "replicates": set(), "batches": set(), "sample_ids": []},
        )
        item["n_samples"] = int(item["n_samples"]) + 1
        if rep_col and row.get(rep_col):
            item["replicates"].add(row[rep_col])
        if batch_col and row.get(batch_col):
            item["batches"].add(row[batch_col])
        if sample_col and row.get(sample_col):
            item["sample_ids"].append(row[sample_col])
    out = []
    for item in counts.values():
        out.append(
            {
                "mark_or_factor": item["mark_or_factor"],
                "stage_or_condition": item["stage_or_condition"],
                "n_samples": item["n_samples"],
                "replicates": ";".join(sorted(item["replicates"])),
                "batches": ";".join(sorted(item["batches"])),
                "sample_ids": ";".join(item["sample_ids"]),
            }
        )
    return sorted(out, key=lambda r: (str(r["stage_or_condition"]), str(r["mark_or_factor"])))


def expression_context_lookup() -> dict[tuple[str, str], dict[str, str]]:
    header, rows = read_table(env("RNA_EXPRESSION_CONTEXT"))
    if not header:
        return {}
    gene_col = first_col(header, ["gene_id", "matched_gene_id", "gene"])
    stage_col = first_col(header, ["stage", "condition", "stage_class", "context"])
    lookup = {}
    for row in rows:
        gid = row.get(gene_col, "") if gene_col else ""
        stage = row.get(stage_col, "") if stage_col else ""
        if gid and stage:
            lookup[(gid, stage.lower())] = row
    return lookup


def supplemental_rna_evidence() -> dict[str, dict[str, str]]:
    specs = [
        ("wgcna", env("RNA_WGCNA_HITS"), ["module_color", "module_label", "is_hub", "kME", "abs_kME"]),
        ("mfuzz", env("RNA_MFUZZ_HITS"), ["cluster", "membership"]),
        ("dtu", env("RNA_DTU_HITS"), ["transcript_id", "variable", "delta_usage", "padj"]),
        ("splicing", env("RNA_SPLICING_HITS"), ["event_type", "splicing_contrast", "FDR", "IncLevelDifference"]),
    ]
    evidence: dict[str, dict[str, str]] = defaultdict(dict)
    for label, path, cols in specs:
        header, rows = read_table(path)
        gene_col = first_col(header, ["gene_id", "matched_gene_id", "gene"]) if header else ""
        if not gene_col:
            continue
        for row in rows:
            gid = row.get(gene_col, "")
            if not gid:
                continue
            evidence[gid][f"{label}_hit"] = "true"
            values = []
            for col in cols:
                if col in row and row[col]:
                    values.append(f"{col}={row[col]}")
            if values:
                prev = evidence[gid].get(f"{label}_summary", "")
                joined = ";".join(values)
                evidence[gid][f"{label}_summary"] = joined if not prev else f"{prev}|{joined}"
    return evidence


def integration_class(deg_status: str, chip: dict[str, str]) -> str:
    total = as_int(chip.get("total_associated_peaks", "0"))
    promoter = as_int(chip.get("promoter_peaks", "0"))
    body = as_int(chip.get("gene_body_peaks", "0"))
    distal = as_int(chip.get("distal_peaks", "0"))
    diff = as_int(chip.get("n_differential_peaks", "0"))
    is_deg = deg_status in {"up", "down"}
    if is_deg and diff:
        return "DEG_with_differential_peak"
    if is_deg and promoter:
        return "DEG_with_promoter_peak"
    if is_deg and body:
        return "DEG_with_gene_body_peak"
    if is_deg and distal:
        return "DEG_with_distal_peak"
    if is_deg:
        return "DEG_only"
    if total:
        return "ChIP_only"
    return "unchanged"


def write_gene_mark_stage_tables(
    gene_master: dict[str, dict[str, str]],
    rna: dict[str, dict[str, str]],
    deg_by_gene: dict[str, list[dict[str, str]]],
    peak_links: list[dict[str, str]],
    chip: dict[str, dict[str, str]],
) -> None:
    mark_config = load_mark_config()
    expr_context = expression_context_lookup()
    relations = []
    for link in peak_links:
        gid = link.get("associated_gene_id", "")
        if not gid:
            continue
        mark = link.get("mark_or_factor", "") or "unknown"
        stage = link.get("condition", "") or "unknown"
        gene = gene_master.get(gid, {"gene_id": gid, "gene_name": gid})
        expr = expr_context.get((gid, stage.lower()), {})
        deg_hits = deg_by_gene.get(gid, [])
        status = "not_significant"
        if any(d.get("deg_status") == "up" for d in deg_hits):
            status = "up"
        elif any(d.get("deg_status") == "down" for d in deg_hits):
            status = "down"
        mark_row = mark_config.get(mark, mark_config.get("unknown", {}))
        relations.append(
            {
                "gene_id": gid,
                "gene_name": gene.get("gene_name", gid),
                "is_epigenetic_machinery": gene.get("is_epigenetic_machinery", "false"),
                "mark_or_factor": mark,
                "stage_or_condition": stage,
                "regulatory_class": mark_row.get("regulatory_class", "unknown"),
                "expected_effect": mark_row.get("expected_effect", "unknown"),
                "peak_id": link.get("peak_id", ""),
                "peak_location": link.get("genomic_annotation", ""),
                "promoter_flag": link.get("promoter_flag", "false"),
                "distance_to_tss": link.get("distance_to_tss", ""),
                "rna_mean_TPM_in_stage": expr.get("mean_TPM", ""),
                "rna_mean_log2TPM_in_stage": expr.get("mean_log2TPM", ""),
                "rna_fraction_expressed_in_stage": expr.get("fraction_expressed", ""),
                "representative_deg_status": status,
                "max_abs_log2FC": max([as_float(d.get("log2FoldChange")) for d in deg_hits] or [0]),
                "min_padj": min([as_float(d.get("padj"), 1.0) for d in deg_hits] or [1]),
                "source_file": link.get("source_file", ""),
            }
        )
    rel_header = [
        "gene_id",
        "gene_name",
        "is_epigenetic_machinery",
        "mark_or_factor",
        "stage_or_condition",
        "regulatory_class",
        "expected_effect",
        "peak_id",
        "peak_location",
        "promoter_flag",
        "distance_to_tss",
        "rna_mean_TPM_in_stage",
        "rna_mean_log2TPM_in_stage",
        "rna_fraction_expressed_in_stage",
        "representative_deg_status",
        "max_abs_log2FC",
        "min_padj",
        "source_file",
    ]
    write_table(outdir("070-integrated-tables") / "gene_mark_stage_links.tsv", relations, rel_header)

    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in relations:
        grouped[(str(row["gene_id"]), str(row["mark_or_factor"]), str(row["stage_or_condition"]))].append(row)
    summary = []
    for (gid, mark, stage), items in sorted(grouped.items()):
        gene = gene_master.get(gid, {"gene_id": gid, "gene_name": gid})
        chip_row = chip.get(gid, {})
        summary.append(
            {
                "gene_id": gid,
                "gene_name": gene.get("gene_name", gid),
                "is_epigenetic_machinery": gene.get("is_epigenetic_machinery", "false"),
                "mark_or_factor": mark,
                "stage_or_condition": stage,
                "n_peaks": len(items),
                "n_promoter_peaks": sum(str(x.get("promoter_flag", "")).lower() == "true" for x in items),
                "peak_locations": ";".join(sorted({str(x.get("peak_location", "")) for x in items if x.get("peak_location")})),
                "regulatory_class": items[0].get("regulatory_class", "unknown"),
                "expected_effect": items[0].get("expected_effect", "unknown"),
                "rna_mean_TPM_in_stage": items[0].get("rna_mean_TPM_in_stage", ""),
                "representative_deg_status": items[0].get("representative_deg_status", ""),
                "gene_total_associated_peaks": chip_row.get("total_associated_peaks", "0"),
                "gene_n_differential_peaks": chip_row.get("n_differential_peaks", "0"),
            }
        )
    sum_header = [
        "gene_id",
        "gene_name",
        "is_epigenetic_machinery",
        "mark_or_factor",
        "stage_or_condition",
        "n_peaks",
        "n_promoter_peaks",
        "peak_locations",
        "regulatory_class",
        "expected_effect",
        "rna_mean_TPM_in_stage",
        "representative_deg_status",
        "gene_total_associated_peaks",
        "gene_n_differential_peaks",
    ]
    write_table(outdir("070-integrated-tables") / "gene_mark_stage_summary.tsv", summary, sum_header)

    mark_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in summary:
        mark_groups[str(row["mark_or_factor"])].append(row)
    by_mark = [
        {
            "mark_or_factor": mark,
            "n_linked_genes": len({str(x["gene_id"]) for x in items}),
            "n_epigenetic_machinery_genes": len({str(x["gene_id"]) for x in items if str(x.get("is_epigenetic_machinery", "")).lower() == "true"}),
            "stages_or_conditions": ";".join(sorted({str(x["stage_or_condition"]) for x in items})),
            "n_promoter_gene_links": sum(as_int(x.get("n_promoter_peaks", "0")) > 0 for x in items),
            "n_deg_linked_genes": len({str(x["gene_id"]) for x in items if x.get("representative_deg_status") in {"up", "down"}}),
        }
        for mark, items in sorted(mark_groups.items())
    ]
    write_table(
        outdir("070-integrated-tables") / "mark_to_gene_catalog.tsv",
        by_mark,
        ["mark_or_factor", "n_linked_genes", "n_epigenetic_machinery_genes", "stages_or_conditions", "n_promoter_gene_links", "n_deg_linked_genes"],
    )


def command_integrate(_args: argparse.Namespace) -> None:
    gene_master = load_gene_master()
    _h, rna_rows = read_table(str(outdir("050-rnaseq-summary") / "rna_gene_summary.tsv"))
    _h, deg_rows = read_table(str(outdir("050-rnaseq-summary") / "rna_deg_long.tsv"))
    _h, chip_rows = read_table(str(outdir("060-chipseq-summary") / "chip_gene_summary.tsv"))
    _h, peak_links = read_table(str(outdir("040-peak-gene-mapping") / "peak_to_gene.tsv"))
    rna = {r["gene_id"]: r for r in rna_rows if r.get("gene_id")}
    chip = {r["gene_id"]: r for r in chip_rows if r.get("gene_id")}
    rna_evidence = supplemental_rna_evidence()
    deg_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in deg_rows:
        if row.get("gene_id"):
            deg_by_gene[row["gene_id"]].append(row)
    all_genes = sorted(set(gene_master) | set(rna) | set(chip) | set(deg_by_gene))
    gene_rows = []
    contrast_rows = []
    for gid in all_genes:
        base = dict(gene_master.get(gid, {"gene_id": gid, "gene_name": gid}))
        rna_row = rna.get(gid, {})
        chip_row = chip.get(gid, {})
        deg_hits = deg_by_gene.get(gid, [])
        best_status = "not_significant"
        if any(d.get("deg_status") == "up" for d in deg_hits):
            best_status = "up"
        elif any(d.get("deg_status") == "down" for d in deg_hits):
            best_status = "down"
        base.update({f"rna_{k}": v for k, v in rna_row.items() if k != "gene_id"})
        base.update({f"chip_{k}": v for k, v in chip_row.items() if k != "gene_id"})
        base.update(rna_evidence.get(gid, {}))
        base["representative_deg_status"] = best_status
        base["integrative_class"] = integration_class(best_status, chip_row)
        gene_rows.append(base)
        if deg_hits:
            for d in deg_hits:
                row = dict(base)
                row.update(
                    {
                        "contrast_id": d.get("contrast_id", "default"),
                        "log2FoldChange": d.get("log2FoldChange", "0"),
                        "padj": d.get("padj", "1"),
                        "deg_status": d.get("deg_status", "not_significant"),
                        "integrative_class": integration_class(d.get("deg_status", "not_significant"), chip_row),
                    }
                )
                contrast_rows.append(row)
        else:
            row = dict(base)
            row.update({"contrast_id": "not_tested", "log2FoldChange": "0", "padj": "1", "deg_status": "not_significant"})
            contrast_rows.append(row)
    gene_header = sorted({k for r in gene_rows for k in r})
    for key in ["gene_id", "gene_name", "integrative_class"]:
        if key in gene_header:
            gene_header.remove(key)
    gene_header = ["gene_id", "gene_name", "integrative_class"] + gene_header
    write_table(outdir("070-integrated-tables") / "integrated_gene_table.tsv", gene_rows, gene_header)
    contrast_header = ["contrast_id"] + [h for h in gene_header if h != "contrast_id"]
    for extra in ["log2FoldChange", "padj", "deg_status"]:
        if extra not in contrast_header:
            contrast_header.append(extra)
    write_table(outdir("070-integrated-tables") / "integrated_by_contrast.tsv", contrast_rows, contrast_header)
    counts = [{"integrative_class": k, "n_genes": v} for k, v in sorted(Counter(r["integrative_class"] for r in gene_rows).items())]
    write_table(outdir("070-integrated-tables") / "integrative_class_counts.tsv", counts, ["integrative_class", "n_genes"])
    write_gene_mark_stage_tables(gene_master, rna, deg_by_gene, peak_links, chip)


def command_score(_args: argparse.Namespace) -> None:
    _h, rows = read_table(str(outdir("070-integrated-tables") / "integrated_gene_table.tsv"))
    scored = []
    for row in rows:
        score = 0.0
        padj = as_float(row.get("rna_min_padj", "1"), 1.0)
        score += min(10.0, -math.log10(max(padj, 1e-300))) if padj < 1 else 0.0
        score += min(5.0, as_float(row.get("rna_max_abs_log2FC", "0")))
        score += 2.0 if as_int(row.get("chip_promoter_peaks", "0")) > 0 else 0.0
        score += 2.0 if as_int(row.get("chip_n_differential_peaks", "0")) > 0 else 0.0
        score += 1.0 if row.get("is_gene_of_interest", "false").lower() == "true" else 0.0
        score += 2.0 if row.get("is_epigenetic_machinery", "false").lower() == "true" else 0.0
        score += min(3.0, as_int(row.get("rna_n_significant_contrasts", "0")) * 0.5)
        score += min(2.0, len([x for x in row.get("chip_marks_or_factors", "").split(";") if x]) * 0.5)
        item = dict(row)
        item["candidate_score"] = f"{score:.4f}"
        item["score_components"] = "deg_significance;rna_log2fc;promoter_peak;differential_peak;gene_interest;epigenetic_machinery;multi_contrast;multi_mark"
        scored.append(item)
    scored.sort(key=lambda r: as_float(r.get("candidate_score")), reverse=True)
    header = ["candidate_score", "gene_id", "gene_name", "integrative_class", "score_components"] + [
        h for h in sorted({k for r in scored for k in r}) if h not in {"candidate_score", "gene_id", "gene_name", "integrative_class", "score_components"}
    ]
    write_table(outdir("080-candidate-scoring") / "candidate_gene_scores.tsv", scored, header)
    top_n = as_int(env("TOP_CANDIDATES_N", "100"), 100)
    write_table(outdir("080-candidate-scoring") / "top_candidates.tsv", scored[:top_n], header)
    _h, contrast_rows = read_table(str(outdir("070-integrated-tables") / "integrated_by_contrast.tsv"))
    by_contrast = []
    score_by_gene = {r.get("gene_id"): r.get("candidate_score", "0") for r in scored}
    for row in contrast_rows:
        item = dict(row)
        item["candidate_score"] = score_by_gene.get(row.get("gene_id"), "0")
        by_contrast.append(item)
    by_contrast.sort(key=lambda r: (r.get("contrast_id", ""), -as_float(r.get("candidate_score"))))
    write_table(outdir("080-candidate-scoring") / "ranked_candidates_by_contrast.tsv", by_contrast, ["candidate_score"] + [h for h in sorted({k for r in by_contrast for k in r}) if h != "candidate_score"])
    by_mark = []
    for row in scored:
        marks = [m for m in row.get("chip_marks_or_factors", "").split(";") if m] or ["no_chip_mark"]
        for mark in marks:
            item = dict(row)
            item["mark_or_factor"] = mark
            by_mark.append(item)
    write_table(outdir("080-candidate-scoring") / "ranked_candidates_by_mark.tsv", by_mark, ["mark_or_factor"] + header)


def command_visualize(_args: argparse.Namespace) -> None:
    vis_dir = outdir("090-visualizations")
    _h, scores = read_table(str(outdir("080-candidate-scoring") / "candidate_gene_scores.tsv"))
    _h, classes = read_table(str(outdir("070-integrated-tables") / "integrative_class_counts.tsv"))
    _h, mark_stage = read_table(str(outdir("060-chipseq-summary") / "chip_mark_stage_metadata.tsv"))
    _h, gene_mark = read_table(str(outdir("070-integrated-tables") / "gene_mark_stage_summary.tsv"))
    _h, epi_catalog = read_table(str(outdir("030-id-harmonization") / "epigenetic_machinery_catalog.tsv"))
    manifest = []

    def svg_text(x: float, y: float, text: object, size: int = 12, anchor: str = "start", weight: str = "400") -> str:
        return (
            f"<text x='{x}' y='{y}' font-size='{size}' font-family='Arial, Helvetica, sans-serif' "
            f"text-anchor='{anchor}' font-weight='{weight}' fill='#111827'>{html.escape(str(text))}</text>"
        )

    def write_svg(name: str, width: int, height: int, body: str) -> None:
        svg = (
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>"
            "<rect width='100%' height='100%' fill='#ffffff'/>"
            f"{body}</svg>\n"
        )
        write_text(vis_dir / f"{name}.svg", svg)
        manifest.append({"figure": f"{name}.svg", "status": "created"})

    def fallback_barh(name: str, title: str, rows: list[tuple[str, float]], color: str = "#2563eb") -> None:
        rows = rows[:30] or [("no_data", 0)]
        width = 900
        height = max(260, 70 + len(rows) * 24)
        left = 220
        max_value = max([v for _label, v in rows] or [1]) or 1
        parts = [svg_text(width / 2, 32, title, 18, "middle", "700")]
        for i, (label, value) in enumerate(rows):
            y = 60 + i * 24
            bar_width = int((width - left - 80) * (value / max_value))
            parts.append(svg_text(12, y + 14, label[:34], 11))
            parts.append(f"<rect x='{left}' y='{y}' width='{bar_width}' height='16' fill='{color}' rx='2'/>")
            parts.append(svg_text(left + bar_width + 8, y + 13, f"{value:g}", 11))
        write_svg(name, width, height, "".join(parts))

    def fallback_matrix(name: str, title: str, rows: list[str], cols: list[str], values: list[list[int]], empty_message: str) -> None:
        width = max(720, 180 + len(cols) * 120)
        height = max(280, 95 + len(rows) * 38)
        parts = [svg_text(width / 2, 32, title, 18, "middle", "700")]
        if not rows or not cols:
            parts.append(svg_text(width / 2, height / 2, empty_message, 16, "middle", "700"))
            write_svg(name, width, height, "".join(parts))
            return
        max_value = max([max(row) for row in values] or [1]) or 1
        left = 170
        top = 72
        cell_w = 105
        cell_h = 30
        for j, col in enumerate(cols):
            parts.append(svg_text(left + j * cell_w + cell_w / 2, top - 12, col, 10, "middle", "700"))
        for i, row_label in enumerate(rows):
            y = top + i * cell_h
            parts.append(svg_text(12, y + 20, row_label, 10))
            for j, _col in enumerate(cols):
                value = values[i][j]
                alpha = 0.15 + 0.75 * (value / max_value)
                x = left + j * cell_w
                parts.append(f"<rect x='{x}' y='{y}' width='{cell_w - 4}' height='{cell_h - 4}' fill='#0891b2' fill-opacity='{alpha:.3f}'/>")
                parts.append(svg_text(x + cell_w / 2 - 2, y + 18, value, 11, "middle", "700"))
        write_svg(name, width, height, "".join(parts))

    def generate_fallback_svgs(reason: Exception) -> None:
        write_text(vis_dir / "visualization_warning.txt", f"Matplotlib unavailable; generated SVG fallback figures. Reason: {reason}\n")
        class_rows = [(r.get("integrative_class", "unknown"), float(as_int(r.get("n_genes", "0")))) for r in classes]
        fallback_barh("barplot_integrative_classes", "Integrative classes", class_rows, "#2563eb")
        top_rows = [(r.get("gene_id", ""), as_float(r.get("candidate_score"))) for r in scores[:30]]
        fallback_barh("top_candidate_scores", "Top candidate genes", top_rows, "#16a34a")
        group_counts = Counter(r.get("machinery_group", "unknown") for r in epi_catalog)
        fallback_barh("epigenetic_catalog_groups", "Epigenetic machinery catalog", [(k, float(v)) for k, v in group_counts.items()], "#0f766e")
        marks = sorted({r.get("mark_or_factor", "unknown") for r in mark_stage})
        stages = sorted({r.get("stage_or_condition", "unknown") for r in mark_stage})
        mark_idx = {m: i for i, m in enumerate(marks)}
        stage_idx = {s: i for i, s in enumerate(stages)}
        matrix = [[0 for _ in stages] for _ in marks]
        for row in mark_stage:
            if row.get("mark_or_factor") in mark_idx and row.get("stage_or_condition") in stage_idx:
                matrix[mark_idx[row["mark_or_factor"]]][stage_idx[row["stage_or_condition"]]] = as_int(row.get("n_samples", "0"))
        fallback_matrix("chip_mark_stage_matrix", "ChIP-seq marks by life-cycle stage", marks, stages, matrix, "No ChIP metadata available")
        relation_counts = Counter((r.get("mark_or_factor", "unknown"), r.get("stage_or_condition", "unknown")) for r in gene_mark)
        rel_marks = sorted({m for m, _s in relation_counts})
        rel_stages = sorted({s for _m, s in relation_counts})
        rel_matrix = [[relation_counts[(m, s)] for s in rel_stages] for m in rel_marks]
        fallback_matrix("gene_mark_stage_matrix", "Gene-mark-stage links", rel_marks, rel_stages, rel_matrix, "No gene-mark-stage links yet")
        workflow = (
            svg_text(450, 34, "Integrative analysis workflow", 18, "middle", "700")
            + "<rect x='45' y='80' width='180' height='70' rx='8' fill='#dbeafe' stroke='#334155'/>"
            + "<rect x='45' y='210' width='180' height='70' rx='8' fill='#dcfce7' stroke='#334155'/>"
            + "<rect x='360' y='145' width='190' height='80' rx='8' fill='#fef3c7' stroke='#334155'/>"
            + "<rect x='680' y='145' width='180' height='80' rx='8' fill='#ede9fe' stroke='#334155'/>"
            + svg_text(135, 110, "RNA-seq", 15, "middle", "700")
            + svg_text(135, 132, "expression + DEG", 11, "middle")
            + svg_text(135, 240, "ChIP-seq", 15, "middle", "700")
            + svg_text(135, 262, "marks + peaks", 11, "middle")
            + svg_text(455, 176, "Integration", 15, "middle", "700")
            + svg_text(455, 199, "gene x mark x stage", 11, "middle")
            + svg_text(770, 176, "Outputs", 15, "middle", "700")
            + svg_text(770, 199, "HTML + figures", 11, "middle")
            + "<path d='M225 115 L360 175' stroke='#334155' stroke-width='2' marker-end='url(#arrow)'/>"
            + "<path d='M225 245 L360 195' stroke='#334155' stroke-width='2' marker-end='url(#arrow)'/>"
            + "<path d='M550 185 L680 185' stroke='#334155' stroke-width='2' marker-end='url(#arrow)'/>"
            + "<defs><marker id='arrow' markerWidth='10' markerHeight='10' refX='9' refY='3' orient='auto' markerUnits='strokeWidth'><path d='M0,0 L0,6 L9,3 z' fill='#334155'/></marker></defs>"
        )
        write_svg("integrative_workflow_overview", 900, 330, workflow)

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch

        def save_all(fig, stem: str, svg: bool = True) -> None:
            exts = ["png", "pdf"] + (["svg"] if svg else [])
            for ext in exts:
                fig.savefig(vis_dir / f"{stem}.{ext}", bbox_inches="tight")
                manifest.append({"figure": f"{stem}.{ext}", "status": "created"})
            plt.close(fig)

        labels = [r["integrative_class"] for r in classes] or ["no_data"]
        values = [as_int(r.get("n_genes", "0")) for r in classes] or [0]
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.bar(labels, values, color=["#2563eb", "#16a34a", "#f59e0b", "#dc2626", "#7c3aed"][: len(labels)])
        ax.set_ylabel("Genes")
        ax.set_title("Integrative classes")
        ax.tick_params(axis="x", labelrotation=30)
        save_all(fig, "barplot_integrative_classes")

        top = scores[: min(30, len(scores))]
        fig, ax = plt.subplots(figsize=(8, max(4, len(top) * 0.25)))
        ax.barh([r.get("gene_id", "") for r in reversed(top)], [as_float(r.get("candidate_score")) for r in reversed(top)], color="#16a34a")
        ax.set_xlabel("Candidate score")
        ax.set_title("Top candidate genes")
        save_all(fig, "top_candidate_scores")

        group_counts = Counter(r.get("machinery_group", "unknown") for r in epi_catalog)
        fig, ax = plt.subplots(figsize=(8, max(4, len(group_counts) * 0.35)))
        groups = list(group_counts.keys()) or ["no_catalog"]
        counts = [group_counts[g] for g in groups] or [0]
        ax.barh(list(reversed(groups)), list(reversed(counts)), color="#0f766e")
        ax.set_xlabel("Genes")
        ax.set_title("Epigenetic machinery catalog")
        save_all(fig, "epigenetic_catalog_groups")

        marks = sorted({r.get("mark_or_factor", "unknown") for r in mark_stage})
        stages = sorted({r.get("stage_or_condition", "unknown") for r in mark_stage})
        matrix = [[0 for _ in stages] for _ in marks]
        mark_index = {m: i for i, m in enumerate(marks)}
        stage_index = {s: i for i, s in enumerate(stages)}
        for row in mark_stage:
            if row.get("mark_or_factor") in mark_index and row.get("stage_or_condition") in stage_index:
                matrix[mark_index[row["mark_or_factor"]]][stage_index[row["stage_or_condition"]]] = as_int(row.get("n_samples", "0"))
        fig, ax = plt.subplots(figsize=(max(6, len(stages) * 1.1), max(4, len(marks) * 0.55)))
        im = ax.imshow(matrix or [[0]], cmap="YlGnBu")
        ax.set_xticks(range(len(stages)))
        ax.set_xticklabels(stages, rotation=30, ha="right")
        ax.set_yticks(range(len(marks)))
        ax.set_yticklabels(marks)
        for i, mark in enumerate(marks):
            for j, _stage in enumerate(stages):
                ax.text(j, i, str(matrix[i][j]), ha="center", va="center", color="#111827")
        ax.set_title("ChIP-seq marks by life-cycle stage")
        fig.colorbar(im, ax=ax, label="Samples")
        save_all(fig, "chip_mark_stage_matrix")

        relation_counts = Counter((r.get("mark_or_factor", "unknown"), r.get("stage_or_condition", "unknown")) for r in gene_mark)
        rel_marks = sorted({m for m, _s in relation_counts})
        rel_stages = sorted({s for _m, s in relation_counts})
        fig, ax = plt.subplots(figsize=(max(6, len(rel_stages) * 1.1), max(4, len(rel_marks) * 0.55)))
        if relation_counts:
            rel_matrix = [[relation_counts[(m, s)] for s in rel_stages] for m in rel_marks]
            im = ax.imshow(rel_matrix, cmap="PuBuGn")
            ax.set_xticks(range(len(rel_stages)))
            ax.set_xticklabels(rel_stages, rotation=30, ha="right")
            ax.set_yticks(range(len(rel_marks)))
            ax.set_yticklabels(rel_marks)
            for i, _mark in enumerate(rel_marks):
                for j, _stage in enumerate(rel_stages):
                    ax.text(j, i, str(rel_matrix[i][j]), ha="center", va="center", color="#111827")
            fig.colorbar(im, ax=ax, label="Linked genes")
        else:
            ax.text(
                0.5,
                0.55,
                "Awaiting server-side ChIP peak annotations",
                ha="center",
                va="center",
                fontsize=14,
                fontweight="bold",
            )
            ax.text(
                0.5,
                0.42,
                "Run on the HPC outputs to populate gene-mark-stage links.",
                ha="center",
                va="center",
                fontsize=10,
            )
            ax.set_axis_off()
        ax.set_title("Gene-mark-stage links")
        save_all(fig, "gene_mark_stage_matrix")

        fig, ax = plt.subplots(figsize=(10, 4.8))
        ax.set_axis_off()
        boxes = [
            (0.04, 0.58, "RNA-seq", "TPM/counts, DEG,\nWGCNA/Mfuzz evidence", "#dbeafe"),
            (0.04, 0.16, "ChIP-seq", "Metadata, peaks,\npeak-gene links", "#dcfce7"),
            (0.38, 0.36, "Integration", "gene x mark x stage\nexpression + chromatin", "#fef3c7"),
            (0.72, 0.36, "Outputs", "candidate genes,\nHTML report, figures", "#ede9fe"),
        ]
        for x, y, title, body, color in boxes:
            patch = FancyBboxPatch((x, y), 0.23, 0.25, boxstyle="round,pad=0.02,rounding_size=0.02", fc=color, ec="#374151", lw=1.2)
            ax.add_patch(patch)
            ax.text(x + 0.115, y + 0.17, title, ha="center", va="center", fontsize=13, fontweight="bold")
            ax.text(x + 0.115, y + 0.08, body, ha="center", va="center", fontsize=9)
        for y in [0.705, 0.285]:
            ax.annotate("", xy=(0.38, 0.485), xytext=(0.27, y), arrowprops=dict(arrowstyle="->", lw=1.6, color="#374151"))
        ax.annotate("", xy=(0.72, 0.485), xytext=(0.61, 0.485), arrowprops=dict(arrowstyle="->", lw=1.6, color="#374151"))
        ax.set_title("Integrative analysis workflow", fontsize=15, fontweight="bold")
        save_all(fig, "integrative_workflow_overview")
    except Exception as exc:  # pragma: no cover - depends on local plotting stack
        generate_fallback_svgs(exc)
    write_table(vis_dir / "visualization_manifest.tsv", manifest, ["figure", "status"])


def command_functional(_args: argparse.Namespace) -> None:
    _h, top = read_table(str(outdir("080-candidate-scoring") / "top_candidates.tsv"))
    f_header, f_rows = read_table(env("FUNCTIONAL_ANNOTATION"))
    if not f_header:
        write_table(outdir("100-functional-analysis") / "functional_enrichment.tsv", [], ["term", "n_selected", "n_background", "note"])
        write_text(outdir("100-functional-analysis") / "functional_analysis_skipped.txt", "No functional annotation file was configured or found.\n")
        return
    gene_col = first_col(f_header, ["gene_id", "gene", "id"])
    term_col = first_col(f_header, ["term", "go", "GO", "kegg", "pathway", "description", "functional_annotation"])
    selected = {r.get("gene_id", "") for r in top if r.get("gene_id")}
    bg_by_term: dict[str, set[str]] = defaultdict(set)
    sel_by_term: dict[str, set[str]] = defaultdict(set)
    for row in f_rows:
        gid = row.get(gene_col, "") if gene_col else ""
        raw_terms = row.get(term_col, "") if term_col else ""
        for term in re.split(r"[;,|]", raw_terms):
            term = term.strip()
            if not gid or not term:
                continue
            bg_by_term[term].add(gid)
            if gid in selected:
                sel_by_term[term].add(gid)
    rows = [
        {
            "term": term,
            "n_selected": len(sel_by_term.get(term, set())),
            "n_background": len(bg),
            "selected_genes": ";".join(sorted(sel_by_term.get(term, set()))),
            "note": "descriptive_count_offline",
        }
        for term, bg in bg_by_term.items()
        if sel_by_term.get(term)
    ]
    rows.sort(key=lambda r: (-r["n_selected"], r["term"]))
    write_table(outdir("100-functional-analysis") / "functional_enrichment.tsv", rows, ["term", "n_selected", "n_background", "selected_genes", "note"])


def command_report(_args: argparse.Namespace) -> None:
    report_dir = outdir("110-reports")
    _h, classes = read_table(str(outdir("070-integrated-tables") / "integrative_class_counts.tsv"))
    _h, top = read_table(str(outdir("080-candidate-scoring") / "top_candidates.tsv"))
    _h, validation = read_table(str(outdir("010-input-validation") / "validation_report.tsv"))
    _h, mark_stage = read_table(str(outdir("060-chipseq-summary") / "chip_mark_stage_metadata.tsv"))
    _h, gene_mark = read_table(str(outdir("070-integrated-tables") / "gene_mark_stage_summary.tsv"))
    _h, epi_catalog = read_table(str(outdir("030-id-harmonization") / "epigenetic_machinery_catalog.tsv"))
    _h, figure_manifest = read_table(str(outdir("090-visualizations") / "visualization_manifest.tsv"))
    n_genes = sum(as_int(r.get("n_genes", "0")) for r in classes)
    generated = dt.datetime.now().isoformat(timespec="seconds")
    md = [
        "# Integrative RNA-seq + ChIP-seq report",
        "",
        f"Generated: {generated}",
        "",
        "## Inputs",
        "",
    ]
    md.extend(f"- {r.get('item')}: {r.get('status')} {r.get('message')}".strip() for r in validation)
    md.extend(
        [
            "",
            "## Summary",
            "",
            f"- Genes in integrated table: {n_genes}",
            f"- Epigenetic machinery genes in catalog: {len(epi_catalog)}",
            f"- ChIP mark/stage metadata combinations: {len(mark_stage)}",
            f"- Gene-mark-stage links from annotated peaks: {len(gene_mark)}",
            "",
        ]
    )
    md.append("### ChIP marks by stage/condition")
    md.append("")
    if mark_stage:
        for r in mark_stage[:50]:
            md.append(f"- {r.get('stage_or_condition')}: {r.get('mark_or_factor')} ({r.get('n_samples')} sample[s])")
    else:
        md.append("- No ChIP metadata rows were available.")
    md.append("")
    md.append("### Integrative classes")
    md.append("")
    for r in classes:
        md.append(f"- {r.get('integrative_class')}: {r.get('n_genes')}")
    md.extend(["", "### Top candidates", ""])
    for r in top[:20]:
        md.append(f"- {r.get('gene_id')} score={r.get('candidate_score')} class={r.get('integrative_class')}")
    md.extend(
        [
            "",
            "### Gene-mark-stage outputs",
            "",
            "- `070-integrated-tables/gene_mark_stage_links.tsv`: one row per peak-gene-mark-stage association.",
            "- `070-integrated-tables/gene_mark_stage_summary.tsv`: collapsed catalog by gene, mark, and stage.",
            "- `070-integrated-tables/mark_to_gene_catalog.tsv`: marks with linked genes and epigenetic machinery counts.",
            "",
            "## Limitations",
            "",
            "- ChIP-seq links are associations between peaks and genes, not proof of causality.",
            "- Concordance depends on mark/factor biology configured in config/chip_marks_config.tsv.",
            "- Offline functional analysis uses supplied annotation only.",
            "",
            "## Key generated files",
            "",
            "- 030-id-harmonization/gene_master_table.tsv",
            "- 040-peak-gene-mapping/peak_to_gene.tsv",
            "- 070-integrated-tables/integrated_gene_table.tsv",
            "- 080-candidate-scoring/candidate_gene_scores.tsv",
            "- 090-visualizations/visualization_manifest.tsv",
            "- 100-functional-analysis/functional_enrichment.tsv",
        ]
    )
    md_text = "\n".join(md) + "\n"
    write_text(report_dir / "integrative_report.md", md_text)

    def table_html(rows: list[dict[str, str]], columns: list[str], limit: int = 12) -> str:
        if not rows:
            return "<p class='muted'>No rows available yet.</p>"
        head = "".join(f"<th>{html.escape(col)}</th>" for col in columns)
        body = []
        for row in rows[:limit]:
            body.append("<tr>" + "".join(f"<td>{html.escape(str(row.get(col, '')))}</td>" for col in columns) + "</tr>")
        more = f"<p class='muted'>Showing {min(limit, len(rows))} of {len(rows)} rows.</p>" if len(rows) > limit else ""
        return f"<div class='table-wrap'><table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>{more}"

    display_figures = []
    seen_figures = set()
    preferred = {}
    for row in figure_manifest:
        fig = row.get("figure", "")
        stem = str(Path(fig).with_suffix("")) if fig else ""
        if not stem:
            continue
        if fig.endswith(".png") or (fig.endswith(".svg") and stem not in preferred):
            preferred[stem] = fig
    for fig in preferred.values():
        if fig not in seen_figures:
            seen_figures.add(fig)
            display_figures.append(fig)
    figure_cards = "".join(
        f"<figure><img src='../090-visualizations/{html.escape(fig)}' alt='{html.escape(fig)}'><figcaption>{html.escape(Path(fig).stem.replace('_', ' '))}</figcaption></figure>"
        for fig in display_figures
    )
    if not figure_cards:
        figure_cards = "<p class='muted'>Run the visualize step to populate figure panels.</p>"

    class_items = "".join(
        f"<div class='metric small'><span>{html.escape(r.get('integrative_class', 'unknown'))}</span><strong>{html.escape(str(r.get('n_genes', '0')))}</strong></div>"
        for r in classes
    )
    html_body = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Integrative RNA-seq + ChIP-seq Report</title>
  <style>
    body {{ margin: 0; font-family: Arial, Helvetica, sans-serif; color: #111827; background: #f8fafc; }}
    header {{ background: #111827; color: white; padding: 28px 36px; }}
    header h1 {{ margin: 0 0 8px; font-size: 28px; }}
    header p {{ margin: 0; color: #cbd5e1; }}
    main {{ max-width: 1180px; margin: 0 auto; padding: 24px; }}
    section {{ margin: 0 0 26px; }}
    h2 {{ font-size: 20px; margin: 0 0 14px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 14px; }}
    .metric, .panel, figure {{ background: white; border: 1px solid #e5e7eb; border-radius: 8px; padding: 16px; box-shadow: 0 1px 2px rgba(15, 23, 42, 0.05); }}
    .metric span {{ display: block; color: #64748b; font-size: 13px; margin-bottom: 8px; }}
    .metric strong {{ font-size: 28px; }}
    .metric.small strong {{ font-size: 22px; }}
    .figures {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(310px, 1fr)); gap: 16px; }}
    figure {{ margin: 0; }}
    figure img {{ width: 100%; height: auto; display: block; }}
    figcaption {{ color: #475569; font-size: 13px; margin-top: 8px; text-transform: capitalize; }}
    .table-wrap {{ overflow-x: auto; border: 1px solid #e5e7eb; border-radius: 8px; background: white; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
    th, td {{ padding: 8px 10px; border-bottom: 1px solid #e5e7eb; text-align: left; white-space: nowrap; }}
    th {{ background: #f1f5f9; font-weight: 700; }}
    .muted {{ color: #64748b; }}
    .status-ok {{ color: #166534; font-weight: 700; }}
    .status-warning {{ color: #92400e; font-weight: 700; }}
    .status-error {{ color: #991b1b; font-weight: 700; }}
    code {{ background: #eef2ff; padding: 2px 5px; border-radius: 4px; }}
  </style>
</head>
<body>
  <header>
    <h1>Integrative RNA-seq + ChIP-seq Report</h1>
    <p>Generated {html.escape(generated)}. Focus: genes, epigenetic marks, life-cycle stages, and expression evidence.</p>
  </header>
  <main>
    <section class="grid">
      <div class="metric"><span>Integrated genes</span><strong>{n_genes}</strong></div>
      <div class="metric"><span>Epigenetic machinery catalog entries</span><strong>{len(epi_catalog)}</strong></div>
      <div class="metric"><span>ChIP mark/stage combinations</span><strong>{len(mark_stage)}</strong></div>
      <div class="metric"><span>Gene-mark-stage links</span><strong>{len(gene_mark)}</strong></div>
    </section>
    <section>
      <h2>Figures</h2>
      <div class="figures">{figure_cards}</div>
    </section>
    <section>
      <h2>Integrative Classes</h2>
      <div class="grid">{class_items}</div>
    </section>
    <section>
      <h2>ChIP Marks By Stage</h2>
      {table_html(mark_stage, ["stage_or_condition", "mark_or_factor", "n_samples", "replicates", "batches"], 30)}
    </section>
    <section>
      <h2>Top Candidate Genes</h2>
      {table_html(top, ["candidate_score", "gene_id", "gene_name", "integrative_class", "is_epigenetic_machinery"], 20)}
    </section>
    <section>
      <h2>Epigenetic Machinery Catalog</h2>
      {table_html(epi_catalog, ["gene_id", "gene_name", "machinery_group", "description", "catalog_source"], 25)}
    </section>
    <section>
      <h2>Gene-Mark-Stage Links</h2>
      {table_html(gene_mark, ["gene_id", "gene_name", "mark_or_factor", "stage_or_condition", "n_peaks", "n_promoter_peaks"], 25)}
      <p class="muted">If this table is empty locally, run the pipeline on the server where ChIP annotated peaks and count outputs are stored.</p>
    </section>
    <section>
      <h2>Input Validation</h2>
      {table_html(validation, ["item", "status", "message", "path"], 40)}
    </section>
    <section class="panel">
      <h2>Key Output Files</h2>
      <p><code>070-integrated-tables/gene_mark_stage_summary.tsv</code> answers the central gene-mark-stage question.</p>
      <p><code>080-candidate-scoring/candidate_gene_scores.tsv</code> ranks potential regulators of parasite plasticity.</p>
      <p><code>030-id-harmonization/epigenetic_machinery_catalog.tsv</code> consolidates the S. mansoni epigenetic machinery catalog.</p>
    </section>
  </main>
</body>
</html>
"""
    write_text(report_dir / "integrative_report.html", html_body)


COMMANDS = {
    "validate": command_validate,
    "prepare": command_prepare,
    "harmonize": command_harmonize,
    "map-peaks": command_map_peaks,
    "summarize-rna": command_summarize_rna,
    "summarize-chip": command_summarize_chip,
    "integrate": command_integrate,
    "score": command_score,
    "visualize": command_visualize,
    "functional": command_functional,
    "report": command_report,
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=sorted(COMMANDS))
    args = parser.parse_args()
    COMMANDS[args.command](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
