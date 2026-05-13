#!/usr/bin/env python3
"""Collect and summarize ichorCNA results across samples.

Produces:
  - ichorCNA_summary.tsv   — one row per sample with key metrics
  - ichorCNA_report.html   — self-contained HTML report with embedded figures
"""

import sys
import argparse
import base64
from pathlib import Path

import jinja2
import pandas as pd


# ---------------------------------------------------------------------------
# ichorCNA event vocabulary
# ---------------------------------------------------------------------------

# Human-readable expansion of ichorCNA event codes (source: ichorCNA docs)
EVENT_LABELS: dict[str, str] = {
    "HOMD": "homozygous_deletion",
    "DLOH": "deletion_loh",
    "NLOH": "neutral_loh",
    "HETD": "heterozygous_deletion",
    "NEUT": "neutral",
    "GAIN": "gain",
    "AMP": "amplification",
    "HLAMP": "high_level_amplification",
    "ASCNA": "allele_specific_cna",
    "UBCNA": "unbalanced_cna",
    "BCNA": "balanced_cna",
}

ALTERED_EVENTS: frozenset[str] = frozenset(k for k in EVENT_LABELS if k != "NEUT")

# Suffixes (after the sample-ID prefix) that identify the key ichorCNA figures.
# Per-chromosome, per-solution, and bias plots are excluded from the report.
IMPORTANT_PDF_SUFFIXES: frozenset[str] = frozenset({
    "_genomeWide",
    "_genomeWideCorrection",
    "_genomeWide_all_sols",
    "_tpdf",
})

# Display order within a sample's figures section
PDF_PRIORITY: list[str] = [
    "_genomeWide",
    "_tpdf",
    "_genomeWideCorrection",
    "_genomeWide_all_sols",
]


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def parse_params(path: Path) -> dict:
    """Parse one ichorCNA params.txt.

    The file has a mixed format: a small TSV summary at the top, then a
    free-text 'Key:    value' block, then a tabular solution grid at the
    bottom.  All useful metrics live in the free-text block, so we scan
    that directly instead of using pandas.
    """
    kv_wanted: dict[str, str] = {
        "Tumor Fraction":        "tumor_fraction",
        "Ploidy":                "ploidy",
        "Subclone Fraction":     "subclone_fraction",
        "Coverage":              "coverage_median",
        "GC-Map correction MAD": "mad_cor",
        "Log-likelihood":        "loglik",
    }
    record: dict = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" not in line:
            continue
        key_part, _, val_part = line.partition(":")
        key = key_part.strip()
        if key not in kv_wanted:
            continue
        # Take first whitespace-delimited token; skips "NA", multi-value lines, etc.
        tokens = val_part.strip().split()
        if not tokens:
            continue
        try:
            record[kv_wanted[key]] = float(tokens[0])
        except ValueError:
            pass  # NA or non-numeric — field simply omitted
    return record


# QC thresholds from ichorCNA FAQ:
#   PASS  — MAD < 0.15 (high quality)
#   WARN  — 0.15 ≤ MAD ≤ 0.30 (usable but elevated noise)
#   FAIL  — MAD > 0.30 (likely too noisy)
def qc_status(mad_cor: float | None) -> str:
    if mad_cor is None:
        return "N/A"
    if mad_cor < 0.15:
        return "PASS"
    if mad_cor <= 0.30:
        return "WARN"
    return "FAIL"


def fraction_genome_altered(path: Path) -> float | None:
    """Fraction of genomic bases bearing a copy-number alteration."""
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()

    start_col = _find_col(df, ["start", "loc.start", "Start"])
    end_col   = _find_col(df, ["end",   "loc.end",   "End"])
    event_col = _find_col(df, ["Corrected_Call", "event", "call", "Event", "Call"])
    cn_col    = _find_col(df, ["Corrected_Copy_Number", "copy.number", "copy_number", "CN", "cn"])

    if start_col is None or end_col is None:
        return None

    sizes    = df[end_col] - df[start_col]
    total_bp = sizes.sum()
    if total_bp == 0:
        return None

    if event_col is not None:
        altered_mask = df[event_col].isin(ALTERED_EVENTS)
    elif cn_col is not None:
        altered_mask = df[cn_col] != 2
    else:
        return None

    return float(sizes[altered_mask].sum()) / float(total_bp)


def segment_counts(path: Path) -> dict:
    """Count segments per event type, keyed by human-readable names."""
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()

    event_col = _find_col(df, ["Corrected_Call", "event", "call", "Event", "Call"])
    if event_col is None:
        return {}

    result = {}
    for code, n in df[event_col].value_counts().items():
        label = EVENT_LABELS.get(str(code).upper(), str(code).lower())
        result[f"n_seg_{label}"] = int(n)
    return result


# ---------------------------------------------------------------------------
# File-matching utilities
# ---------------------------------------------------------------------------

def strip_suffix(name: str, suffixes: list[str]) -> str:
    for s in suffixes:
        if name.endswith(s):
            return name[: -len(s)]
    return name


def pdf_sort_key(path: Path) -> int:
    stem = path.stem
    for i, suffix in enumerate(PDF_PRIORITY):
        if stem.endswith(suffix):
            return i
    return len(PDF_PRIORITY)


def index_pdfs_by_sample(
    pdf_paths: list[Path], sample_ids: list[str]
) -> dict[str, list[Path]]:
    """Match each PDF to a sample ID by longest-prefix match; keep only important figures."""
    result: dict[str, list[Path]] = {sid: [] for sid in sample_ids}
    sorted_ids = sorted(sample_ids, key=len, reverse=True)
    for pdf_path in pdf_paths:
        stem = pdf_path.stem
        matched = next((sid for sid in sorted_ids if stem.startswith(sid)), None)
        if matched and stem[len(matched):] in IMPORTANT_PDF_SUFFIXES:
            result[matched].append(pdf_path)
    return result


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

def _b64_pdf(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode("ascii")


def render_html_report(
    summary_df: pd.DataFrame,
    records: dict[str, dict],
    pdfs_by_sample: dict[str, list[Path]],
    template_path: Path,
) -> str:
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(template_path.parent)),
        autoescape=jinja2.select_autoescape(["html"]),
    )
    tmpl = env.get_template(template_path.name)

    # Build per-sample data for the template
    samples_data = []
    for _, row in summary_df.iterrows():
        sid = row["sample"]
        pdfs = sorted(pdfs_by_sample.get(sid, []), key=pdf_sort_key)
        samples_data.append({
            "id": sid,
            "metrics": records[sid],
            "figures": [
                {"label": p.stem, "b64": _b64_pdf(p)}
                for p in pdfs
            ],
        })

    # Column display config passed into the template
    metric_display = {
        "qc_status":            "QC Status",
        "mad_cor":              "GC-Map MAD",
        "tumor_fraction":       "Tumor Fraction",
        "ploidy":               "Ploidy",
        "subclone_fraction":    "Subclone Fraction",
        "coverage_median":      "Coverage (Median)",
        "loglik":               "Log-Likelihood",
        "fraction_genome_altered": "Fraction Genome Altered",
    }

    return tmpl.render(
        summary_columns=list(summary_df.columns),
        summary_rows=summary_df.to_dict(orient="records"),
        samples=samples_data,
        metric_display=metric_display,
        n_samples=len(summary_df),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize ichorCNA results and produce a TSV + self-contained HTML report."
    )
    parser.add_argument("--params", nargs="+", required=True, metavar="PARAMS_TXT",
                        help="*.params.txt files from ichorCNA")
    parser.add_argument("--segs",   nargs="*", default=[],  metavar="SEG_TXT",
                        help="*.seg.txt files (matched to samples by filename)")
    parser.add_argument("--pdfs",   nargs="*", default=[],  metavar="PDF",
                        help="PDF figure files from ichorCNA (matched by filename prefix)")
    parser.add_argument("--output-tsv",  default="ichorCNA_summary.tsv")
    parser.add_argument("--output-html", default="ichorCNA_report.html")
    args = parser.parse_args()

    template_path = Path(__file__).parent / "ichorCNA_report.html.j2"

    seg_by_sample: dict[str, Path] = {
        strip_suffix(Path(p).name, [".seg.txt", ".seg"]): Path(p)
        for p in args.segs
    }

    rows: list[dict] = []
    records: dict[str, dict] = {}

    for params_str in args.params:
        params_path = Path(params_str)
        sample_id   = strip_suffix(params_path.name, [".params.txt"])
        record: dict = {"sample": sample_id}

        try:
            record.update(parse_params(params_path))
        except Exception as exc:
            print(f"WARNING: could not parse {params_path}: {exc}", file=sys.stderr)
        record["qc_status"] = qc_status(record.get("mad_cor"))

        seg_path = seg_by_sample.get(sample_id)
        if seg_path is not None:
            try:
                fga = fraction_genome_altered(seg_path)
                if fga is not None:
                    record["fraction_genome_altered"] = round(fga, 6)
                record.update(segment_counts(seg_path))
            except Exception as exc:
                print(f"WARNING: could not parse {seg_path}: {exc}", file=sys.stderr)

        rows.append(record)
        records[sample_id] = record

    summary = pd.DataFrame(rows)
    fixed = [
        "sample", "qc_status", "mad_cor", "tumor_fraction", "ploidy",
        "subclone_fraction", "coverage_median", "loglik",
        "fraction_genome_altered",
    ]
    seg_cols = sorted(c for c in summary.columns if c.startswith("n_seg_"))
    summary  = summary[[c for c in fixed if c in summary.columns] + seg_cols]

    summary.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Wrote TSV summary for {len(rows)} sample(s) to {args.output_tsv}")

    pdfs_by_sample = index_pdfs_by_sample(
        [Path(p) for p in args.pdfs], list(records.keys())
    )
    html_content = render_html_report(summary, records, pdfs_by_sample, template_path)
    Path(args.output_html).write_text(html_content, encoding="utf-8")
    print(f"Wrote HTML report to {args.output_html}")


if __name__ == "__main__":
    main()