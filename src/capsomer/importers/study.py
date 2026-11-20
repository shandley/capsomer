"""Start a bundle from a study's sample table, and record read quality control."""

import json
from pathlib import Path
from typing import Any, cast

from capsomer.builder import Builder
from capsomer.bundle import quote
from capsomer.importers.common import tsv


def init_bundle(
    bundle_dir: Path,
    sample_table: Path,
    *,
    bundle_id: str | None = None,
    title: str | None = None,
) -> int:
    """Create a bundle holding the samples and libraries of a study.

    The table is tab-separated with one row per library. It needs a sample_id
    column. library_id defaults to sample_id and control_type to not_recorded.
    Columns named as the catalog names them go to the sample or library table.
    Every other column is kept on the sample as study metadata.
    """
    with Builder(bundle_dir, bundle_id=bundle_id, title=title) as builder:
        builder.con.execute(f"CREATE TEMP TABLE _study AS SELECT * FROM {tsv(sample_table)}")
        found = [str(r[0]) for r in builder.con.execute("DESCRIBE _study").fetchall()]
        if "sample_id" not in found:
            raise ValueError(f"{sample_table} has no sample_id column")
        library_columns = builder.catalog.tables["library"].column_names
        sample_columns = builder.catalog.tables["sample"].column_names
        on_library = [c for c in found if c in library_columns and c != "sample_id"]
        metadata = [c for c in found if c not in library_columns and c not in sample_columns]
        control = "control_type" if "control_type" in found else "'not_recorded'"
        on_sample = [c for c in found if c in sample_columns and c != "control_type"] + metadata
        builder.add(
            "sample",
            f"SELECT DISTINCT {', '.join(quote(c) for c in on_sample)}, "
            f"coalesce({control}, 'not_recorded') AS control_type FROM _study",
        )
        library = "library_id" if "library_id" in found else "sample_id"
        kept = [quote(c) for c in on_library if c != "library_id"]
        builder.add(
            "library",
            f"SELECT {library} AS library_id, sample_id{''.join(', ' + c for c in kept)} FROM _study",
        )
        builder.describe(
            "sample", {c: f"Study metadata: column {c} of the sample table." for c in metadata}
        )
        return int(builder.scalar("SELECT count(*) FROM library"))


def add_fastp(
    builder: Builder,
    report: Path,
    *,
    run_id: str,
    library_id: str,
    parent: str | None = None,
) -> None:
    """Record a fastp run and the read totals it reports for one library."""
    data = cast(dict[str, Any], json.loads(report.read_text()))
    summary = cast(dict[str, Any], data["summary"])
    before = int(summary["before_filtering"]["total_reads"])
    after = int(summary["after_filtering"]["total_reads"])
    paired = "paired" in str(summary.get("sequencing", ""))
    builder.add_run(
        run_id,
        tool="fastp",
        tool_version=str(summary.get("fastp_version") or "") or None,
        step="read_qc",
        description="Adapter and quality trimming of reads.",
        parent=parent,
        parameters={"command": data.get("command")},
        libraries=[library_id],
    )
    per_fragment = 2 if paired else 1
    builder.con.execute(
        "UPDATE library SET reads_raw = coalesce(reads_raw, ?), "
        "fragments_raw = coalesce(fragments_raw, ?), reads_after_qc = ?, fragments_after_qc = ? "
        "WHERE library_id = ?",
        [before, before // per_fragment, after, after // per_fragment, library_id],
    )
