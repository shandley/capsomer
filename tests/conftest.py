"""A minimal conforming bundle, built once per test."""

from pathlib import Path

import duckdb
import pytest

from capsomer import ids
from capsomer.bundle import build_manifest, write_table

READ = "ACGTACGTAC"
CONTIG = "ACGTACGTACGGTTAACC"


def minimal_tables() -> dict[str, str]:
    read, contig = ids.sequence_id(READ), ids.sequence_id(CONTIG)
    return {
        "sample": "SELECT 's1' AS sample_id, 'not_control' AS control_type, 'vaccine' AS treatment",
        "library": "SELECT 'l1' AS library_id, 's1' AS sample_id, 5 AS fragments_after_qc",
        "run": (
            "SELECT * FROM (VALUES ('p', 'demo', 'pipeline', NULL), ('p.cluster', 'demo', "
            "'read_clustering', 'p'), ('p.search', 'demo', 'search', 'p'), ('p.map', 'demo', "
            "'read_mapping', 'p')) AS t(run_id, tool, step, parent_run_id)"
        ),
        "run_parameter": (
            "SELECT * FROM (VALUES ('p.cluster', 'identity', '1.0'), ('p.search', 'hits_kept', "
            "'best'), ('p.map', 'mapped_unit', 'read clusters')) AS t(run_id, name, value)"
        ),
        "reference_release": (
            "SELECT 'tax' AS reference_release_id, 'demo taxonomy' AS name, "
            "'formal_taxonomy' AS kind, 'abc' AS checksum"
        ),
        "taxon": "SELECT 'tax:1' AS taxon_id, 'tax' AS reference_release_id, 'Demovirus' AS name",
        "sequence": (
            f"SELECT * FROM (VALUES ('{read}', 'read_cluster', 10, '{READ}'), "
            f"('{contig}', 'contig', 18, '{CONTIG}')) AS t(sequence_id, kind, length, bases)"
        ),
        "read_count": (
            f"SELECT '{read}' AS sequence_id, 'l1' AS library_id, 'p.cluster' AS run_id, "
            "'single' AS read_form, 'not_applicable' AS read_end, 5 AS count, "
            "'fragments' AS count_unit"
        ),
        "taxon_assignment": (
            f"SELECT * FROM (VALUES ('{read}', 'p.search', 'assigned', 'tax:1', true), "
            f"('{contig}', 'p.search', 'none', NULL, false)) "
            "AS t(sequence_id, run_id, outcome, taxon_id, is_preferred)"
        ),
        "alignment": (
            f"SELECT '{read}' AS query_sequence_id, '{contig}' AS target_sequence_id, "
            "'p.map' AS run_id, 'nucleotide_to_nucleotide' AS search_type, 1 AS target_start, "
            "10 AS target_end, true AS is_primary"
        ),
    }


def write_bundle(out_dir: Path, tables: dict[str, str]) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    for name, query in tables.items():
        write_table(con, name, query, out_dir)
    build_manifest(
        out_dir,
        bundle_id="demo",
        title="Demo bundle",
        extra_columns={"sample": {"treatment": "Treatment group."}},
    )
    return out_dir


@pytest.fixture
def bundle(tmp_path: Path) -> Path:
    return write_bundle(tmp_path / "bundle", minimal_tables())
