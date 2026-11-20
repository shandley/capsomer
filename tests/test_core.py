import json
import subprocess
import sys
from importlib import resources
from pathlib import Path

from capsomer import ids, load_catalog
from capsomer.schema import SHAPES
from capsomer.validate import conforms, validate
from tests.conftest import minimal_tables, write_bundle


def test_sequence_id_matches_ga4gh_vector() -> None:
    assert ids.sequence_id("ACGT") == "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
    assert ids.sequence_id("acgt") == ids.sequence_id("ACGT")


def test_no_strand_rule() -> None:
    assert ids.sequence_id("AAC") != ids.sequence_id("GTT")


def test_read_pair_id_covers_both_reads_in_order() -> None:
    pair = ids.read_pair_id("AAC", "GGT")
    assert pair.startswith("RP.")
    assert pair != ids.read_pair_id("GGT", "AAC")
    assert pair != ids.read_pair_id("AAC", "GGA")


def test_genome_id_ignores_member_order() -> None:
    a, b = ids.sequence_id("AAC"), ids.sequence_id("GGT")
    assert ids.genome_id([a, b]) == ids.genome_id([b, a])
    assert ids.genome_id([a]).startswith("GN.")


def test_catalog_shapes_and_references() -> None:
    catalog = load_catalog()
    assert len(catalog.tables) == 22
    required = {t.name for t in catalog.tables.values() if t.catalog_required}
    assert required == {
        "sample",
        "library",
        "run",
        "reference_release",
        "sequence",
        "taxon",
        "taxon_assignment",
    }
    for table in catalog.tables.values():
        assert table.shape in SHAPES
        assert table.row
        for column in table.columns:
            assert column.description, f"{table.name}.{column.name} has no description"
            if column.references:
                assert column.references in catalog.tables
                assert column.name.endswith(f"{column.references}_id")
        for group in table.exactly_one:
            assert all(table.column(c) for c in group)
        assert all(table.column(k) for k in table.key)


def test_schema_is_valid_linkml() -> None:
    path = resources.files("capsomer").joinpath("schema/capsomer.yaml")
    done = subprocess.run(
        [sys.executable, "-m", "linkml.generators.jsonschemagen", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert done.returncode == 0, done.stderr[-2000:]
    assert "ReadCount" in json.loads(done.stdout)["$defs"]


def test_minimal_bundle_conforms(bundle: Path) -> None:
    findings = validate(bundle)
    assert conforms(findings), findings
    manifest = json.loads((bundle / "manifest.json").read_text())
    assert manifest["tables"]["read_count"]["shape"] == "count"
    assert manifest["tables"]["sample"]["columns"]["treatment"]["extra"] is True
    assert (
        manifest["tables"]["alignment"]["columns"]["target_sequence_id"]["references"] == "sequence"
    )


def checks(findings: list[object]) -> set[tuple[str, str]]:
    return {(f.severity, f.check) for f in findings}  # type: ignore[attr-defined]


def test_broken_bundles_are_caught(tmp_path: Path) -> None:
    tables = minimal_tables()
    tables["library"] = "SELECT 'l1' AS library_id, 'missing' AS sample_id, 3 AS fragments_after_qc"
    found = checks(list(validate(write_bundle(tmp_path / "a", tables))))
    assert ("error", "reference_exists") in found
    assert ("error", "counts_within_total") in found

    tables = minimal_tables()
    tables["sample"] = (
        "SELECT * FROM (VALUES ('s1', 'not_control'), ('s1', 'bogus')) t(sample_id, control_type)"
    )
    found = checks(list(validate(write_bundle(tmp_path / "b", tables))))
    assert ("critical", "identifier_unique") in found
    assert ("error", "allowed_value") in found

    tables = minimal_tables()
    del tables["taxon"]
    tables["sequence"] = tables["sequence"].replace("ACGTACGTAC'", "ACGTACGTAA'", 1)
    found = checks(list(validate(write_bundle(tmp_path / "c", tables))))
    assert ("critical", "required_table") in found
    assert ("error", "identifier_matches_bases") in found
