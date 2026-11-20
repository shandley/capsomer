"""Import a tiny ViroForge-style dataset and check the truth bundle."""

import json
from pathlib import Path

from capsomer.bundle import connect
from capsomer.importers.viroforge import import_viroforge
from capsomer.validate import conforms, validate

VIRUS = "ACGTACGTACGTAAACGGTTACGATCGA"
HOST = "TTTTGGGGCCCCAAAATTGGCCAATTGG"


def fastq(names: list[str], bases: str) -> str:
    return "".join(f"@{name}\n{bases[:10]}\n+\n{'I' * 10}\n" for name in names)


def make_dataset(folder: Path) -> Path:
    for sub in ("fasta", "fastq", "metadata"):
        (folder / sub).mkdir(parents=True)
    # Two contaminant records share one sequence, as real ViroForge output does.
    (folder / "fasta" / "demo.fasta").write_text(
        f">GCF_000000001.1 Demo virus | Microviridae | abundance=0.9\n{VIRUS}\n"
        f">host_human_0000 host\n{HOST}\n>host_human_0001 host\n{HOST}\n"
    )
    r1 = [
        "GCF_000000001.1_0_0/1 source=viral",
        "GCF_000000001.1_1_0/1 source=viral adapter_injected=truseq:24bp",
        "GCF_000000001.1_1_0/1_dup1 source=viral",
        "host_human_0001_0_0/1 source=host",
    ]
    (folder / "fastq" / "demo_R1.fastq").write_text(fastq(r1, VIRUS))
    (folder / "fastq" / "demo_R2.fastq").write_text(
        fastq([n.replace("/1", "/2") for n in r1], HOST)
    )
    (folder / "metadata" / "demo_duplicate_manifest.tsv").write_text(
        "template_read_id\tcopy_read_id\tcopy_number\ttotal_copies\n"
        "GCF_000000001.1_1_0/1\tGCF_000000001.1_1_0/1_dup1\t1\t1\n"
        "GCF_000000001.1_1_0/2\tGCF_000000001.1_1_0/2_dup1\t1\t1\n"
    )
    metadata = {
        "generation_info": {"viroforge_version": "0.0.0", "random_seed": 1},
        "collection": {"name": "Demo collection", "description": "Two sequences."},
        "configuration": {
            "coverage": 5.0,
            "read_length": 10,
            "platform": "novaseq",
            "molecule_type": "dna",
            "vlp_protocol": "none",
            "amplification": "none",
        },
        "sequences": [
            {"genome_id": "GCF_000000001.1", "genome_name": "Demo virus", "sequence_type": "viral"},
            {"genome_id": "host_human_0000", "genome_name": "host", "sequence_type": "contaminant"},
        ],
        "benchmarking": {
            "taxonomy": {
                "GCF_000000001.1": {
                    "ncbi_taxid": 12345,
                    "realm": "Monodnaviria",
                    "family": "Microviridae",
                    "genus": None,
                    "species": "Demo virus",
                    "is_known": True,
                }
            },
            "contamination_manifest": {
                "contaminants": [
                    {
                        "genome_id": "host_human_0000",
                        "contaminant_type": "host_dna",
                        "organism": "Homo sapiens",
                    },
                    {
                        "genome_id": "host_human_0001",
                        "contaminant_type": "host_dna",
                        "organism": "Homo sapiens",
                    },
                ]
            },
            "expected_coverage": {
                "per_genome": [
                    {
                        "genome_id": "GCF_000000001.1",
                        "expected_coverage": 5.0,
                        "relative_abundance": 0.9,
                    }
                ]
            },
        },
    }
    (folder / "metadata" / "demo_metadata.json").write_text(json.dumps(metadata))
    return folder


def test_truth_bundle_conforms_and_joins(tmp_path: Path) -> None:
    dataset = make_dataset(tmp_path / "dataset")
    bundle = tmp_path / "truth" / "bundle"
    summary = import_viroforge(dataset, bundle, sample_id="demo")
    assert summary.reads_without_source == 0
    findings = validate(bundle)
    assert conforms(findings), findings
    con = connect(bundle)

    def one(sql: str) -> object:
        row = con.execute(sql).fetchone()
        assert row is not None
        return row[0]

    # Three records, two distinct sequences.
    assert one("SELECT count(*) FROM sequence") == 2
    assert one("SELECT count(*) FROM x_sequence_origin") == 3
    assert one("SELECT fragments_raw FROM library") == 4
    assert (
        one("SELECT native_id FROM taxon WHERE name = 'Demo virus' AND native_id IS NOT NULL")
        == "12345"
    )
    # Every read joins to the family of the genome it was simulated from.
    assert (
        one(
            "SELECT count(*) FROM alignment a JOIN genome_member m ON m.sequence_id = a.target_sequence_id "
            "JOIN taxon_assignment ta USING (genome_id) JOIN taxon t USING (taxon_id) "
            "WHERE t.family_name = 'Microviridae'"
        )
        == 6
    )
    assert one("SELECT direct_count FROM taxon_profile") == 3
    artifacts = con.execute(
        "SELECT artifact, count(*) FROM x_read_artifact GROUP BY 1 ORDER BY 1"
    ).fetchall()
    assert artifacts == [("adapter", 2), ("pcr_duplicate", 2)]
    assert one(
        "SELECT detail FROM x_read_artifact WHERE artifact = 'pcr_duplicate' ORDER BY read_id LIMIT 1"
    ) == ("copy of GCF_000000001.1_1_0/1")
