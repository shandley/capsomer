"""Importers for the contig side: assembly, geNomad, CheckV, clusters, coverage, mappings.

The fixture files copy the column headers of real tool output (geNomad 1.13,
CheckV 1.1, vclust 1.3, MMseqs2, CoverM 0.8). The values are made up.
"""

import random
from pathlib import Path

import pytest

from capsomer import ids
from capsomer.builder import Builder
from capsomer.bundle import connect
from capsomer.importers.checkv import add_checkv
from capsomer.importers.clusters import add_clusters
from capsomer.importers.contigs import add_contigs
from capsomer.importers.genomad import add_genomad
from capsomer.importers.mapping import add_coverm, add_paf, paf_records, sam_records
from capsomer.importers.study import add_fastp, init_bundle
from capsomer.validate import conforms, validate

RNG = random.Random(7)
SEQ = {name: "".join(RNG.choice("ACGT") for _ in range(n)) for name, n in (("a", 60), ("b", 80), ("c", 50), ("d", 40))}  # fmt: skip


def tsv(path: Path, header: str, *rows: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [line.replace(" ", "\t").replace("~", " ") for line in (header, *rows)]
    path.write_text("\n".join(lines) + "\n")
    return path


def fasta(path: Path, records: dict[str, str]) -> Path:
    path.write_text("".join(f">{name}\n{bases}\n" for name, bases in records.items()))
    return path


@pytest.fixture
def tools(tmp_path: Path) -> dict[str, Path]:
    """Made-up outputs of one small analysis: two assemblies, geNomad, CheckV, clusters, coverage."""
    out: dict[str, Path] = {}
    out["samples"] = tsv(tmp_path / "samples.tsv", "sample_id layout body_site", "s1 paired gut")
    out["first"] = fasta(
        tmp_path / "first.fa",
        {
            "k141_1 flag=1 multi=12.5 len=60": SEQ["a"],
            "k141_2 flag=0 multi=3.0 len=80": SEQ["b"],
            "k141_3 flag=0 multi=2.0 len=50": SEQ["c"],
        },
    )
    out["second"] = fasta(
        tmp_path / "second.fa",
        {"NODE_1_length_60_cov_9.5": SEQ["a"], "NODE_2_length_40_cov_4.0": SEQ["d"]},
    )
    g = tmp_path / "genomad"
    tsv(
        g / "x_aggregated_classification/x_aggregated_classification.tsv",
        "seq_name chromosome_score plasmid_score virus_score",
        "k141_1 0.02 0.03 0.95",
        "k141_2 0.80 0.10 0.10",
        "k141_3 0.10 0.85 0.05",
    )
    tsv(
        g / "x_summary/x_virus_summary.tsv",
        "seq_name length topology coordinates n_genes genetic_code virus_score fdr n_hallmarks marker_enrichment taxonomy",
        "k141_1 60 DTR NA 3 11 0.95 NA 2 4.5 Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;",
        "k141_2|provirus_11_40 30 Provirus 11-40 2 11 0.91 NA 1 2.0 Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;",
    )
    tsv(
        g / "x_summary/x_plasmid_summary.tsv",
        "seq_name length topology n_genes genetic_code plasmid_score fdr n_hallmarks marker_enrichment conjugation_genes amr_genes",
        "k141_3 50 No~terminal~repeats 2 11 0.85 NA 1 1.5 NA NA",
    )
    (g / "x_summary/x_summary.json").write_text(
        '{"module": "summary", "parameters": {"min_score": 0.7}}'
    )
    tsv(
        g / "x_annotate/x_taxonomy.tsv",
        "seq_name n_genes_with_taxonomy agreement taxid lineage",
        "k141_1 3 1.0000 1246 Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;",
    )
    tsv(
        g / "x_find_proviruses/x_provirus.tsv",
        "seq_name source_seq start end length n_genes v_vs_c_score in_seq_edge integrases",
        "k141_2|provirus_11_40 k141_2 11 40 30 2 5.0 False NA",
    )
    tsv(
        g / "x_find_proviruses/x_provirus_taxonomy.tsv",
        "seq_name n_genes_with_taxonomy agreement taxid lineage",
        "k141_2|provirus_11_40 2 1.0000 12762 Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;Herelleviridae",
    )
    out["genomad"] = g
    c = tmp_path / "checkv"
    tsv(
        c / "quality_summary.tsv",
        "contig_id contig_length provirus proviral_length gene_count viral_genes host_genes checkv_quality miuvig_quality completeness completeness_method contamination kmer_freq warnings",
        "k141_1 60 No NA 3 2 0 Complete High-quality 100.0 DTR~(high-confidence) 0.0 1.0 ",
    )  # fmt: skip
    tsv(
        c / "completeness.tsv",
        "contig_id contig_length viral_length aai_expected_length",
        "k141_1 60 NA 61.0",
    )
    tsv(
        c / "complete_genomes.tsv",
        "contig_id contig_length kmer_freq prediction_type confidence_level",
        "k141_1 60 1.0 DTR high",
    )
    out["checkv"] = c
    out["viruses"] = fasta(tmp_path / "viruses.fa", {"k141_1": SEQ["a"]})
    out["pooled"] = fasta(
        tmp_path / "pooled.fa", {"k141_1": SEQ["a"], "NODE_1": SEQ["a"], "NODE_2": SEQ["d"]}
    )
    out["vclust"] = tsv(
        tmp_path / "clusters.tsv", "object cluster", "k141_1 k141_1", "NODE_1 k141_1", "NODE_2 NODE_2"
    )  # fmt: skip
    out["mmseqs"] = tmp_path / "mmseqs.tsv"
    out["mmseqs"].write_text("k141_1\tk141_1\nk141_1\tNODE_1\nk141_1\tNODE_2\n")
    out["coverm"] = tmp_path / "coverm.tsv"
    out["coverm"].write_text(
        "Contig\tmap Read Count\tmap Mean\tmap Covered Fraction\tmap Covered Bases\n"
        "k141_1\t40\t12.5\t1\t60\nNODE_2\t0\t0\t0\t0\n"
    )
    out["fastp"] = tmp_path / "fastp.json"
    out["fastp"].write_text(
        '{"summary": {"fastp_version": "1.0", "sequencing": "paired end (100 cycles + 100 cycles)", '
        '"before_filtering": {"total_reads": 200}, "after_filtering": {"total_reads": 180}}, '
        '"command": "fastp -i a -I b"}'
    )
    out["paf"] = tmp_path / "x.paf"
    out["paf"].write_text("NODE_2\t40\t0\t40\t-\tk141_1\t60\t10\t50\t38\t40\t60\ttp:A:P\n")
    return out


@pytest.fixture
def bundle(tmp_path: Path, tools: dict[str, Path]) -> Path:
    """A bundle built with no pipeline in particular: a sample table, then one tool at a time."""
    out = tmp_path / "study" / "bundle"
    init_bundle(out, tools["samples"], title="Test study")
    with Builder(out) as b:
        add_fastp(b, tools["fastp"], run_id="qc", library_id="s1")
        add_contigs(b, tools["first"], run_id="asm1", tool="megahit", libraries=["s1"], inputs=["qc"])  # fmt: skip
        add_contigs(b, tools["second"], run_id="asm2", tool="metaspades", libraries=["s1"], inputs=["qc"])  # fmt: skip
    with Builder(out) as b:
        add_genomad(b, tools["genomad"], tools["first"], run_id="genomad", database_version="1.9", inputs=["asm1"])  # fmt: skip
        add_checkv(b, tools["checkv"], tools["viruses"], run_id="checkv", database_version="1.5", inputs=["genomad"])  # fmt: skip
        add_clusters(b, tools["vclust"], tools["pooled"], run_id="votu", tool="vclust", level="species", parameters={"min_ani": 0.95})  # fmt: skip
        add_clusters(b, tools["mmseqs"], tools["pooled"], run_id="loose", tool="mmseqs", level="genus", layout="mmseqs", parameters={"min_seq_id": 0.7})  # fmt: skip
        add_coverm(b, tools["coverm"], tools["pooled"], run_id="cov", libraries={"map": "s1"}, parameters={"mapped_unit": "reads"})  # fmt: skip
        add_paf(b, tools["paf"], tools["pooled"], tools["pooled"], run_id="aln", tool="minimap2", parameters={"hits_kept": "all"})  # fmt: skip
    return out


def test_bundle_built_tool_by_tool_conforms(bundle: Path) -> None:
    findings = validate(bundle)
    assert conforms(findings), findings
    con = connect(bundle)
    assert con.execute("SELECT body_site, control_type FROM sample").fetchall() == [
        ("gut", "not_recorded")
    ]
    assert con.execute("SELECT reads_raw, fragments_after_qc FROM library").fetchall() == [(200, 90)]  # fmt: skip


def test_identical_contigs_from_two_assemblies_are_one_sequence(bundle: Path) -> None:
    con = connect(bundle)
    shared = ids.sequence_id(SEQ["a"])
    sources = con.execute(
        "SELECT run_id, name, assembler_depth FROM sequence_source WHERE sequence_id = ? ORDER BY 1",
        [shared],
    ).fetchall()
    assert sources == [("asm1", "k141_1", 12.5), ("asm2", "NODE_1_length_60_cov_9.5", 9.5)]
    assert con.execute("SELECT count(*) FROM sequence WHERE sequence_id = ?", [shared]).fetchall() == [(1,)]  # fmt: skip
    assert con.execute("SELECT run_id, library_id FROM run_library WHERE run_id = 'asm2'").fetchall() == [("asm2", "s1")]  # fmt: skip


def test_genomad_calls_scores_and_taxa(bundle: Path) -> None:
    con = connect(bundle)
    calls = dict(
        con.execute(
            "SELECT s.name, v.call FROM viral_call v JOIN sequence_source s USING (sequence_id) "
            "WHERE s.run_id IN ('asm1', 'genomad')"
        ).fetchall()
    )
    assert calls == {
        "k141_1": "virus",
        "k141_2": "none",
        "k141_3": "plasmid",
        "k141_2|provirus_11_40": "virus",
    }
    none = con.execute(
        "SELECT cellular_score, virus_score FROM viral_call WHERE sequence_id = ?",
        [ids.sequence_id(SEQ["b"])],
    ).fetchall()
    assert none == [(0.8, 0.1)]
    taxa = con.execute(
        "SELECT t.name, t.rank, t.native_id, a.outcome FROM taxon_assignment a "
        "LEFT JOIN taxon t USING (taxon_id) WHERE a.run_id = 'genomad' ORDER BY 1 NULLS LAST"
    ).fetchall()
    assert taxa[:2] == [
        ("Caudoviricetes", "class", "1246", "assigned"),
        ("Herelleviridae", "family", "12762", "assigned"),
    ]
    assert [row[3] for row in taxa[2:]] == ["none", "none"]
    assert con.execute(
        "SELECT value FROM run_parameter WHERE run_id = 'genomad' AND name = 'summary.min_score'"
    ).fetchall() == [("0.7",)]


def test_provirus_is_a_sequence_cut_from_its_contig(bundle: Path) -> None:
    con = connect(bundle)
    region = ids.sequence_id(SEQ["b"][10:40])
    assert con.execute(
        "SELECT s.kind, s.length, r.to_sequence_id, r.relation, r.region_start, r.region_end "
        "FROM sequence s JOIN sequence_relation r ON r.from_sequence_id = s.sequence_id "
        "WHERE s.sequence_id = ?",
        [region],
    ).fetchall() == [("region", 30, ids.sequence_id(SEQ["b"]), "part_of", 11, 40)]


def test_quality_rows_from_two_tools(bundle: Path) -> None:
    con = connect(bundle)
    rows = con.execute(
        "SELECT run_id, completeness, contamination, tier, miuvig_tier, terminal_repeat, expected_length "
        "FROM quality WHERE sequence_id = ? ORDER BY 1",
        [ids.sequence_id(SEQ["a"])],
    ).fetchall()
    assert rows == [
        ("checkv", 1.0, 0.0, "Complete", "high_quality", "direct", 61.0),
        ("genomad", None, None, None, None, "direct", None),
    ]


def test_clusters_coverage_and_alignments(bundle: Path) -> None:
    con = connect(bundle)
    sizes = con.execute(
        "SELECT c.run_id, c.level, count(*) FROM cluster c JOIN cluster_member m USING (cluster_id) "
        "GROUP BY c.cluster_id, c.run_id, c.level ORDER BY 1, 3"
    ).fetchall()
    # k141_1 and NODE_1 are the same bases, so they are one member.
    assert sizes == [("loose", "genus", 2), ("votu", "species", 1), ("votu", "species", 1)]
    assert con.execute("SELECT mapped_count, mean_depth, count_unit FROM coverage").fetchall() == [
        (40, 12.5, "reads")
    ]
    assert con.execute(
        "SELECT query_start, query_end, query_strand, target_start, target_end, identity, is_primary "
        "FROM alignment"
    ).fetchall() == [(1, 40, "minus", 11, 50, 0.95, True)]


def test_a_run_can_be_added_once(bundle: Path, tools: dict[str, Path]) -> None:
    builder = Builder(bundle)
    with pytest.raises(ValueError, match="already in the bundle"):
        add_contigs(builder, tools["first"], run_id="asm1", tool="megahit")


def test_results_for_sequences_not_in_the_bundle_are_refused(
    tmp_path: Path, tools: dict[str, Path]
) -> None:
    out = tmp_path / "other" / "bundle"
    init_bundle(out, tools["samples"])
    builder = Builder(out)
    add_contigs(builder, tools["second"], run_id="asm2", tool="metaspades")
    with pytest.raises(ValueError, match="not in the bundle"):
        add_genomad(builder, tools["genomad"], tools["first"], run_id="genomad")


def test_sam_positions_follow_the_read_as_sequenced() -> None:
    forward = "r1\t99\tctg\t101\t60\t5S90M5I10M2D10M\t=\t300\t0\t*\t*\tNM:i:9"
    reverse = "r1\t147\tctg\t300\t60\t100M20S\t=\t101\t0\t*\t*\tNM:i:0"
    unmapped = "r2\t77\t*\t0\t0\t*\t*\t0\t0\t*\t*"
    first, second = sam_records(["@HD\tVN:1.6", forward, reverse, unmapped])
    assert first[:8] == ("r1", "first", "ctg", 6, 120, "plus", 101, 212)
    assert first[9] == 117 and first[10] == 2 and first[8] == pytest.approx(108 / 117)
    assert second[:8] == ("r1", "second", "ctg", 21, 120, "minus", 300, 399)


def test_paf_coordinates_become_one_based() -> None:
    (row,) = paf_records(["q\t100\t0\t50\t+\tt\t200\t20\t70\t49\t50\t60\ttp:A:S\n"])
    assert row[3:8] == (1, 50, "plus", 21, 70)
    assert row[12] is False and row[13] == 0.5 and row[14] == 0.25


def test_sequences_with_no_terminal_repeats_say_so(bundle: Path) -> None:
    con = connect(bundle)
    assert con.execute(
        "SELECT run_id, terminal_repeat FROM quality WHERE sequence_id = ?",
        [ids.sequence_id(SEQ["c"])],
    ).fetchall() == [("genomad", "none")]


def test_virsorter2_calls_scores_and_partial_sequences(
    tmp_path: Path, bundle: Path, tools: dict[str, Path]
) -> None:
    """Headers copied from VirSorter2 2.2.4 output. The values are made up."""
    from capsomer.importers.virsorter2 import add_virsorter2

    out = tmp_path / "vs2"
    tsv(
        out / "final-viral-score.tsv",
        "seqname dsDNAphage ssDNA max_score max_score_group length hallmark viral cellular",
        "k141_1||full 0.98 0.40 0.98 dsDNAphage 58 3 100.000 0.000",
        "k141_2||0_partial 0.95 0.10 0.95 dsDNAphage 30 1 80.000 20.000",
        "k141_3||lt2gene nan nan nan ssDNA 50 1 nan nan",
    )
    tsv(out / "iter-0/all-fullseq-proba.tsv", "seqname dsDNAphage ssDNA", "k141_2 0.30 0.20")
    fasta(out / "final-viral-combined.fa", {"k141_1||full": SEQ["a"][1:59], "k141_2||0_partial": SEQ["b"][20:50]})  # fmt: skip
    (out / "config.yaml").write_text("PROBA_CUTOFF: 0.5\nGROUPS: [dsDNAphage, ssDNA]\nSEQFILE: /somewhere/private\n")  # fmt: skip
    with Builder(bundle) as b:
        add_virsorter2(b, out, tools["first"], run_id="vs2", tool_version="2.2.4")
    assert conforms(validate(bundle))
    con = connect(bundle)
    rows = {
        row[0]: row[1:]
        for row in con.execute(
            "SELECT s.sequence_id, v.call, v.virus_score, v.n_hallmark_genes, v.fewer_than_two_genes "
            "FROM viral_call v JOIN sequence s USING (sequence_id) WHERE v.run_id = 'vs2'"
        ).fetchall()
    }
    region = ids.sequence_id(SEQ["b"][20:50])
    assert rows == {
        ids.sequence_id(SEQ["a"]): ("virus", 0.98, 3, False),
        ids.sequence_id(SEQ["b"]): ("none", 0.3, None, None),
        ids.sequence_id(SEQ["c"]): ("virus", None, 1, True),
        region: ("virus", 0.95, 1, False),
    }
    assert con.execute(
        "SELECT to_sequence_id, region_start, region_end FROM sequence_relation WHERE from_sequence_id = ?",
        [region],
    ).fetchall() == [(ids.sequence_id(SEQ["b"]), 21, 50)]
    settings = dict(
        con.execute("SELECT name, value FROM run_parameter WHERE run_id = 'vs2'").fetchall()
    )
    assert settings == {"proba_cutoff": "0.5", "groups": "dsDNAphage,ssDNA"}
    # Two classifiers now describe the same contig.
    assert con.execute(
        "SELECT run_id, call FROM viral_call WHERE sequence_id = ? ORDER BY 1",
        [ids.sequence_id(SEQ["b"])],
    ).fetchall() == [("genomad", "none"), ("vs2", "none")]


def test_checkm2_quality_describes_genomes(
    tmp_path: Path, bundle: Path, tools: dict[str, Path]
) -> None:
    """Header copied from CheckM2 1.1.0 output. The values are made up."""
    from capsomer.importers.checkm2 import add_checkm2

    genomes = tmp_path / "bins"
    genomes.mkdir()
    fasta(genomes / "bin1.fna", {"x": SEQ["a"], "y": SEQ["d"]})
    fasta(genomes / "bin2.fna", {"z": SEQ["b"]})
    report = tsv(
        tmp_path / "quality_report.tsv",
        "Name Completeness Contamination Completeness_Model_Used Translation_Table_Used Coding_Density "
        "Contig_N50 Average_Gene_Length Genome_Size GC_Content Total_Coding_Sequences Total_Contigs "
        "Max_Contig_Length Additional_Notes",
        "bin1 91.5 135.58 Neural~Network~(Specific~Model) 11 0.9 60 300.0 100 0.5 2 2 60 None",
        "bin2 40.0 0.5 Gradient~Boost~(General~Model) 11 0.8 80 300.0 80 0.5 1 1 80 None",
    )
    with Builder(bundle) as b:
        add_checkm2(b, report, genomes, run_id="checkm2", binning_run="asm1")
    assert conforms(validate(bundle))
    con = connect(bundle)
    two = ids.genome_id([ids.sequence_id(SEQ["a"]), ids.sequence_id(SEQ["d"])])
    assert con.execute(
        "SELECT g.genome_id, g.kind, g.n_sequences, g.run_id, q.completeness, q.contamination, q.completeness_method "
        "FROM quality q JOIN genome g USING (genome_id) WHERE q.run_id = 'checkm2' ORDER BY g.n_sequences DESC"
    ).fetchall() == [
        (two, "bin", 2, "asm1", 0.915, 1.3558, "Neural Network (Specific Model)"),
        (ids.genome_id([ids.sequence_id(SEQ["b"])]), "single_sequence", 1, None, 0.4, 0.005, "Gradient Boost (General Model)"),
    ]  # fmt: skip
    assert con.execute("SELECT count(*) FROM genome_member WHERE genome_id = ?", [two]).fetchall() == [(2,)]  # fmt: skip
