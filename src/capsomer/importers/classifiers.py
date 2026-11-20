"""Import viral calls from Cenote-Taker 3, PhaBOX2 (PhaMer) and Jaeger.

Each tool reports by sequence name. The names are joined to the bundle through
the FASTA file the tool was given. A sequence the tool was given and did not
call viral gets a row with call none, so "examined and not called" is recorded.
"""

from collections.abc import Mapping, Sequence
from pathlib import Path

from capsomer.builder import Builder
from capsomer.bundle import literal
from capsomer.importers.common import load_fasta, tsv

JAEGER_CALL = {
    "phage": "virus",
    "virus": "virus",
    "plasmid": "plasmid",
    "bacteria": "cellular",
    "archaea": "cellular",
    "eukarya": "cellular",
}
JAEGER_CLASSES = ("bacteria", "phage", "eukarya", "archaea", "plasmid", "virus")


def _start(
    builder: Builder,
    fasta: Path,
    *,
    run_id: str,
    tool: str,
    tool_version: str | None,
    description: str,
    inputs: Sequence[str],
    parent: str | None,
    parameters: Mapping[str, object] | None,
) -> None:
    load_fasta(builder.con, fasta)
    missing = builder.scalar(
        "SELECT count(*) FROM _fasta WHERE sequence_id NOT IN (SELECT sequence_id FROM sequence)"
    )
    if missing:
        raise ValueError(f"{missing} sequences of {fasta} are not in the bundle; add them first")
    builder.add_run(
        run_id,
        tool=tool,
        tool_version=tool_version,
        step="viral_identification",
        description=description,
        parent=parent,
        parameters=parameters,
        inputs=inputs,
    )


def _one_per_sequence(select: str, order: str = "f.name") -> str:
    return f"{select} QUALIFY row_number() OVER (PARTITION BY f.sequence_id ORDER BY {order}) = 1"


def add_cenotetaker3(
    builder: Builder,
    results: Path,
    fasta: Path,
    *,
    run_id: str,
    tool_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    parameters: Mapping[str, object] | None = None,
) -> dict[str, int]:
    """Add one Cenote-Taker 3 run made without prophage pruning.

    Cenote-Taker 3 lists the sequences it calls viral and gives no score. Its
    terminal-repeat finding goes to quality.
    """
    summary = sorted(results.glob("*_virus_summary.tsv"))
    if not summary:
        raise FileNotFoundError(f"no *_virus_summary.tsv in {results}")
    con = builder.con
    con.execute(f"CREATE OR REPLACE TEMP TABLE _ct3 AS SELECT * FROM {tsv(summary[0])}")
    pruned = sorted(results.glob("*_prune_summary.tsv"))
    if pruned:
        chunks = builder.scalar(
            f"SELECT count(*) FROM {tsv(pruned[0])} WHERE chunk_name IS NOT NULL AND chunk_name <> 'NaN'"
        )
        if chunks:
            raise NotImplementedError("Cenote-Taker 3 runs with prophage pruning are not read yet")
    _start(
        builder,
        fasta,
        run_id=run_id,
        tool="Cenote-Taker 3",
        tool_version=tool_version,
        description="Identification of viral sequences by hallmark genes.",
        inputs=inputs,
        parent=parent,
        parameters=parameters,
    )
    run = literal(run_id)
    hallmarks = (
        "CAST(c.virion_hallmark_count AS BIGINT) + CAST(c.rep_hallmark_count AS BIGINT) "
        "+ CAST(c.RDRP_hallmark_count AS BIGINT)"
    )
    added = {
        "viral_call": builder.add(
            "viral_call",
            _one_per_sequence(
                f"SELECT f.sequence_id, {run} AS run_id, "
                "CASE WHEN c.input_name IS NULL THEN 'none' ELSE 'virus' END AS call, "
                f"{hallmarks} AS n_hallmark_genes, c.organism AS organism_guess "
                "FROM _fasta f LEFT JOIN _ct3 c ON c.input_name = f.name"
            ),
        ),
        "quality": builder.add(
            "quality",
            _one_per_sequence(
                f"SELECT f.sequence_id, {run} AS run_id, CASE c.end_feature WHEN 'DTR' THEN 'direct' "
                "WHEN 'ITR' THEN 'inverted' WHEN 'None' THEN 'none' END AS terminal_repeat "
                "FROM _fasta f JOIN _ct3 c ON c.input_name = f.name "
                "WHERE c.end_feature IN ('DTR', 'ITR', 'None')"
            ),
        ),
    }
    builder.describe(
        "viral_call",
        {"organism_guess": "The organism name Cenote-Taker 3 proposes for the sequence."},
    )
    builder.note(
        "Cenote-Taker 3 gives no score. It calls a sequence viral when it carries enough hallmark "
        "genes; the minimum is a setting of the run."
    )
    return added


def add_phabox2(
    builder: Builder,
    results: Path,
    fasta: Path,
    *,
    run_id: str,
    tool_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    parameters: Mapping[str, object] | None = None,
) -> int:
    """Add the virus identification step (PhaMer) of one PhaBOX2 run."""
    table = results / "final_prediction" / "phamer_prediction.tsv"
    if not table.exists():
        table = results / "phamer_prediction.tsv"
    con = builder.con
    con.execute(f"CREATE OR REPLACE TEMP TABLE _phamer AS SELECT * FROM {tsv(table)}")
    _start(
        builder,
        fasta,
        run_id=run_id,
        tool="PhaBOX2 PhaMer",
        tool_version=tool_version,
        description="Identification of viral sequences.",
        inputs=inputs,
        parent=parent,
        parameters=parameters,
    )
    added = builder.add(
        "viral_call",
        _one_per_sequence(
            f"SELECT f.sequence_id, {literal(run_id)} AS run_id, "
            "CASE WHEN p.Pred = 'virus' THEN 'virus' ELSE 'none' END AS call, "
            "TRY_CAST(p.PhaMerScore AS DOUBLE) AS virus_score, p.PhaMerConfidence AS tool_confidence, "
            "TRY_CAST(p.Proportion AS DOUBLE) AS matched_protein_fraction, p.Pred AS tool_label "
            "FROM _fasta f LEFT JOIN _phamer p ON p.Accession = f.name"
        ),
    )
    builder.describe(
        "viral_call",
        {
            "tool_confidence": "The confidence statement the tool attaches to its call, in its own words.",
            "matched_protein_fraction": "PhaBOX2: fraction of the sequence's proteins that matched its reference.",
            "tool_label": "The class the tool itself reports, before it is mapped to call.",
        },
    )
    return added


def add_jaeger(
    builder: Builder,
    table: Path,
    fasta: Path,
    *,
    run_id: str,
    tool_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    parameters: Mapping[str, object] | None = None,
) -> int:
    """Add one Jaeger run from its per-sequence table.

    Jaeger sorts sequences into six classes. Phage and virus become call
    virus, plasmid stays plasmid, and bacteria, archaea and eukarya become
    cellular. The six scores are kept in detail.
    """
    con = builder.con
    con.execute(f"CREATE OR REPLACE TEMP TABLE _jaeger AS SELECT * FROM {tsv(table)}")
    present = {str(r[0]) for r in con.execute("DESCRIBE _jaeger").fetchall()}
    classes = [c for c in JAEGER_CLASSES if f"{c}_score" in present]
    _start(
        builder,
        fasta,
        run_id=run_id,
        tool="Jaeger",
        tool_version=tool_version,
        description="Classification of sequences by a neural network, without reference searches.",
        inputs=inputs,
        parent=parent,
        parameters=parameters,
    )

    def best(wanted: Sequence[str]) -> str:
        found = [f"TRY_CAST(j.{c}_score AS DOUBLE)" for c in wanted if c in classes]
        return f"greatest({', '.join(found)})" if found else "NULL"

    mapping = " ".join(f"WHEN {literal(k)} THEN {literal(v)}" for k, v in JAEGER_CALL.items())
    detail = (
        "json_object("
        + ", ".join(f"{literal(c)}, TRY_CAST(j.{c}_score AS DOUBLE)" for c in classes)
        + ")"
    )
    added = builder.add(
        "viral_call",
        _one_per_sequence(
            f"SELECT f.sequence_id, {literal(run_id)} AS run_id, "
            f"CASE j.prediction {mapping} ELSE 'none' END AS call, "
            f"{best(('phage', 'virus'))} AS virus_score, {best(('plasmid',))} AS plasmid_score, "
            f"{best(('bacteria', 'archaea', 'eukarya'))} AS cellular_score, "
            f"CASE WHEN j.contig_id IS NOT NULL THEN {detail} END AS detail, j.prediction AS tool_label, "
            "TRY_CAST(j.reliability_score AS DOUBLE) AS reliability_score "
            "FROM _fasta f LEFT JOIN _jaeger j ON j.contig_id = f.name"
        ),
    )
    builder.describe(
        "viral_call",
        {
            "tool_label": "The class the tool itself reports, before it is mapped to call.",
            "reliability_score": "Jaeger's reliability score for its prediction.",
        },
    )
    return added
