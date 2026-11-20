"""Import clusters of sequences from a two-column table.

Two layouts are read. vclust writes a header and then member, cluster. MMseqs2
writes representative, member with no header. The level of the clusters
(species, genus, or any other) is whatever the researcher says it is.
"""

from collections.abc import Mapping, Sequence
from pathlib import Path

from capsomer.builder import Builder
from capsomer.bundle import literal
from capsomer.importers.common import load_fasta

LAYOUTS = {"vclust": ("member", "cluster"), "mmseqs": ("cluster", "member")}


def add_clusters(
    builder: Builder,
    table: Path,
    fasta: Path,
    *,
    run_id: str,
    tool: str,
    level: str,
    layout: str = "vclust",
    tool_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    parameters: Mapping[str, object] | None = None,
) -> int:
    """Add the clusters of one clustering run and return how many there are.

    When the table names a cluster by one of its members, that member is the
    representative. Otherwise the longest member is.
    """
    con = builder.con
    load_fasta(con, fasta)
    first, second = LAYOUTS[layout]
    con.execute(
        f"CREATE OR REPLACE TEMP TABLE _pairs AS SELECT column0 AS {first}, column1 AS {second} "
        f"FROM read_csv({literal(str(table))}, delim='\\t', header=false, all_varchar=true, "
        "columns={'column0': 'VARCHAR', 'column1': 'VARCHAR'}) "
        "WHERE NOT (column0 = 'object' AND column1 = 'cluster')"
    )
    unknown = builder.scalar(
        "SELECT count(*) FROM _pairs WHERE member NOT IN (SELECT name FROM _fasta)"
    )
    if unknown:
        raise ValueError(f"{unknown} members in {table} are not in {fasta}")
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _members AS "
        "SELECT p.cluster, f.sequence_id, f.name, f.length, p.cluster = p.member AS named "
        "FROM _pairs p JOIN _fasta f ON f.name = p.member"
    )
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _representative AS SELECT cluster, sequence_id FROM _members "
        "QUALIFY row_number() OVER (PARTITION BY cluster ORDER BY named DESC, length DESC, name) = 1"
    )
    builder.add_run(
        run_id,
        tool=tool,
        tool_version=tool_version,
        step="sequence_clustering",
        description=f"Clustering of sequences at the {level} level.",
        parent=parent,
        parameters=parameters,
        inputs=inputs,
    )
    prefix = literal(run_id + ":")
    added = builder.add(
        "cluster",
        f"SELECT {prefix} || sequence_id AS cluster_id, {literal(level)} AS level, "
        f"{literal(run_id)} AS run_id, sequence_id AS representative_sequence_id FROM _representative",
    )
    builder.add(
        "cluster_member",
        f"SELECT DISTINCT {prefix} || r.sequence_id AS cluster_id, m.sequence_id, "
        "m.sequence_id = r.sequence_id AS is_representative "
        "FROM _members m JOIN _representative r USING (cluster)",
    )
    return added
