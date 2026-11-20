"""Import read mappings, coverage tables and sequence-to-sequence alignments.

Three inputs are read: a CoverM contig table (coverage per library), a SAM or
BAM file of reads mapped to sequences (one alignment row per mapped read, when
the bundle holds the reads), and a PAF file of sequences aligned to sequences.
"""

import re
import subprocess
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path

from capsomer.builder import Builder
from capsomer.bundle import literal, quote
from capsomer.importers.common import load_fasta

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COVERM = {
    " Read Count": ("mapped_count", "BIGINT"),
    " Mean": ("mean_depth", "DOUBLE"),
    " Covered Fraction": ("covered_fraction", "DOUBLE"),
    " Covered Bases": ("covered_bases", "BIGINT"),
}
ALIGNMENT_TYPES = (
    "query VARCHAR, query_end_of_pair VARCHAR, target VARCHAR, query_start BIGINT, query_end BIGINT, "
    "query_strand VARCHAR, target_start BIGINT, target_end BIGINT, identity DOUBLE, "
    "alignment_length BIGINT, mismatches BIGINT, mapping_quality BIGINT, is_primary BOOLEAN, "
    "query_coverage DOUBLE, target_coverage DOUBLE"
)


def _require_sequences(builder: Builder, table: str, fasta: Path) -> None:
    missing = builder.scalar(
        f"SELECT count(*) FROM {table} WHERE sequence_id NOT IN (SELECT sequence_id FROM sequence)"
    )
    if missing:
        raise ValueError(f"{missing} sequences of {fasta} are not in the bundle; add them first")


def add_coverm(
    builder: Builder,
    table: Path,
    fasta: Path,
    *,
    run_id: str,
    libraries: Mapping[str, str],
    tool_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    parameters: Mapping[str, object] | None = None,
) -> int:
    """Add coverage from a CoverM contig table and return the rows added.

    libraries maps each sample name in the table's header to a library in the
    bundle. Sequences with no mapped reads in a library get no row.
    """
    con = builder.con
    load_fasta(con, fasta)
    _require_sequences(builder, "_fasta", fasta)
    header = table.read_text().split("\n", 1)[0].split("\t")
    con.execute(
        f"CREATE OR REPLACE TEMP TABLE _coverm AS SELECT * FROM read_csv({literal(str(table))}, "
        "delim='\\t', header=true, all_varchar=true)"
    )
    builder.add_run(
        run_id,
        tool="CoverM",
        tool_version=tool_version,
        step="read_mapping",
        description="Coverage of sequences by mapped reads.",
        parent=parent,
        parameters=parameters,
        libraries=sorted(set(libraries.values())),
        inputs=inputs,
    )
    added = 0
    for sample, library in libraries.items():
        found = {
            column: target
            for suffix, target in COVERM.items()
            if (column := sample + suffix) in header
        }
        if not found:
            raise ValueError(f"{table} has no columns for sample {sample}")
        select = ", ".join(
            f"CAST(c.{quote(column)} AS {sql_type}) AS {name}"
            for column, (name, sql_type) in found.items()
        )
        observed = (
            f"WHERE CAST(c.{quote(sample + ' Read Count')} AS BIGINT) > 0"
            if sample + " Read Count" in header
            else ""
        )
        added += builder.add(
            "coverage",
            f"SELECT f.sequence_id, {literal(library)} AS library_id, {literal(run_id)} AS run_id, "
            f"{select}, 'reads' AS count_unit FROM _coverm c JOIN _fasta f ON f.name = c.{quote(header[0])} "
            f"{observed} QUALIFY row_number() OVER (PARTITION BY f.sequence_id ORDER BY f.name) = 1",
        )
    return added


def sam_records(lines: Iterable[str]) -> Iterator[tuple[object, ...]]:
    """One tuple per mapped record of SAM text, in the column order of ALIGNMENT_TYPES.

    Query positions count from the first base of the read as sequenced, so a
    read that mapped to the minus strand has its soft clips swapped.
    """
    for line in lines:
        if line.startswith("@"):
            continue
        field = line.rstrip("\n").split("\t")
        flag = int(field[1])
        if flag & 4 or field[2] == "*":
            continue
        ops = [(int(n), op) for n, op in CIGAR.findall(field[5])]
        on_target = sum(n for n, op in ops if op in "MDN=X")
        on_query = sum(n for n, op in ops if op in "MI=X")
        columns = sum(n for n, op in ops if op in "MID=X")
        gaps = sum(n for n, op in ops if op in "ID")
        lead = ops[0][0] if ops and ops[0][1] in "SH" else 0
        tail = ops[-1][0] if len(ops) > 1 and ops[-1][1] in "SH" else 0
        reverse = bool(flag & 16)
        start = (tail if reverse else lead) + 1
        edits = next((int(tag[5:]) for tag in field[11:] if tag.startswith("NM:i:")), None)
        position = int(field[3])
        yield (
            field[0],
            "first" if flag & 64 else "second" if flag & 128 else "",
            field[2],
            start,
            start + on_query - 1,
            "minus" if reverse else "plus",
            position,
            position + on_target - 1,
            None if edits is None or not columns else (columns - edits) / columns,
            columns,
            None if edits is None else edits - gaps,
            int(field[4]),
            not flag & 0x900,
            on_query / (lead + on_query + tail) if on_query else None,
            None,
        )


def paf_records(lines: Iterable[str]) -> Iterator[tuple[object, ...]]:
    """One tuple per PAF line, in the column order of ALIGNMENT_TYPES."""
    for line in lines:
        field = line.rstrip("\n").split("\t")
        if len(field) < 12:
            continue
        q_len, q_start, q_end = int(field[1]), int(field[2]), int(field[3])
        t_len, t_start, t_end = int(field[6]), int(field[7]), int(field[8])
        matches, columns = int(field[9]), int(field[10])
        kind = next((tag[5:] for tag in field[12:] if tag.startswith("tp:A:")), None)
        yield (
            field[0],
            "",
            field[5],
            q_start + 1,
            q_end,
            "minus" if field[4] == "-" else "plus",
            t_start + 1,
            t_end,
            matches / columns if columns else None,
            columns,
            None,
            int(field[11]),
            None if kind is None else kind == "P",
            (q_end - q_start) / q_len if q_len else None,
            (t_end - t_start) / t_len if t_len else None,
        )


def _load(builder: Builder, records: Iterator[tuple[object, ...]]) -> int:
    """Fill the temp table _aln from parsed records, by way of a file in the work folder."""
    path = builder.work / "alignments.tsv"
    count = 0
    with path.open("w") as out:
        for record in records:
            out.write("\t".join("" if v is None else str(v) for v in record) + "\n")
            count += 1
    builder.con.execute(f"CREATE OR REPLACE TEMP TABLE _aln ({ALIGNMENT_TYPES})")
    if count:
        types = {part.split()[0]: part.split()[1] for part in ALIGNMENT_TYPES.split(", ")}
        builder.con.execute(
            f"INSERT INTO _aln SELECT * FROM read_csv({literal(str(path))}, delim='\\t', header=false, "
            f"quote='', nullstr='', columns={types!r})"
        )
    path.unlink()
    return count


def _sam_lines(path: Path) -> Iterator[str]:
    if path.suffix == ".sam":
        with path.open() as handle:
            yield from handle
        return
    try:
        process = subprocess.Popen(
            ["samtools", "view", str(path)], stdout=subprocess.PIPE, text=True
        )
    except FileNotFoundError as problem:
        raise RuntimeError("reading a BAM file needs samtools on the PATH") from problem
    assert process.stdout is not None
    yield from process.stdout
    if process.wait():
        raise RuntimeError(f"samtools view failed on {path}")


def add_read_mapping(
    builder: Builder,
    fasta: Path,
    *,
    run_id: str,
    tool: str,
    library_id: str,
    alignments: Path | None = None,
    tool_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    parameters: Mapping[str, object] | None = None,
) -> dict[str, int]:
    """Record a run that mapped one library's reads to the sequences of a FASTA file.

    With a SAM or BAM file, and when the bundle holds the library's reads,
    every mapped read becomes an alignment row. A read is matched by its name
    and its end of the pair. Returns the alignments added and the mapped
    records whose read the bundle does not hold.
    """
    con = builder.con
    load_fasta(con, fasta)
    _require_sequences(builder, "_fasta", fasta)
    builder.add_run(
        run_id,
        tool=tool,
        tool_version=tool_version,
        step="read_mapping",
        description="Mapping of reads to sequences.",
        parent=parent,
        parameters=parameters,
        libraries=[library_id],
        inputs=inputs,
    )
    if alignments is None or "read" not in builder.tables:
        return {"alignment": 0, "unmatched": 0}
    records = _load(builder, sam_records(_sam_lines(alignments)))
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _reads AS SELECT read_id, name, "
        "regexp_replace(name, '/[12]$', '') AS fragment, read_end FROM read WHERE library_id = ?",
        [library_id],
    )
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _target AS SELECT DISTINCT name, sequence_id FROM _fasta"
    )
    matched = (
        "SELECT r.read_id, a.* FROM _aln a JOIN _reads r "
        "ON r.fragment = a.query AND r.read_end = a.query_end_of_pair WHERE a.query_end_of_pair IS NOT NULL "
        "UNION ALL SELECT r.read_id, a.* FROM _aln a JOIN _reads r ON r.name = a.query "
        "WHERE a.query_end_of_pair IS NULL"
    )
    added = builder.add(
        "alignment",
        f"SELECT m.read_id AS query_read_id, t.sequence_id AS target_sequence_id, {literal(run_id)} AS run_id, "
        "'nucleotide_to_nucleotide' AS search_type, m.query_start, m.query_end, m.query_strand, "
        "m.target_start, m.target_end, m.identity, m.alignment_length, m.mismatches, m.query_coverage, "
        f"m.mapping_quality, m.is_primary FROM ({matched}) m JOIN _target t ON t.name = m.target",
    )
    return {"alignment": added, "unmatched": records - added}


def add_paf(
    builder: Builder,
    paf: Path,
    query_fasta: Path,
    target_fasta: Path,
    *,
    run_id: str,
    tool: str,
    tool_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    parameters: Mapping[str, object] | None = None,
) -> int:
    """Add alignments of sequences to sequences from a PAF file and return the rows added."""
    con = builder.con
    load_fasta(con, query_fasta, "_query")
    load_fasta(con, target_fasta, "_fasta")
    _require_sequences(builder, "_query", query_fasta)
    _require_sequences(builder, "_fasta", target_fasta)
    with paf.open() as handle:
        _load(builder, paf_records(handle))
    builder.add_run(
        run_id,
        tool=tool,
        tool_version=tool_version,
        step="search",
        description="Alignment of sequences to sequences.",
        parent=parent,
        parameters=parameters,
        inputs=inputs,
    )
    return builder.add(
        "alignment",
        f"SELECT q.sequence_id AS query_sequence_id, t.sequence_id AS target_sequence_id, {literal(run_id)} AS run_id, "
        "'nucleotide_to_nucleotide' AS search_type, a.query_start, a.query_end, a.query_strand, "
        "a.target_start, a.target_end, a.identity, a.alignment_length, a.query_coverage, a.target_coverage, "
        "a.mapping_quality, a.is_primary FROM _aln a "
        "JOIN (SELECT DISTINCT name, sequence_id FROM _query) q ON q.name = a.query "
        "JOIN (SELECT DISTINCT name, sequence_id FROM _fasta) t ON t.name = a.target",
    )
