"""Helpers shared by importers."""

import re
from collections.abc import Iterator
from pathlib import Path

import duckdb

from capsomer import ids
from capsomer.bundle import literal


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def tsv(path: Path) -> str:
    return (
        f"read_csv({literal(str(path))}, delim='\\t', header=true, all_varchar=true, "
        "quote='', escape='', nullstr='NA', null_padding=true, strict_mode=false)"
    )


def values(rows: list[tuple[str | None, ...]], columns: tuple[str, ...]) -> str:
    body = ", ".join(
        "(" + ", ".join("NULL" if v is None else literal(v) for v in row) + ")" for row in rows
    )
    return f"SELECT * FROM (VALUES {body}) AS t({', '.join(columns)})"


def snake(name: str) -> str:
    return re.sub(r"[^0-9a-z]+", "_", name.lower()).strip("_")


def load_fasta(con: duckdb.DuckDBPyConnection, fasta: Path, table: str = "_fasta") -> int:
    """Load a FASTA file as a temp table and return the number of records.

    Columns: name (the header up to the first space), header, sequence_id,
    length, gc_content, bases. Tools report results by name, and the
    identifier computed here from the bases is what joins them to the bundle.
    """
    con.execute(
        f"CREATE OR REPLACE TEMP TABLE {table} (name VARCHAR, header VARCHAR, "
        "sequence_id VARCHAR, length BIGINT, gc_content DOUBLE, bases VARCHAR)"
    )
    rows = [
        (
            header.split()[0],
            header,
            ids.sequence_id(bases.upper()),
            len(bases),
            ids.gc_content(bases.upper()),
            bases.upper(),
        )
        for header, bases in fasta_records(fasta)
    ]
    if rows:
        con.executemany(f"INSERT INTO {table} VALUES (?, ?, ?, ?, ?, ?)", rows)
    return len(rows)
