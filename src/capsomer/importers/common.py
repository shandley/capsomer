"""Helpers shared by importers."""

import re
from collections.abc import Iterator
from pathlib import Path

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
