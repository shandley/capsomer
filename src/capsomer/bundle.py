"""Write catalog tables as Parquet and build the bundle manifest."""

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import duckdb

from capsomer.schema import Catalog, Table, load_catalog

MANIFEST_NAME = "manifest.json"
NUMERIC_TYPES = ("BIGINT", "DOUBLE", "INTEGER", "FLOAT", "HUGEINT", "UBIGINT")


def quote(identifier: str) -> str:
    return '"' + identifier.replace('"', '""') + '"'


def literal(text: str) -> str:
    return "'" + text.replace("'", "''") + "'"


def relation_columns(con: duckdb.DuckDBPyConnection, relation: str) -> dict[str, str]:
    """Column name -> SQL type, for a table, view or parquet scan."""
    rows = con.execute(f"DESCRIBE SELECT * FROM {relation}").fetchall()
    return {str(row[0]): str(row[1]) for row in rows}


def write_table(
    con: duckdb.DuckDBPyConnection,
    name: str,
    source: str,
    out_dir: Path,
    catalog: Catalog | None = None,
) -> int:
    """Write one catalog table from a query.

    Catalog columns come first, in catalog order and cast to the catalog type.
    Columns the query lacks are written as nulls. Columns the catalog does not
    know are kept after them as extra columns.
    """
    table = (catalog or load_catalog()).tables[name]
    con.execute(f"CREATE OR REPLACE TEMP VIEW _capsomer_src AS {source}")
    present = relation_columns(con, "_capsomer_src")
    select = [
        f"CAST({quote(c.name) if c.name in present else 'NULL'} AS {c.sql_type}) AS {quote(c.name)}"
        for c in table.columns
    ]
    select += [quote(extra) for extra in present if extra not in table.column_names]
    order = ", ".join(quote(k) for k in table.key) if table.key else "ALL"
    path = out_dir / f"{name}.parquet"
    con.execute(
        f"COPY (SELECT {', '.join(select)} FROM _capsomer_src ORDER BY {order}) "
        f"TO {literal(str(path))} (FORMAT parquet, COMPRESSION zstd)"
    )
    count = con.execute(f"SELECT count(*) FROM read_parquet({literal(str(path))})").fetchone()
    return int(count[0]) if count else 0


def connect(bundle_dir: Path) -> duckdb.DuckDBPyConnection:
    """Open a bundle for querying: one view per table, named after the table."""
    con = duckdb.connect()
    for path in sorted(bundle_dir.glob("*.parquet")):
        scan = f"read_parquet({literal(str(path))})"
        con.execute(f"CREATE VIEW {quote(path.stem)} AS SELECT * FROM {scan}")
    return con


def run_sql(bundle_dir: Path, sql: str, max_rows: int = 40) -> str:
    """Run one or more statements against a bundle and return the results as text.

    Statements are separated by a semicolon at the end of a line. Comment lines
    that precede a statement are printed above its result.
    """
    con = connect(bundle_dir)
    blocks: list[str] = []
    for statement in (s.strip() for s in sql.split(";\n")):
        lines = statement.rstrip(";").splitlines()
        comments = [line[2:].strip() for line in lines if line.startswith("--")]
        body = "\n".join(line for line in lines if not line.startswith("--")).strip()
        if not body:
            continue
        result = con.execute(body)
        out = list(comments)
        if result.description:
            header = [str(d[0]) for d in result.description]
            rows = [
                [("" if v is None else str(v)) for v in row] for row in result.fetchmany(max_rows)
            ]
            widths = [
                max(len(h), *(len(r[i]) for r in rows)) if rows else len(h)
                for i, h in enumerate(header)
            ]
            out.append("  ".join(h.ljust(w) for h, w in zip(header, widths, strict=True)))
            out += ["  ".join(v.ljust(w) for v, w in zip(r, widths, strict=True)) for r in rows]
        blocks.append("\n".join(out))
    return "\n\n".join(blocks)


def _column_entry(
    con: duckdb.DuckDBPyConnection, scan: str, name: str, sql_type: str, table: Table | None
) -> dict[str, Any]:
    known = table.column(name) if table else None
    entry: dict[str, Any] = {"type": sql_type, "description": known.description if known else ""}
    if known:
        if known.references:
            entry["references"] = known.references
        if known.enum:
            entry["allowed_values"] = list(known.enum)
        entry.update(known.annotations)
    col = quote(name)
    if sql_type in NUMERIC_TYPES:
        row = con.execute(f"SELECT count({col}), min({col}), max({col}) FROM {scan}").fetchone()
        if row:
            entry["non_null"], entry["min"], entry["max"] = int(row[0]), row[1], row[2]
    else:
        row = con.execute(f"SELECT count({col}) FROM {scan}").fetchone()
        if row:
            entry["non_null"] = int(row[0])
    return entry


def build_manifest(
    bundle_dir: Path,
    bundle_id: str,
    title: str,
    extra_columns: dict[str, dict[str, str]] | None = None,
    access: dict[str, str] | None = None,
    notes: list[str] | None = None,
) -> dict[str, Any]:
    """Describe every Parquet table in a bundle and write manifest.json.

    extra_columns gives descriptions for columns outside the catalog, as
    {table: {column: description}}. access gives an access tier per table and
    defaults to "open".
    """
    catalog = load_catalog()
    con = duckdb.connect()
    tables: dict[str, Any] = {}
    for path in sorted(bundle_dir.glob("*.parquet")):
        name = path.stem
        scan = f"read_parquet({literal(str(path))})"
        known = catalog.tables.get(name)
        described = (extra_columns or {}).get(name, {})
        columns: dict[str, Any] = {}
        for column, sql_type in relation_columns(con, scan).items():
            entry = _column_entry(con, scan, column, sql_type, known)
            if not entry["description"] and column in described:
                entry["description"] = described[column]
            if known and known.column(column) is None:
                entry["extra"] = True
            columns[column] = entry
        rows = con.execute(f"SELECT count(*) FROM {scan}").fetchone()
        tables[name] = {
            "shape": known.shape if known else None,
            "catalog_entry": f"capsomer:{name}@{catalog.version}" if known else None,
            "one_row_is": known.row if known else None,
            "description": known.description if known else "",
            "key": list(known.key) if known else [],
            "row_count": int(rows[0]) if rows else 0,
            "access": (access or {}).get(name, "open"),
            "columns": columns,
        }
    manifest: dict[str, Any] = {
        "standard": "capsomer",
        "standard_version": catalog.version,
        "bundle_id": bundle_id,
        "title": title,
        "created": datetime.now(UTC).isoformat(timespec="seconds"),
        "conventions": {
            "joins": "A column named <table>_id, or ending in it, refers to that table.",
            "coordinates": "1-based, inclusive.",
            "missing_count_row": "A thing with no row for a library was not observed there.",
            "evalue": "Comparable within one run only. Bit scores compare across runs.",
        },
        "notes": notes or [],
        "tables": tables,
    }
    (bundle_dir / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    return manifest
