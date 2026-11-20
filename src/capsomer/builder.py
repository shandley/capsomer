"""Add rows to a bundle, new or existing, one tool at a time.

An importer opens a Builder on a bundle folder, adds the rows its tool
produced, and writes the bundle back. Rows already in the bundle are kept.
"""

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import TracebackType
from typing import Any, Self, cast

import duckdb

from capsomer.bundle import MANIFEST_NAME, build_manifest, literal, quote, relation_columns
from capsomer.schema import load_catalog


class Builder:
    def __init__(self, bundle_dir: Path, *, bundle_id: str | None = None, title: str | None = None):
        self.dir = bundle_dir
        self.dir.mkdir(parents=True, exist_ok=True)
        self.work = bundle_dir.parent / "work"
        self.work.mkdir(parents=True, exist_ok=True)
        self.catalog = load_catalog()
        self.con = duckdb.connect()
        self.tables: set[str] = set()
        path = self.dir / MANIFEST_NAME
        old = cast(dict[str, Any], json.loads(path.read_text())) if path.exists() else {}
        self.bundle_id = bundle_id or str(old.get("bundle_id") or bundle_dir.resolve().parent.name)
        self.title = title or str(old.get("title") or self.bundle_id)
        self.notes: list[str] = [str(n) for n in old.get("notes", [])]
        self.extra: dict[str, dict[str, str]] = {}
        self.access: dict[str, str] = {}
        self.custom: dict[str, dict[str, Any]] = {}
        for name, body in cast(dict[str, dict[str, Any]], old.get("tables", {})).items():
            columns = cast(dict[str, dict[str, Any]], body.get("columns", {}))
            self.access[name] = str(body.get("access", "open"))
            self.extra[name] = {
                c: str(e["description"]) for c, e in columns.items() if e.get("description")
            }
            if name not in self.catalog.tables:
                self.custom[name] = {
                    "shape": body.get("shape"),
                    "one_row_is": body.get("one_row_is"),
                    "description": body.get("description", ""),
                    "key": body.get("key", []),
                }
        for parquet in sorted(self.dir.glob("*.parquet")):
            self.con.execute(
                f"CREATE TABLE {quote(parquet.stem)} AS SELECT * FROM read_parquet({literal(str(parquet))})"
            )
            self.tables.add(parquet.stem)

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        kind: type[BaseException] | None,
        value: BaseException | None,
        trace: TracebackType | None,
    ) -> None:
        if kind is None:
            self.write()
        self.con.close()

    def scalar(self, sql: str, parameters: Sequence[Any] = ()) -> Any:
        row = self.con.execute(sql, list(parameters)).fetchone()
        return row[0] if row else None

    def add(self, name: str, query: str) -> int:
        """Append the rows of a query to a table and return how many were added.

        Catalog columns are cast to their catalog types, and columns the query
        lacks are null. Other columns are kept as extra columns. In a thing
        table, a row whose identifier the bundle already holds is skipped.
        """
        table = self.catalog.tables.get(name)
        self.con.execute(f"CREATE OR REPLACE TEMP VIEW _new AS {query}")
        present = relation_columns(self.con, "_new")
        select = "*"
        if table:
            parts = [
                f"CAST({quote(c.name) if c.name in present else 'NULL'} AS {c.sql_type}) AS {quote(c.name)}"
                for c in table.columns
            ]
            parts += [quote(extra) for extra in present if extra not in table.column_names]
            select = ", ".join(parts)
        source = f"SELECT {select} FROM _new"
        if table and table.identifier:
            ident = quote(table.identifier)
            source += f" QUALIFY row_number() OVER (PARTITION BY {ident}) = 1"
            if name in self.tables:
                source = f"SELECT * FROM ({source}) WHERE {ident} NOT IN (SELECT {ident} FROM {quote(name)})"
        self.con.execute(f"CREATE OR REPLACE TEMP TABLE _rows AS {source}")
        added = int(self.scalar("SELECT count(*) FROM _rows"))
        if name in self.tables:
            self.con.execute(
                f"CREATE OR REPLACE TABLE {quote(name)} AS "
                f"SELECT * FROM {quote(name)} UNION ALL BY NAME SELECT * FROM _rows"
            )
        else:
            self.con.execute(f"CREATE TABLE {quote(name)} AS SELECT * FROM _rows")
            self.tables.add(name)
        return added

    def add_rows(
        self, name: str, columns: Sequence[str], rows: Sequence[Sequence[Any]], types: str = ""
    ) -> int:
        """Append rows given as Python values. types is a column list for CREATE TABLE."""
        spec = types or ", ".join(f"{quote(c)} VARCHAR" for c in columns)
        self.con.execute(f"CREATE OR REPLACE TEMP TABLE _given ({spec})")
        if rows:
            marks = ", ".join("?" for _ in columns)
            self.con.executemany(f"INSERT INTO _given VALUES ({marks})", [list(r) for r in rows])
        return self.add(name, "SELECT * FROM _given")

    def add_run(
        self,
        run_id: str,
        *,
        tool: str,
        step: str,
        tool_version: str | None = None,
        description: str | None = None,
        parent: str | None = None,
        parameters: Mapping[str, object] | None = None,
        libraries: Sequence[str] = (),
        inputs: Sequence[str] = (),
        references: Sequence[str] = (),
    ) -> None:
        """Record a run with its settings and what it read. A run can be added once."""
        if "run" in self.tables and self.scalar(
            "SELECT count(*) FROM run WHERE run_id = ?", [run_id]
        ):
            raise ValueError(f"run {run_id} is already in the bundle")
        self.add_rows(
            "run",
            ("run_id", "tool", "tool_version", "step", "parent_run_id", "description"),
            [(run_id, tool, tool_version, step, parent, description)],
        )
        settings = [(run_id, k, str(v)) for k, v in (parameters or {}).items() if v is not None]
        if settings:
            self.add_rows("run_parameter", ("run_id", "name", "value"), settings)
        for table, column, found in (
            ("run_library", "library_id", libraries),
            ("run_input", "upstream_run_id", inputs),
            ("run_reference", "reference_release_id", references),
        ):
            if found:
                self.add_rows(table, ("run_id", column), [(run_id, value) for value in found])

    def describe(self, table: str, columns: Mapping[str, str]) -> None:
        """Give descriptions for extra columns of a table."""
        self.extra.setdefault(table, {}).update(columns)

    def declare(self, table: str, *, shape: str, one_row_is: str, description: str, key: list[str]):
        """Declare a table outside the catalog."""
        self.custom[table] = {
            "shape": shape,
            "one_row_is": one_row_is,
            "description": description,
            "key": key,
        }

    def note(self, text: str) -> None:
        if text not in self.notes:
            self.notes.append(text)

    def write(self) -> dict[str, Any]:
        """Write every table as Parquet, and the manifest."""
        for name in sorted(self.tables):
            table = self.catalog.tables.get(name)
            order = ", ".join(quote(k) for k in table.key) if table and table.key else "ALL"
            path = self.dir / f"{name}.parquet"
            self.con.execute(
                f"COPY (SELECT * FROM {quote(name)} ORDER BY {order}) TO {literal(str(path))} "
                "(FORMAT parquet, COMPRESSION zstd)"
            )
        return build_manifest(
            self.dir,
            bundle_id=self.bundle_id,
            title=self.title,
            extra_columns=self.extra,
            access=self.access,
            notes=self.notes,
            custom_tables=self.custom,
        )
