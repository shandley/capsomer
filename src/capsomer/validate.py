"""Check a bundle against the catalog. Each check has a severity.

critical: the bundle cannot be read reliably.
error:    the bundle can be read, but an analysis of it would be wrong.
warning:  the bundle is usable but incomplete.

The checks concern structure only. None judges an analytical choice.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import duckdb

from capsomer import ids
from capsomer.bundle import MANIFEST_NAME, literal, quote, relation_columns
from capsomer.schema import Table, load_catalog

SEVERITIES = ("critical", "error", "warning")
COUNT_TABLES = ("read_count", "coverage", "taxon_profile")
DIGEST_SAMPLE = 2000


@dataclass(frozen=True)
class Finding:
    severity: str
    check: str
    table: str
    message: str
    failing_rows: int = 0


class _Checker:
    def __init__(self, bundle_dir: Path) -> None:
        self.dir = bundle_dir
        self.catalog = load_catalog()
        self.con = duckdb.connect()
        self.findings: list[Finding] = []
        self.present: dict[str, dict[str, str]] = {}
        for path in sorted(bundle_dir.glob("*.parquet")):
            scan = f"read_parquet({literal(str(path))})"
            self.con.execute(f"CREATE VIEW {quote(path.stem)} AS SELECT * FROM {scan}")
            self.present[path.stem] = relation_columns(self.con, quote(path.stem))

    def add(self, severity: str, check: str, table: str, message: str, rows: int = 0) -> None:
        self.findings.append(Finding(severity, check, table, message, rows))

    def count(self, sql: str) -> int:
        row = self.con.execute(sql).fetchone()
        return int(row[0]) if row and row[0] is not None else 0

    def has(self, table: str, *columns: str) -> bool:
        return table in self.present and all(c in self.present[table] for c in columns)

    def flag(self, severity: str, check: str, table: str, where: str, message: str) -> None:
        """Record a finding when any row of the table satisfies the condition."""
        rows = self.count(f"SELECT count(*) FROM {quote(table)} WHERE {where}")
        if rows:
            self.add(severity, check, table, message, rows)

    # ------------------------------------------------------------ structure

    def manifest(self) -> dict[str, Any]:
        path = self.dir / MANIFEST_NAME
        if not path.exists():
            self.add("critical", "manifest_present", "", "The manifest file is absent.")
            return {}
        try:
            manifest = cast(dict[str, Any], json.loads(path.read_text()))
        except json.JSONDecodeError as problem:
            self.add(
                "critical", "manifest_present", "", f"The manifest is not valid JSON: {problem}"
            )
            return {}
        listed = set(cast(dict[str, Any], manifest.get("tables", {})))
        for name in sorted(set(self.present) ^ listed):
            where = "the manifest" if name in listed else "the bundle folder"
            self.add("error", "manifest_matches_files", name, f"The table is only in {where}.")
        for name, body in cast(dict[str, Any], manifest.get("tables", {})).items():
            for column, entry in cast(dict[str, Any], body.get("columns", {})).items():
                if not entry.get("description"):
                    self.add(
                        "warning", "column_described", name, f"Column {column} has no description."
                    )
        return manifest

    def tables_present(self) -> None:
        for table in self.catalog.tables.values():
            if table.catalog_required and table.name not in self.present:
                self.add("critical", "required_table", table.name, "A required table is missing.")
        if not any(name in self.present for name in COUNT_TABLES):
            self.add("critical", "count_table", "", "The bundle holds no count table.")
        for name in self.present:
            if name not in self.catalog.tables and not name.startswith("x_"):
                self.add(
                    "error", "table_name", name, "A table outside the catalog must start with x_."
                )

    def columns_present(self, table: Table) -> bool:
        """Critical checks on one table. False means row checks cannot run."""
        found, usable = self.present[table.name], True
        for column in table.columns:
            if column.name not in found:
                if column.required:
                    self.add(
                        "critical",
                        "required_column",
                        table.name,
                        f"Column {column.name} is missing.",
                    )
                    usable = False
            elif found[column.name] != column.sql_type:
                self.add(
                    "critical",
                    "column_type",
                    table.name,
                    f"Column {column.name} is {found[column.name]}, not {column.sql_type}.",
                )
                usable = False
        if usable and table.identifier:
            ident = quote(table.identifier)
            rows = self.count(
                f"SELECT coalesce(sum(n), 0) FROM (SELECT count(*) AS n FROM {quote(table.name)} "
                f"GROUP BY {ident} HAVING count(*) > 1 OR {ident} IS NULL)"
            )
            if rows:
                self.add(
                    "critical",
                    "identifier_unique",
                    table.name,
                    f"Identifier {table.identifier} is empty or appears more than once.",
                    rows,
                )
        return usable

    # ------------------------------------------------------------ row checks

    def rows(self, table: Table) -> None:
        name, found = table.name, self.present[table.name]
        for column in table.columns:
            if column.name not in found:
                continue
            col = quote(column.name)
            if column.required:
                self.flag(
                    "error",
                    "required_value",
                    name,
                    f"{col} IS NULL",
                    f"Column {column.name} is empty.",
                )
            if column.enum:
                allowed = ", ".join(literal(v) for v in column.enum)
                self.flag(
                    "error",
                    "allowed_value",
                    name,
                    f"{col} IS NOT NULL AND {col} NOT IN ({allowed})",
                    f"Column {column.name} holds a value outside its allowed list.",
                )
            for bound, op in ((column.minimum, "<"), (column.maximum, ">")):
                if bound is not None:
                    self.flag(
                        "error",
                        "value_range",
                        name,
                        f"{col} {op} {bound}",
                        f"Column {column.name} holds a value {op} {bound:g}.",
                    )
            if column.references:
                self.reference(name, column.name, column.references)
        if table.key and not table.identifier and all(k in found for k in table.key):
            keys = ", ".join(quote(k) for k in table.key)
            rows = self.count(
                f"SELECT coalesce(sum(n), 0) FROM (SELECT count(*) AS n FROM {quote(name)} "
                f"GROUP BY {keys} HAVING count(*) > 1)"
            )
            if rows:
                self.add(
                    "error",
                    "key_unique",
                    name,
                    f"Rows repeat the key ({', '.join(table.key)}).",
                    rows,
                )
        for group in table.exactly_one:
            filled = " + ".join(
                f"({quote(c)} IS NOT NULL)::INT" if c in found else "0" for c in group
            )
            self.flag(
                "error",
                "exactly_one",
                name,
                f"({filled}) <> 1",
                f"Exactly one of {', '.join(group)} must be filled.",
            )

    def reference(self, name: str, column: str, target: str) -> None:
        col = quote(column)
        target_table = self.catalog.tables[target]
        if not target_table.identifier or not self.has(target, target_table.identifier):
            self.flag(
                "error",
                "reference_exists",
                name,
                f"{col} IS NOT NULL",
                f"Column {column} refers to table {target}, which the bundle does not hold.",
            )
            return
        self.flag(
            "error",
            "reference_exists",
            name,
            f"{col} IS NOT NULL AND {col} NOT IN "
            f"(SELECT {quote(target_table.identifier)} FROM {quote(target)})",
            f"Column {column} points at a row of {target} that does not exist.",
        )

    # ------------------------------------------------------------ specific checks

    def taxon_assignment(self) -> None:
        name = "taxon_assignment"
        if not self.has(name, "outcome", "taxon_id", "run_id"):
            return
        self.flag(
            "error",
            "outcome_consistent",
            name,
            "(outcome = 'none' AND taxon_id IS NOT NULL) OR (outcome = 'assigned' AND taxon_id IS NULL)",
            "Outcome and taxon disagree: none must have no taxon, assigned must have one.",
        )
        if not self.has(name, "is_preferred") or not self.has("run", "run_id", "parent_run_id"):
            return
        kinds = [
            c for c in ("sequence_id", "genome_id", "cluster_id", "read_id") if self.has(name, c)
        ]
        subject = "coalesce(" + ", ".join(f"'{c}:' || {c}" for c in kinds) + ")"
        rows = self.count(
            "WITH RECURSIVE up(run_id, root) AS ("
            " SELECT run_id, run_id FROM run WHERE parent_run_id IS NULL"
            " UNION ALL SELECT r.run_id, up.root FROM run r JOIN up ON r.parent_run_id = up.run_id) "
            f"SELECT count(*) FROM (SELECT {subject} AS subject, root FROM {name} t "
            "JOIN up USING (run_id) WHERE is_preferred GROUP BY 1, 2 HAVING count(*) > 1)"
        )
        if rows:
            self.add(
                "error",
                "one_preferred",
                name,
                "A thing has two preferred rows from one pipeline.",
                rows,
            )

    def alignment(self) -> None:
        for side in ("query", "target"):
            if self.has("alignment", f"{side}_start", f"{side}_end"):
                self.flag(
                    "error",
                    "coordinate_order",
                    "alignment",
                    f"{side}_end < {side}_start",
                    f"{side}_end is less than {side}_start.",
                )

    def read_totals(self) -> None:
        if not self.has("read_count", "library_id", "run_id", "count", "count_unit"):
            return
        if not self.has("library", "library_id", "fragments_after_qc", "reads_after_qc"):
            return
        totals = (
            "SELECT c.library_id, c.run_id, sum(c.count) AS counted, any_value("
            "CASE c.count_unit WHEN 'fragments' THEN l.fragments_after_qc ELSE l.reads_after_qc END"
            ") AS total FROM read_count c JOIN library l USING (library_id) GROUP BY ALL"
        )
        for severity, check, where, message in (
            ("error", "counts_within_total", "counted > total", "Counts exceed the library total."),
            (
                "warning",
                "counts_reach_total",
                "counted < total",
                "Counts fall short of the library total.",
            ),
            (
                "warning",
                "library_total_known",
                "total IS NULL",
                "The library records no total to check counts against.",
            ),
        ):
            rows = self.count(f"SELECT count(*) FROM ({totals}) WHERE {where}")
            if rows:
                self.add(severity, check, "read_count", f"{message} Rows are libraries.", rows)

    def digests(self) -> None:
        if not self.has("sequence", "sequence_id", "bases", "kind"):
            return
        mate = "mate_bases" if self.has("sequence", "mate_bases") else "NULL"
        rows = self.con.execute(
            f"SELECT sequence_id, bases, {mate} FROM sequence WHERE bases IS NOT NULL "
            f"USING SAMPLE {DIGEST_SAMPLE} ROWS"
        ).fetchall()
        wrong = 0
        for found, bases, mate_bases in rows:
            expected = ids.read_pair_id(bases, mate_bases) if mate_bases else ids.sequence_id(bases)
            wrong += found != expected
        if wrong:
            self.add(
                "error",
                "identifier_matches_bases",
                "sequence",
                f"Identifier does not match the digest of the bases, in a sample of {len(rows)} rows.",
                wrong,
            )

    def completeness(self) -> None:
        if self.has("sequence", "bases", "kind"):
            self.flag(
                "warning",
                "bases_present",
                "sequence",
                "bases IS NULL AND kind <> 'reference'",
                "A read cluster or contig has no bases.",
            )
        if self.has("run", "run_id", "step"):
            params = (
                "SELECT run_id FROM run_parameter"
                if self.has("run_parameter", "run_id")
                else "SELECT NULL WHERE false"
            )
            self.flag(
                "warning",
                "run_has_parameters",
                "run",
                f"step <> 'pipeline' AND run_id NOT IN ({params})",
                "A run lists no parameters.",
            )
            kept = (
                "SELECT run_id FROM run_parameter WHERE name = 'hits_kept'"
                if self.has("run_parameter", "run_id", "name")
                else "SELECT NULL WHERE false"
            )
            self.flag(
                "warning",
                "hits_kept_declared",
                "run",
                f"step = 'search' AND run_id NOT IN ({kept})",
                "A search run does not declare which hits were kept.",
            )
        if self.has("reference_release", "checksum"):
            self.flag(
                "warning",
                "release_checksum",
                "reference_release",
                "checksum IS NULL",
                "A reference release has no checksum.",
            )

    def run(self) -> list[Finding]:
        self.manifest()
        self.tables_present()
        for name in self.present:
            table = self.catalog.tables.get(name)
            if table and self.columns_present(table):
                self.rows(table)
        self.taxon_assignment()
        self.alignment()
        self.read_totals()
        self.digests()
        self.completeness()
        return sorted(self.findings, key=lambda f: (SEVERITIES.index(f.severity), f.table, f.check))


def validate(bundle_dir: Path) -> list[Finding]:
    """Run every check on a bundle and return the findings, worst first."""
    return _Checker(bundle_dir).run()


def conforms(findings: list[Finding]) -> bool:
    """A bundle conforms when it has no critical failures and no errors."""
    return not any(f.severity in ("critical", "error") for f in findings)


def format_report(findings: list[Finding]) -> str:
    lines = [
        f"{f.severity.upper():8} {f.table or '(bundle)':18} {f.message}"
        + (f" [{f.failing_rows:,} rows]" if f.failing_rows else "")
        for f in findings
    ]
    tally = {s: sum(f.severity == s for f in findings) for s in SEVERITIES}
    verdict = "conforms" if conforms(findings) else "does not conform"
    summary = ", ".join(f"{n} {s}" for s, n in tally.items())
    return "\n".join([*lines, f"Result: {verdict} ({summary})."])
