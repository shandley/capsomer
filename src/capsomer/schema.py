"""Load the Capsomer catalog from its LinkML schema file."""

import re
from dataclasses import dataclass
from functools import cache
from importlib import resources
from typing import Any, cast

import yaml

SHAPES = ("thing", "count", "label", "link")
SQL_TYPES = {
    "string": "VARCHAR",
    "integer": "BIGINT",
    "float": "DOUBLE",
    "boolean": "BOOLEAN",
    "date": "DATE",
}


@dataclass(frozen=True)
class Column:
    name: str
    sql_type: str
    required: bool
    description: str
    references: str | None
    enum: tuple[str, ...] | None
    minimum: float | None
    maximum: float | None
    annotations: dict[str, str]


@dataclass(frozen=True)
class Table:
    name: str
    shape: str
    row: str
    description: str
    catalog_required: bool
    identifier: str | None
    key: tuple[str, ...]
    exactly_one: tuple[tuple[str, ...], ...]
    columns: tuple[Column, ...]

    def column(self, name: str) -> Column | None:
        return next((c for c in self.columns if c.name == name), None)

    @property
    def column_names(self) -> tuple[str, ...]:
        return tuple(c.name for c in self.columns)


@dataclass(frozen=True)
class Catalog:
    version: str
    tables: dict[str, Table]


def table_name(class_name: str) -> str:
    """ReadCount -> read_count."""
    return re.sub(r"(?<!^)(?=[A-Z])", "_", class_name).lower()


def schema_text() -> str:
    return resources.files("capsomer").joinpath("schema/capsomer.yaml").read_text("utf-8")


def _annotations(raw: dict[str, Any]) -> dict[str, str]:
    found = cast(dict[str, Any], raw.get("annotations") or {})
    return {key: str(value) for key, value in found.items()}


def _column(name: str, raw: dict[str, Any], enums: dict[str, Any], classes: set[str]) -> Column:
    range_name = str(raw.get("range", "string"))
    references: str | None = None
    enum: tuple[str, ...] | None = None
    if range_name in enums:
        values = cast(dict[str, Any], enums[range_name]["permissible_values"])
        enum, sql_type = tuple(str(v) for v in values), "VARCHAR"
    elif range_name in classes:
        references, sql_type = table_name(range_name), "VARCHAR"
    else:
        sql_type = SQL_TYPES[range_name]
    minimum, maximum = raw.get("minimum_value"), raw.get("maximum_value")
    return Column(
        name=name,
        sql_type=sql_type,
        required=bool(raw.get("required", False)),
        description=str(raw.get("description", "")),
        references=references,
        enum=enum,
        minimum=None if minimum is None else float(minimum),
        maximum=None if maximum is None else float(maximum),
        annotations=_annotations(raw),
    )


def _table(class_name: str, raw: dict[str, Any], enums: dict[str, Any], classes: set[str]) -> Table:
    notes = _annotations(raw)
    attributes = cast(dict[str, dict[str, Any] | None], raw["attributes"])
    columns = tuple(_column(n, a or {}, enums, classes) for n, a in attributes.items())
    identifier = next((n for n, a in attributes.items() if a and a.get("identifier")), None)
    unique = cast(dict[str, Any], raw.get("unique_keys") or {})
    key: tuple[str, ...] = (identifier,) if identifier else ()
    if "primary" in unique:
        key = tuple(str(s) for s in unique["primary"]["unique_key_slots"])
    groups = tuple(
        tuple(group.split()) for group in notes.get("exactly_one", "").split(";") if group.strip()
    )
    if notes.get("shape") not in SHAPES:
        raise ValueError(f"{class_name}: shape must be one of {SHAPES}")
    return Table(
        name=table_name(class_name),
        shape=notes["shape"],
        row=notes["row"],
        description=" ".join(str(raw.get("description", "")).split()),
        catalog_required=notes.get("catalog_required", "False") == "True",
        identifier=identifier,
        key=key,
        exactly_one=groups,
        columns=columns,
    )


@cache
def load_catalog() -> Catalog:
    """Read the packaged LinkML schema into table definitions."""
    raw = cast(dict[str, Any], yaml.safe_load(schema_text()))
    enums = cast(dict[str, Any], raw["enums"])
    classes = cast(dict[str, dict[str, Any]], raw["classes"])
    names = set(classes)
    tables = [_table(name, body, enums, names) for name, body in classes.items()]
    return Catalog(version=str(raw["version"]), tables={t.name: t for t in tables})
