"""Capsomer: a data standard for virome analysis results."""

from capsomer.schema import Catalog, Column, Table, load_catalog

__version__ = "0.1.0"
__all__ = ["Catalog", "Column", "Table", "__version__", "load_catalog"]
