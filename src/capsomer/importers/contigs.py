"""Import assembled contigs from a FASTA file.

Any assembler's output can be read. The depth an assembler prints in the
contig header is kept when the header follows the MEGAHIT or SPAdes form.
"""

from collections.abc import Mapping, Sequence
from pathlib import Path

from capsomer.builder import Builder
from capsomer.bundle import literal
from capsomer.importers.common import load_fasta

DEPTH_PATTERNS = (r"multi=([0-9.]+)", r"_cov_([0-9.]+)")


def add_contigs(
    builder: Builder,
    fasta: Path,
    *,
    run_id: str,
    tool: str,
    tool_version: str | None = None,
    libraries: Sequence[str] = (),
    inputs: Sequence[str] = (),
    parent: str | None = None,
    parameters: Mapping[str, object] | None = None,
) -> int:
    """Add the contigs of one assembly run and return how many the file holds.

    A contig whose bases the bundle already holds adds no sequence row. It
    still adds a sequence_source row, so the bundle records that this run
    produced it too.
    """
    count = load_fasta(builder.con, fasta)
    builder.add_run(
        run_id,
        tool=tool,
        tool_version=tool_version,
        step="assembly",
        description="Assembly of reads into contigs.",
        parent=parent,
        parameters=parameters,
        libraries=libraries,
        inputs=inputs,
    )
    builder.add(
        "sequence",
        "SELECT sequence_id, 'contig' AS kind, length, bases, gc_content FROM _fasta",
    )
    depth = "coalesce(" + ", ".join(
        f"nullif(regexp_extract(header, {literal(p)}, 1), '')" for p in DEPTH_PATTERNS
    )
    builder.add(
        "sequence_source",
        f"SELECT sequence_id, {literal(run_id)} AS run_id, name, {depth}) AS assembler_depth FROM _fasta",
    )
    return count
