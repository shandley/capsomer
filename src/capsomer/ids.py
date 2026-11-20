"""Content-derived identifiers.

A sequence identifier is the GA4GH sequence digest of the bases as written:
the first 24 bytes of the SHA-512 hash of the upper-case sequence, base64url
encoded, with the prefix ``SQ.``. Capsomer adds no strand or rotation rule.
"""

import base64
import hashlib
from collections.abc import Iterable

SEQUENCE_PREFIX = "SQ."
READ_PAIR_PREFIX = "RP."
GENOME_PREFIX = "GN."


def sha512t24u(data: bytes) -> str:
    """GA4GH truncated digest: SHA-512, first 24 bytes, base64url."""
    return base64.urlsafe_b64encode(hashlib.sha512(data).digest()[:24]).decode("ascii")


def sequence_id(bases: str) -> str:
    """Identifier of one sequence, taken over the bases as written."""
    return SEQUENCE_PREFIX + sha512t24u(bases.upper().encode("ascii"))


def read_pair_id(first: str, second: str) -> str:
    """Identifier of a read pair kept as one unit: both reads, as sequenced."""
    joined = f"{first.upper()},{second.upper()}"
    return READ_PAIR_PREFIX + sha512t24u(joined.encode("ascii"))


def genome_id(member_sequence_ids: Iterable[str]) -> str:
    """Identifier of a genome, from the sorted identifiers of its sequences."""
    members = sorted(set(member_sequence_ids))
    if not members:
        raise ValueError("a genome needs at least one member sequence")
    return GENOME_PREFIX + sha512t24u(",".join(members).encode("ascii"))


def gc_content(bases: str) -> float | None:
    """Fraction of bases that are G or C, or None for an empty sequence."""
    if not bases:
        return None
    upper = bases.upper()
    return (upper.count("G") + upper.count("C")) / len(upper)
