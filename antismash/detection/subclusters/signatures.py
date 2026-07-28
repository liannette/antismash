"""HMM signatures for the subcluster detection module."""

from types import MappingProxyType
from typing import Mapping, Optional

from antismash.common import path
from antismash.common.signature import HmmSignature


class SubclusterHmmSignature(HmmSignature):
    """An HMM signature extended with an optional Pfam accession."""

    def __init__(self, name: str, description: str, cutoff: int,
                 hmm_path: str, seed_count: int = 0, *,
                 accession: Optional[str] = None) -> None:
        super().__init__(name, description, cutoff, hmm_path, seed_count)
        self.accession = accession


_SIGNATURE_CACHE: list[SubclusterHmmSignature] = []
_SIGNATURE_BY_NAME_CACHE: dict[str, SubclusterHmmSignature] = {}


def _ensure_signatures_loaded() -> None:
    """Load the subcluster HMM signatures from disk into the caches, once.

    Both the ordered list and the name-keyed mapping are populated from the
    same read, so they can never drift out of sync.
    """
    if _SIGNATURE_CACHE:
        return
    filename = path.get_full_path(__file__, "data", "hmmdetails.txt")
    _SIGNATURE_CACHE.extend(_read_signatures(filename))
    _SIGNATURE_BY_NAME_CACHE.update((signature.name, signature) for signature in _SIGNATURE_CACHE)


def get_signature_profiles() -> list[SubclusterHmmSignature]:
    """Return all subcluster HMM signatures, loading from disk on first call."""
    _ensure_signatures_loaded()
    return list(_SIGNATURE_CACHE)


def get_signature_profiles_by_name() -> Mapping[str, SubclusterHmmSignature]:
    """Return all subcluster HMM signatures keyed by signature name.

    The mapping is built once (when the signatures are first loaded) and returned
    as a read-only view, so repeated calls stay O(1) regardless of how many
    signatures the module grows to contain.
    """
    _ensure_signatures_loaded()
    return MappingProxyType(_SIGNATURE_BY_NAME_CACHE)


def _read_signatures(detail_file: str) -> list[SubclusterHmmSignature]:
    """Parse a 5-column hmmdetails TSV into signature objects.

    Columns (tab-separated): name  description  cutoff  hmm_file  accession
    """
    bad_lines: list[str] = []
    signatures: list[SubclusterHmmSignature] = []
    with open(detail_file, "r", encoding="utf-8") as data:
        for line in data.read().split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            try:
                name, desc, cutoff, filename, accession = line.split("\t")
            except ValueError:
                bad_lines.append(line)
                continue
            signatures.append(SubclusterHmmSignature(
                name, 
                desc, 
                int(cutoff), 
                path.get_full_path(detail_file, filename),
                accession=accession)
            )

    if bad_lines:
        raise ValueError("Invalid lines in HMM detail file (first 10):\n%s" % "\n".join(bad_lines[:10]))

    return signatures