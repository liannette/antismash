"""HMM signatures for the subcluster detection module."""

from typing import Optional

from antismash.common import path
from antismash.common.signature import HmmSignature


class SubclusterHmmSignature(HmmSignature):
    """An HMM signature extended with an optional Pfam accession."""

    def __init__(self, name: str, description: str, cutoff: int,
                 hmm_path: str, seed_count: int = 0, *,
                 accession: Optional[str] = None) -> None:
        super().__init__(name, description, cutoff, hmm_path, seed_count)
        self.accession = accession


# the description of every signature used by the module
DETAILS_FILE = path.get_full_path(__file__, "data", "hmmdetails.txt")
# the combined profile database built from those signatures by prepare_data()
AGGREGATE_HMM_FILE = path.get_full_path(__file__, "data", "subclusters.hmm")

_SIGNATURE_CACHE: dict[str, SubclusterHmmSignature] = {}


def _ensure_signatures_loaded() -> None:
    """Load the subcluster HMM signatures from disk into the cache, once."""
    if _SIGNATURE_CACHE:
        return
    signatures = _read_signatures(DETAILS_FILE)
    _SIGNATURE_CACHE.update((signature.name, signature) for signature in signatures)
    if len(_SIGNATURE_CACHE) != len(signatures):
        raise ValueError(f"Duplicate signature names in {DETAILS_FILE}")


def get_signatures() -> dict[str, SubclusterHmmSignature]:
    """Return all subcluster HMM signatures, keyed by name, loading from disk
    on first call.

    The returned dict is a copy; modifying it does not affect the cache.
    """
    _ensure_signatures_loaded()
    return dict(_SIGNATURE_CACHE)


def get_signature(name: str) -> SubclusterHmmSignature:
    """Return a single subcluster HMM signature by name, loading from disk on
    first call.

    Arguments:
        name: the name of the signature

    Returns:
        the matching signature

    Raises:
        ValueError: if no signature with that name exists
    """
    _ensure_signatures_loaded()
    if name not in _SIGNATURE_CACHE:
        raise ValueError(f"Unknown subcluster signature {name!r}")
    return _SIGNATURE_CACHE[name]


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