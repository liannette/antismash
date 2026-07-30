"""HMM signatures for the subcluster detection module."""

from antismash.common import path
from antismash.common.signature import HmmSignature


# the description of every signature used by the module
DETAILS_FILE = path.get_full_path(__file__, "data", "hmmdetails.txt")
# the combined profile database built from those signatures by prepare_data()
AGGREGATE_HMM_FILE = path.get_full_path(__file__, "data", "subcluster_seeds.hmm")

_SIGNATURE_CACHE: dict[str, "SubclusterHmmSignature"] = {}


class SubclusterHmmSignature(HmmSignature):
    """An HMM signature extended an accession."""

    def __init__(self, name: str, description: str, cutoff: int,
                 hmm_path: str, seed_count: int = 0, *,
                 accession: str) -> None:
        super().__init__(name, description, cutoff, hmm_path, seed_count, 
                         accession)


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
    """
    _ensure_signatures_loaded()
    return _SIGNATURE_CACHE


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