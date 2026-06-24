"""HMM signature profiles for the subcluster detection module."""

from typing import Optional

from antismash.common import path
from antismash.common.signature import HmmSignature


class SubclusterHmmSignature(HmmSignature):
    """An HMM signature extended with an optional Pfam accession."""

    def __init__(self, name: str, description: str, cutoff: int,
                 hmm_path: str, accession: Optional[str] = None,
                 seed_count: int = 0) -> None:
        super().__init__(name, description, cutoff, hmm_path, seed_count)
        self.accession = accession


_PROFILE_CACHE: dict[str, SubclusterHmmSignature] = {}


def get_subcluster_profiles() -> dict[str, SubclusterHmmSignature]:
    """Return all subcluster HMM profiles, loading from disk on first call."""
    if not _PROFILE_CACHE:
        filename = path.get_full_path(__file__, "data", "hmmdetails.txt")
        _PROFILE_CACHE.update(_read_profiles(filename))
    return dict(_PROFILE_CACHE)


def _read_profiles(detail_file: str) -> dict[str, SubclusterHmmSignature]:
    """Parse a 5-column hmmdetails TSV into SubclusterHmmSignature objects.

    Columns (tab-separated): name  description  cutoff  hmm_file  [accession]
    The accession column is optional; omitted or blank entries get accession=None.
    """
    bad_lines: list[str] = []
    profiles: dict[str, SubclusterHmmSignature] = {}
    with open(detail_file, "r", encoding="utf-8") as data:
        for line in data.read().split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            try:
                parts = line.split("\t")
                name, desc, cutoff, filename = parts[:4]
                if len(parts) > 4 and parts[4].strip():
                    accession = parts[4].strip()
                else:
                    accession = None
            except ValueError:
                bad_lines.append(line)
                continue
            profiles[name] = SubclusterHmmSignature(
                name, desc, int(cutoff), path.get_full_path(detail_file, filename), 
                accession=accession)

    if bad_lines:
        raise ValueError("Invalid lines in HMM detail file (first 10):\n%s" % "\n".join(bad_lines[:10]))

    return profiles
