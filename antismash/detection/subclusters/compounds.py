"""Compound metadata for the subcluster detection module."""

from dataclasses import dataclass
from typing import Optional

from antismash.common import path


@dataclass(frozen=True)
class CompoundInfo:
    """Chemical metadata associated with a subcluster rule."""
    name: str
    smiles: Optional[str]
    classification: list[str]


_COMPOUND_CACHE: dict[str, CompoundInfo] = {}


def get_subcluster_compounds() -> dict[str, CompoundInfo]:
    """Return all subcluster compound metadata, loading from disk on first call."""
    if not _COMPOUND_CACHE:
        filename = path.get_full_path(__file__, "data", "compound_details.txt")
        _COMPOUND_CACHE.update(_read_compounds(filename))
    return dict(_COMPOUND_CACHE)


def _read_compounds(detail_file: str) -> dict[str, CompoundInfo]:
    """Parse compound_details.txt into a dict keyed by rule name.

    Columns (tab-separated): rule_name  compound_name  smiles  classification
    Classification is semicolon-separated. Lines starting with # are comments.
    """
    compounds: dict[str, CompoundInfo] = {}
    with open(detail_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError(f"Invalid line in {detail_file}: {line!r}")
            rule_name, name, smiles, classification_raw = parts[:4]
            compounds[rule_name] = CompoundInfo(
                name=name,
                smiles=smiles or None,
                classification=[c.strip() for c in classification_raw.split(";") if c.strip()],
            )
    return compounds