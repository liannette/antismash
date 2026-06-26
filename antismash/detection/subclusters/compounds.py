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
        filename = path.get_full_path(__file__, "data", "compounds.tsv")
        _COMPOUND_CACHE.update(_read_compounds(filename))
    return dict(_COMPOUND_CACHE)


def _read_compounds(detail_file: str) -> dict[str, CompoundInfo]:
    #TODO: read substructure compound metadata from file
    return NotImplementedError()