# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Construction of sub-regions from detected subclusters.

Sub-regions built here are added to the record during the region-formation step
of the main pipeline, letting subclusters extend existing regions or form new
ones. How far a detected subcluster is allowed to alter region boundaries
depends on the mode used:
 - "clip": only the parts overlapping an area found by another detection module
           are kept, truncated to that area, so regions can never grow
 - "extend": overlapping subclusters are kept in full, so they can extend an
             existing region, but subclusters without any overlap are discarded
 - "create": every subcluster is kept in full, so they can also create new regions
"""
from enum import StrEnum, auto
from typing import Sequence

from antismash.common.secmet import Record, SubRegion
from antismash.common.secmet.features import CDSCollection
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
)

# the label given to every sub-region created by this module
LABEL = "subclusters"


class SubRegionMode(StrEnum):
    CLIP = auto()
    EXTEND = auto()
    CREATE = auto()


def gather_foreign_areas(record: Record) -> list[CDSCollection]:
    """Protoclusters and subregions on the record from other detection modules.

    This module's own protoclusters are never added to the record, so everything
    found here belongs to a module that has already run.
    """
    areas: list[CDSCollection] = list(record.get_protoclusters())
    areas.extend(record.get_subregions())
    return areas


def build_subregions(record: Record, areas: Sequence[CDSCollection], *, tool: str,
                     mode: SubRegionMode) -> list[SubRegion]:
    """Build the sub-region features for a set of detected subclusters.

    Arguments:
        record: the record the subclusters were detected in, which must already
            contain the areas found by any earlier detection module
        areas: the areas detected by this module, i.e. its protoclusters
        tool: the tool name to set on each resulting sub-region
        mode: how far the subclusters may alter region boundaries

    Returns:
        a list of sub-regions, which may be empty
    """
    if mode not in SubRegionMode:
        raise ValueError(f"Unknown subcluster subregion mode: {mode!r}")

    if not areas:
        return []

    subregions = [
        SubRegion(location, tool=tool, label=LABEL)
        for location, _ in record.get_potential_regions(areas)
    ]

    # "create": every subcluster is used for subregion formation
    if mode == SubRegionMode.CREATE:
        return subregions

    # for "extend" and "clip", the merged areas of the other modules are needed
    external = [location for location, _
                in record.get_potential_regions(gather_foreign_areas(record))]

    # in both modes, only subclusters overlapping one of those areas are kept
    subregions = [subregion for subregion in subregions
                  if any(subregion.overlaps_with(area) for area in external)]

    # "extend": those overlapping subclusters are used in full
    if mode == SubRegionMode.EXTEND:
        return subregions

    # "clip": only the sections shared with an area from another module are kept,
    # so the truncated results are locations rather than protoclusters

    def clip_location(location: Location, limit: Location) -> Location | None:
        parts = list(location.parts)
        parts[0] = FeatureLocation(max(parts[0].start, limit.parts[0].start), parts[0].end, location.strand)
        parts[-1] = FeatureLocation(parts[-1].start, min(parts[-1].end, limit.parts[-1].end), location.strand)

        parts = [p for p in parts if len(p) > 1]
        if not parts:
            return None
        if len(parts) > 1:
            return CompoundLocation(parts)
        return parts[0]
        
    clipped = []
    for subregion in subregions:
        for other_location in external:
            if not subregion.overlaps_with(other_location):
                continue
            new_location = clip_location(subregion.location, other_location)
            if new_location is None:
                continue
            clipped.append(SubRegion(new_location, tool=tool, label=LABEL))

    return clipped