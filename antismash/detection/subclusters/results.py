import logging
from collections import defaultdict
from dataclasses import dataclass
from enum import StrEnum, auto
from functools import cached_property
from typing import Any, Optional, Self, Sequence

from antismash.common.hmm_rule_parser.rule_parser import DetectionRule
from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults, RuleDetectionResults
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Record, Region, SubRegion
from antismash.common.secmet.features import CDSCollection
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
    location_contains_other,
    locations_overlap,
)

from .compounds import CompoundInfo, get_compound
from .ruleset import get_ruleset
from .signatures import get_signatures

# how far a detected subcluster is allowed to alter region boundaries when
# subclusters are added to the record as subregions:
#  - "clip": only the parts overlapping an area found by another detection module
#            are kept, truncated to that area, so regions can never grow
#  - "extend": overlapping subclusters are kept in full, so they can extend an
#              existing region, but subclusters without any overlap are discarded
#  - "any": every subcluster is kept in full, so they can also form new regions


class SubRegionMode(StrEnum):
    CLIP = auto()
    EXTEND = auto()
    ANY = auto()


@dataclass(frozen=True)
class CDSDomainHit:
    """A single signature match found in a specific CDS."""
    domain_name: str
    domain_description: Optional[str]
    domain_accession: Optional[str]
    cds_name: str
    evalue: float
    bitscore: float


class SubclusterPrediction:
    """A single predicted subcluster, used as a view object for HTML rendering.

    Attributes:
        rule: The detection rule whose conditions were met
        location: location of the protocluster that produced this prediction,
            i.e. the matching core along with any neighbourhood the rule defines
        cds_results: the per-CDS detection results for every CDS that contributed
            to this prediction, as returned by the detection pipeline
    """

    def __init__(
            self,
            *,
            rule: DetectionRule,
            location: FeatureLocation,
            cds_results: list[CDSResults],
    ) -> None:
        self.rule = rule
        self.location = location
        self.cds_results = cds_results

    @property
    def rule_name(self) -> str:
        """Name of the matching detection rule."""
        return self.rule.name

    @cached_property
    def compound(self) -> CompoundInfo:
        """The compound associated with the detection rule.

        Every rule must have a matching entry in the compound details file.
        """
        return get_compound(self.rule.name)

    @property
    def conditions_str(self) -> str:
        text = str(self.rule.conditions)
        if text.startswith("(") and text.endswith(")"):
            return text[1:-1]
        return text

    @cached_property
    def domain_hits(self) -> list[CDSDomainHit]:
        """A flat list of every domain hit, each paired with its CDS name."""
        hits: list[CDSDomainHit] = []
        for cds_result in self.cds_results:
            cds_name = cds_result.cds.get_name()
            fired = cds_result.definition_domains.get(self.rule_name, set())
            for domain_name in sorted(fired):
                matching = [d for d in cds_result.domains if d.name == domain_name]
                if not matching:
                    continue
                best = max(matching, key=lambda d: d.bitscore)
                signature = get_signatures()[domain_name]
                hits.append(CDSDomainHit(
                    domain_name=signature.name,
                    domain_description=signature.description,
                    domain_accession=signature.accession,
                    cds_name=cds_name,
                    evalue=best.evalue,
                    bitscore=best.bitscore,
                ))
        return hits

    @property
    def domain_hits_by_cds(self) -> dict[str, list[CDSDomainHit]]:
        hits_by_cds: dict[str, list[CDSDomainHit]] = defaultdict(list)
        for hit in self.domain_hits:
            hits_by_cds[hit.cds_name].append(hit)
        return dict(hits_by_cds)

    @property
    def domain_hits_by_domain(self) -> dict[str, list[CDSDomainHit]]:
        hits_by_domain: dict[str, list[CDSDomainHit]] = defaultdict(list)
        for hit in self.domain_hits:
            hits_by_domain[hit.domain_name].append(hit)
        return dict(hits_by_domain)

    def __repr__(self) -> str:
        return (
            f"SubclusterPrediction(rule_name={self.rule_name!r}, "
            f"location={self.location.start}-{self.location.end}, "
            f"cds_count={len(self.cds_results)})"
        )


class SubclusterDetectionResults(DetectionResults):
    """Results class for the Subcluster detection module """

    schema_version = 1  # increment when the JSON format changes

    def __init__(
            self,
            record_id: str,
            rule_results: RuleDetectionResults,
            rule_names: set[str],
            strictness: str,
            record: Record,
            subregion_mode: SubRegionMode = SubRegionMode.CLIP,
    ) -> None:

        if subregion_mode not in SubRegionMode:
            raise ValueError(f"Unknown subcluster subregion mode: {subregion_mode!r}")

        super().__init__(record_id)
        self.rule_results = rule_results
        self.rule_names = rule_names
        self.strictness = strictness
        self.subregion_mode = subregion_mode
        self._record = record

        ruleset = get_ruleset(self.strictness)
        self.predictions = [
            SubclusterPrediction(
                rule=ruleset.get_rule_by_name(protocluster.product),
                location=protocluster.location,
                cds_results=cds_results,
            )
            for protocluster, cds_results in rule_results.cds_by_cluster.items()
        ]

    def get_predictions_for_region(self, region: Region) -> list[SubclusterPrediction]:
        """Return all predictions fully contained within the given region."""
        return [
            prediction for prediction in self.predictions
            if location_contains_other(region.location, prediction.location)
        ]

    def get_predicted_subregions(self) -> list[SubRegion]:
        """Return sub-region features for the detected subclusters.

        These are added to the record during the region-formation step of the
        main pipeline, letting subclusters extend existing regions or form new
        ones. How far they may alter region boundaries depends on the configured
        subregion mode.
        """
        if not self.predictions:
            return []

        current: Sequence[CDSCollection] = self.rule_results.protoclusters
        subregions = [
            SubRegion(location, tool=self.rule_results.tool, label="subclusters")
            for location, _ in self._record.get_potential_regions(current)
        ]

        # "any": use each subcluster prediction for subregion formation
        if self.subregion_mode == SubRegionMode.ANY:
            return subregions

        # get merged locations of other clusters for "extend" and "clip"
        external = [location for location, _ in
                    self._record.get_potential_regions(self._get_foreign_areas())]

        subregions = list(filter(lambda x: any(x.overlaps_with(extern) for extern in external), subregions))

        # "extend": only subcluster predictions that overlap with
        # an area from another module are used for for subregion formation
        if self.subregion_mode == SubRegionMode.EXTEND:
            return subregions

        assert self.subregion_mode == SubRegionMode.CLIP

        # "clip": only the sections shared with an area from another module are
        # kept, so the truncated results are locations rather than protoclusters

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
            areas: list[Location] = []
            for area in external:
                if subregion.overlaps_with(area):
                    areas.append(area)
            for location in areas:
                new_location = clip_location(subregion.location, location)
                if new_location is None:
                    continue
                clipped.append(SubRegion(new_location, tool=self.rule_results.tool,
                                         label="subclusters"))

        return clipped

    def _get_foreign_areas(self) -> list[CDSCollection]:
        """Protoclusters and subregions on the record from other detection modules.
        """
        areas: list[CDSCollection] = list(self._record.get_protoclusters())
        areas.extend(self._record.get_subregions())
        return areas

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "strictness": self.strictness,
            "subregion_mode": self.subregion_mode,
            "rule_names": sorted(self.rule_names),
            "rule_results": self.rule_results.to_json(),
        }

    @classmethod
    def from_json(cls, data: dict[str, Any], record: Record) -> Optional[Self]:
        if data.get("schema_version") != cls.schema_version:
            logging.debug(
                "Discarding subcluster results: schema version %s != %s",
                data.get("schema_version"), cls.schema_version,
            )
            return None
        
        rule_results = RuleDetectionResults.from_json(data["rule_results"], record)
        if rule_results is None:
            logging.debug("Discarding subcluster results: rule detection result schema changed")
            return None
        
        return cls(
            record_id=data["record_id"],
            rule_results=rule_results,
            rule_names=set(data["rule_names"]),
            strictness=data["strictness"],
            subregion_mode=data["subregion_mode"],
            record=record,
        )

    def add_to_record(self, record: Record) -> None:
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        self.rule_results.annotate_cds_features()


def _intersect_locations(outer: Location, inner: Location,
                         wrap_point: Optional[int]) -> Optional[Location]:
    """Truncate a location to the section it shares with another.

    Arguments:
        outer: the location to truncate to
        inner: the location to truncate
        wrap_point: the record length for circular records, otherwise None

    Returns:
        the shared location, or None if the two don't overlap
    """
    parts = []
    for outer_part in outer.parts:
        for inner_part in inner.parts:
            if not locations_overlap(outer_part, inner_part):
                continue
            parts.append(FeatureLocation(max(outer_part.start, inner_part.start),
                                         min(outer_part.end, inner_part.end), 1))
    if not parts:
        return None
    if len(parts) == 1:
        return parts[0]
    # with both locations crossing the origin, the shared section does too, and
    # by convention is ordered with the piece before the origin first
    if wrap_point is not None:
        before = [part for part in parts if part.end == wrap_point]
        after = [part for part in parts if part.start == 0]
        if len(before) == 1 and len(after) == 1:
            return CompoundLocation([before[0], after[0]])
    # otherwise the two only share separate sections, which a single subregion
    # cannot cover, so keep the largest of them
    return max(parts, key=len)
