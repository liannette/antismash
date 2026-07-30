import logging
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Optional, Self

from antismash.common.hmm_rule_parser.rule_parser import DetectionRule
from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults, RuleDetectionResults
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Protocluster, Record, Region, SubRegion
from antismash.common.secmet.locations import (
    FeatureLocation,
    location_contains_other,
    locations_overlap,
)

from .compounds import CompoundInfo, get_compound
from .ruleset import get_ruleset
from .signatures import get_signatures


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
        core_location: location of the protocluster core that produced this prediction
        cds_results: the per-CDS detection results for every CDS that contributed
            to this prediction, as returned by the detection pipeline
    """

    def __init__(
            self,
            *,
            rule: DetectionRule,
            core_location: FeatureLocation,
            cds_results: list[CDSResults],
    ) -> None:
        self.rule = rule
        self.core_location = core_location
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
                best = max(matching, key=lambda d: d.bitscore) if matching else None
                signature = get_signatures()[domain_name]
                hits.append(CDSDomainHit(
                    domain_name=signature.name,
                    domain_description=signature.description,
                    domain_accession=signature.accession,
                    cds_name=cds_name,
                    evalue=best.evalue if best else None,
                    bitscore=best.bitscore if best else None,
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
            f"core={self.core_location.start}-{self.core_location.end}, "
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
            as_subregions: bool = False,
            require_overlap: bool = False,
            record: Optional[Record] = None,
    ) -> None:
        super().__init__(record_id)
        self.rule_results = rule_results
        self.rule_names = rule_names
        self.strictness = strictness
        # when True, each detected subcluster protocluster is exposed as a
        # sub-region feature so that it is part of the region formation step
        self.as_subregions = as_subregions
        # when True, only subclusters that overlap a cluster found by another
        # detection module are emitted, so subclusters can extend existing
        # regions but never form standalone regions (or merge with each other)
        self.require_overlap = require_overlap
        # the record being analysed, needed to look up clusters from other
        # detection modules when require_overlap is enabled
        self._record = record
        ruleset = get_ruleset(self.strictness)
        self.predictions = [
            SubclusterPrediction(
                rule=ruleset.get_rule_by_name(protocluster.product),
                core_location=protocluster.core_location,
                cds_results=cds_results,
            )
            for protocluster, cds_results in rule_results.cds_by_cluster.items()
        ]

    def get_predictions_for_region(self, region: Region) -> list[SubclusterPrediction]:
        """Return all predictions fully contained within the given region."""
        return [
            prediction for prediction in self.predictions
            if location_contains_other(region.location, prediction.core_location)
        ]

    def get_predictions_outside_regions(self, record: Record) -> list[SubclusterPrediction]:
        """Return predictions not fully contained by any region in the record."""
        return [
            prediction for prediction in self.predictions
            if not any(
                location_contains_other(region.location, prediction.core_location)
                for region in record.get_regions()
            )
        ]

    def get_predicted_subregions(self) -> list[SubRegion]:
        """Return each detected subcluster as a sub-region feature.

        When enabled by the corresponding option, these are added to the record
        during the region-formation step of the main pipeline, letting subclusters
        form new regions or extend existing ones. When disabled, an empty list is
        returned and subclusters remain display-only annotations within existing
        regions.
        """
        if not self.as_subregions:
            return []
        subregions = [
            SubRegion(protocluster.location, tool=self.rule_results.tool,
                      label=protocluster.product)
            for protocluster in self.rule_results.protoclusters
        ]
        if not self.require_overlap:
            return subregions

        # keep only subclusters that overlap with a cluster from another detection
        # module or a subregion; the rest are dropped so they neither create new 
        # regions nor merge with one another
        existing = [feature.location for feature in self._get_foreign_clusters()]
        return [
            subregion for subregion in subregions
            if any(locations_overlap(subregion.location, location) for location in existing)
        ]

    def _get_foreign_clusters(self) -> list[Protocluster | SubRegion]:
        """Protoclusters and subregions on the record from other detection modules.

        Used when overlap is required, to decide which subclusters may extend
        an existing region. Region features do not exist yet at the point this
        runs, so overlap is tested against the protoclusters and subregions that
        other detection modules have already added to the record.
        """
        if self._record is None:
            return []
        features: list[Protocluster | SubRegion] = []
        features.extend(self._record.get_protoclusters())
        features.extend(self._record.get_subregions())
        return [feature for feature in features if feature.tool != self.rule_results.tool]

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "strictness": self.strictness,
            "as_subregions": self.as_subregions,
            "require_overlap": self.require_overlap,
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
            as_subregions=data["as_subregions"],
            require_overlap=data["require_overlap"],
            record=record,
        )

    def add_to_record(self, record: Record) -> None:
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        self.rule_results.annotate_cds_features()
