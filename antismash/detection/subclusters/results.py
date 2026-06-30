import logging
from dataclasses import dataclass, field
from  typing import Any, Optional, Self

from antismash.common.hmm_rule_parser.rule_parser import DetectionRule
from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults, RuleDetectionResults
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Record, Region
from antismash.common.secmet.locations import FeatureLocation, location_contains_other

from .compounds import CompoundInfo, get_subcluster_compounds
from .ruleset import get_ruleset
from .signatures import SubclusterHmmSignature, get_subcluster_profiles


@dataclass(frozen=True)
class DomainHit:
    """A single profile match that contributed to a subcluster hit."""
    profile: SubclusterHmmSignature
    cds_name: str
    evalue: float
    bitscore: float


@dataclass
class SubclusterPrediction:
    """A single predicted subcluster.

    After detection, call ``enrich(rule, profiles, compounds)`` once to populate
    the derived fields used by the HTML template.

    Attributes:
        rule_name: Name of the ``DetectionRule`` that fired.
        start: Start of the core location (0-based bp).
        end: End of the core location (0-based bp).
        cds_results: ``CDSResults`` instances for every CDS that contributed
            to this prediction, as returned by the rule-based detection pipeline.
    """
    rule_name: str
    start: int
    end: int
    cds_results: list[CDSResults]

    # Lazy fields — populated on first access; not part of __init__
    _rule: Optional[DetectionRule] = field(default=None, init=False, repr=False)
    _domain_hits: Optional[list[DomainHit]] = field(default=None, init=False, repr=False)
    _compound: Optional[CompoundInfo] = field(default=None, init=False, repr=False)

    @property
    def rule(self) -> DetectionRule:
        """The ``DetectionRule`` for this subcluster, loaded lazily from the ruleset."""
        if self._rule is None:
            self._rule = get_ruleset().get_rule_by_name(self.rule_name)
        return self._rule

    @property
    def conditions_str(self) -> str:
        """Condition string for display, with outer parentheses stripped."""
        text = str(self.rule.conditions)
        if text.startswith("(") and text.endswith(")"):
            return text[1:-1]
        return text

    @property
    def domain_hits(self) -> list[DomainHit]:
        """Domain hits that fired for this subcluster, in deterministic order."""
        if self._domain_hits is None:
            self._enrich()
        return self._domain_hits

    @property
    def compound(self) -> CompoundInfo:
        """Compound metadata."""
        if not self._compound:
            self._compound = get_subcluster_compounds().get(self.rule_name)
        return self._compound

    @property
    def cds_names(self) -> list[str]:
        """Sorted unique CDS names that contributed to this prediction."""
        return sorted({cr.cds.get_name() for cr in self.cds_results})

    def _enrich(self) -> Self:
        """Compute domain hits from ``cds_results`` and store them in ``_domain_hits``.

        Returns ``self`` so callers can write ``hit._enrich()`` inline.

        Note: AA location (query_start/query_end) is not available here.
        HMMerHit carries it, but it is discarded when converted to SecMetQualifier.Domain
        in hmm_rule_parser/cluster_prediction.py. A fix there would be needed to expose it.
        """
        profiles = get_subcluster_profiles()
        hits: list[DomainHit] = []
        for cds_result in self.cds_results:
            fired = cds_result.definition_domains.get(self.rule_name, set())
            for domain_name in sorted(fired):
                matching = [d for d in cds_result.domains if d.name == domain_name]
                best = max(matching, key=lambda d: d.bitscore) if matching else None
                hits.append(DomainHit(
                    profile=profiles[domain_name],
                    cds_name=cds_result.cds.get_name(),
                    evalue=best.evalue if best else None,
                    bitscore=best.bitscore if best else None,
                ))
        self._domain_hits = hits
        return self

    def to_json(self) -> dict[str, Any]:
        """Serialise the raw detection data for this prediction."""
        return {
            "rule_name": self.rule_name,
            "start": self.start,
            "end": self.end,
            "cds_results": [cr.to_json() for cr in self.cds_results],
        }

    @classmethod
    def from_json(cls, data: dict[str, Any], record: Record) -> Self:
        """Reconstruct an un-enriched SubclusterPrediction from a serialised dict.

        The returned prediction must be enriched by calling
        ``enrich(rule, profiles, compounds)`` before its derived properties
        can be accessed.
        """
        cds_results = [CDSResults.from_json(cr, record) for cr in data["cds_results"]]
        return cls(
            rule_name=data["rule_name"],
            start=data["start"],
            end=data["end"],
            cds_results=cds_results,
        )

    def __repr__(self) -> str:
        return (
            f"SubclusterPrediction(rule_name={self.rule_name!r}, "
            f"core={self.start}-{self.end}, "
            f"cds_count={len(self.cds_results)})"
        )


class SubclusterDetectionResults(DetectionResults):
    """Results class for the Subcluster detection module """

    schema_version = "1"  # increment when the JSON format changes

    def __init__(
            self,
            record_id: str,
            rule_results: RuleDetectionResults,
            rule_names: set[str],
            strictness: str = "strict",
    ) -> None:
        super().__init__(record_id)
        self.rule_results = rule_results
        self.rule_names = rule_names
        self.strictness = strictness
        self.hits = [
            SubclusterPrediction(
                rule_name=protocluster.product,
                start=protocluster.core_location.start,
                end=protocluster.core_location.end,
                cds_results=cds_results,
            )
            for protocluster, cds_results in rule_results.cds_by_cluster.items()
        ]

    def get_hits_for_region(self, region: Region) -> list[SubclusterPrediction]:
        """Return all hits fully contained within the given region."""
        return [
            hit for hit in self.hits
            if location_contains_other(region.location, FeatureLocation(hit.start, hit.end))
        ]

    def get_hits_outside_regions(self, record: Record) -> list[SubclusterPrediction]:
        """Return hits not fully contained by any region in the record."""
        return [
            hit for hit in self.hits
            if not any(
                location_contains_other(region.location, FeatureLocation(hit.start, hit.end))
                for region in record.get_regions()
            )
        ]

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "strictness": self.strictness,
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
            logging.debug("Discarding subcluster results: RuleDetectionResults schema changed")
            return None
        return cls(
            record_id=data["record_id"],
            rule_results=rule_results,
            rule_names=set(data.get("rule_names", [])),
            strictness=data.get("strictness", "strict"),
        )

    def add_to_record(self, record: Record) -> None:
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        self.rule_results.annotate_cds_features()