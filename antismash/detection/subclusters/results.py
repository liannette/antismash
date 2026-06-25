import logging
from dataclasses import dataclass, field
from  typing import Any, Optional, Self

from antismash.common.hmm_rule_parser.rule_parser import DetectionRule
from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Record, Region
from antismash.common.secmet.locations import FeatureLocation

from .signatures import SubclusterHmmSignature


@dataclass(frozen=True)
class CompoundInfo:
    """Compound information for a subcluster."""
    name: str
    smiles: Optional[str]
    classification: list[str]


@dataclass(frozen=True)
class HmmHit:
    """A single profile match that contributed to a subcluster hit."""
    profile: SubclusterHmmSignature
    cds_name: str


@dataclass
class SubclusterPrediction:
    """A single detected subcluster.

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

    # Derived fields — populated by enrich(); not part of __init__
    _rule: Optional[DetectionRule] = field(default=None, init=False, repr=False)
    _conditions_str: str = field(default="", init=False, repr=False)
    _domain_hits: list[HmmHit] = field(default_factory=list, init=False, repr=False)
    _compound: CompoundInfo = field(default=None, init=False, repr=False)
    _enriched: bool = field(default=False, init=False, repr=False)

    @property
    def rule(self) -> DetectionRule:
        """The ``DetectionRule`` that fired, available after enrichment."""
        self._require_enriched()
        return self._rule

    @property
    def conditions_str(self) -> str:
        """Condition string for display, with outer parentheses stripped."""
        self._require_enriched()
        text = str(self.rule.conditions)
        if text.startswith("(") and text.endswith(")"):
            return text[1:-1]
        return text

    @property
    def domain_hits(self) -> list[HmmHit]:
        """Domain hits that fired for this subcluster, in deterministic order."""
        self._require_enriched()
        return self._domain_hits

    @property
    def compound(self) -> CompoundInfo:
        """Compound metadata, or ``None`` when metadata is unavailable."""
        self._require_enriched()
        return self._compound

    @property
    def cds_locus_tags(self) -> list[str]:
        """Sorted unique CDS locus tags that contributed to this prediction."""
        return sorted({cr.cds.get_name() for cr in self.cds_results})

    def enrich(self, rule: DetectionRule,
               profiles: dict[str, SubclusterHmmSignature],
               compound: CompoundInfo) -> Self:
        """Populate derived fields from the fired rule and pre-loaded metadata.

        Arguments:
            rule: The ``DetectionRule`` whose name matches ``self.rule_name``.
            profiles: Mapping of profile name to ``SubclusterHmmSignature``,
                as returned by ``get_subcluster_profiles()``.
            compound: Compound metadata for this rule, or ``None``.

        Returns ``self`` so callers can write ``prediction.enrich(...)`` inline.
        """
        self._rule = rule

        hits: list[HmmHit] = []
        for cds_result in self.cds_results:
            for domain_name in sorted(cds_result.definition_domains.get(self.rule_name, set())):
                hits.append(HmmHit(
                    profile=profiles[domain_name],
                    cds_name=cds_result.cds.get_name(),
                ))
        self._domain_hits = hits
        self._compound = compound
        self._enriched = True
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

    def _require_enriched(self) -> None:
        if not self._enriched:
            raise RuntimeError(
                f"{self!r} has not been enriched — call enrich(rule, profiles, compound) first"
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
            rules_version: str,
            hits: list[SubclusterPrediction],
    ) -> None:
        super().__init__(record_id)
        self.rules_version = rules_version
        self.hits = hits

    def get_hits_for_region(self, region: Region) -> list[SubclusterPrediction]:
        """Return all hits that overlap the given region."""
        return [
            hit for hit in self.hits
            if region.overlaps_with(FeatureLocation(hit.start, hit.end))
        ]

    def get_hits_outside_regions(self, record: Record) -> list[SubclusterPrediction]:
        """Return hits that do not overlap any region in the record."""
        return [
            hit for hit in self.hits
            if not any(
                region.overlaps_with(FeatureLocation(hit.start, hit.end))
                for region in record.get_regions()
            )
        ]

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "rules_version": self.rules_version,
            "hits": [hit.to_json() for hit in self.hits],
        }

    @classmethod
    def from_json(cls, data: dict[str, Any], record: Record) -> Optional[Self]:
        if data.get("schema_version") != cls.schema_version:
            logging.debug(
                "Discarding subcluster results: schema version %s != %s",
                data.get("schema_version"), cls.schema_version,
            )
            return None
        return cls(
            record_id=data["record_id"],
            rules_version=data["rules_version"],
            hits=[SubclusterPrediction.from_json(hit_data, record) for hit_data in data["hits"]],
        )

    def add_to_record(self, record: Record) -> None:
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        # any results would be added here
        # for an example of new features, see antismash.modules.tta
        # for an example of qualifiers, see antismash.modules.t2pks
        # any new feature types or qualifiers would be implemented in antismash.common.secmet,
        #   and would need to be able to be converted to and from biopython's SeqFeature without loss
        raise NotImplementedError()  # remove this when completed