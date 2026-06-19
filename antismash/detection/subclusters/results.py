import logging
from dataclasses import dataclass, field
from  typing import Any, Optional, Self

from antismash.common.hmm_rule_parser.rule_parser import DetectionRule
from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults, Ruleset
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Record, Region
from antismash.common.secmet.locations import FeatureLocation



@dataclass(frozen=True)
class CompoundInfo:
    """Compound information for a subcluster."""
    name: str
    smiles: Optional[str]
    classification: list[str]


@dataclass(frozen=True)
class DomainInfo:
    """Data for one phmm domain"""
    name: str
    acc: Optional[str]
    description: Optional[str]


@dataclass(frozen=True)
class DomainHit:
    """Data for one phmm domain hit that contributed to a subcluster hit."""
    cds_locus_tag: str
    domain: DomainInfo


def _read_hmm_accession(hmm_file: str) -> Optional[str]:
    """Parse the bare Pfam accession from an HMM file header.
 
    Pfam accessions are stored in the ``ACC`` field as ``PFxxxxx.N``
    (name + version suffix).  antiSMASH's Python objects do not expose this
    value, so it must be read directly from the file.
 
    Returns ``None`` for custom (non-Pfam) profiles or when ``hmm_file`` is
    empty / unreadable.
    """
    if not hmm_file:
        return None
    try:
        with open(hmm_file, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("ACC"):
                    # "ACC   PF01278.15\n"  →  "PF01278"
                    return line.split()[1].split(".")[0]
                if line.startswith("HMM "):
                    break  # reached model block — no ACC line present
    except OSError:
        logging.warning("Could not read HMM file for accession lookup: %s", hmm_file)
    return None


@dataclass
class SubclusterHit:
    """A single detected subcluster.

    After detection, call ``enrich(ruleset)`` once (while the Ruleset is still
    in scope) to populate the derived fields used by the HTML template.

 
    Attributes:
        rule: The ``DetectionRule`` that fired.
        start: Start of the core location (0-based bp).
        end: End of the core location (0-based bp).
        cds_results: ``CDSResults`` instances for every CDS that contributed 
            to this hit, as returned by the rule-based detection pipeline.
    """
    rule: DetectionRule
    start: int
    end: int
    cds_results: list[CDSResults]

    # Derived fields — populated by enrich(); not part of __init__
    _conditions_str: str = field(default="", init=False, repr=False)
    _domain_hits: list[DomainHit] = field(default_factory=list, init=False, repr=False)
    _compound: Optional[CompoundInfo] = field(default=None, init=False, repr=False)
    _enriched: bool = field(default=False, init=False, repr=False)

    @property
    def conditions_str(self) -> str:
        """Condition string for display, with outer parentheses stripped."""
        text = str(self.rule.conditions)
        if text.startswith("(") and text.endswith(")"):
            text = text[1:-1]
        return text

    @property
    def domain_hits(self) -> list[DomainHit]:
        """Domain hits that fired for this subcluster, in deterministic order."""
        self._require_enriched()
        return self._domain_hits
 
    @property
    def compound(self) -> Optional[CompoundInfo]:
        """Compound metadata, or ``None`` when metadata is unavailable."""
        self._require_enriched()
        return self._compound
 
    # @property
    # def cds_locus_tags(self) -> list[str]:
    #     """Sorted list of unique CDS locus tags that contributed to this hit."""
    #     return sorted({cds_result.cds.get_name() for cds_result in self.cds_results})

    def enrich(self, ruleset: Ruleset, compounds: dict[str, CompoundInfo]) -> "SubclusterHit":
        """Populate derived fields from the detection ruleset.
 
        Must be called while the ``Ruleset`` is still in scope (i.e. inside
        ``run_on_record``), because domain descriptions and Pfam accessions are
        looked up from ``ruleset.all_profiles``.
 
        Returns ``self`` so callers can write ``hit.enrich(ruleset)`` inline.
        """
        # conditions string
        self._conditions_str = self.conditions_str
 
        # domain hits
        hits: list[DomainHit] = []
        for cds_result in self.cds_results:
            locus_tag = cds_result.cds.get_name()
            fired = sorted(cds_result.definition_domains.get(self.rule.name, set()))
            for domain_name in fired:
                profile = ruleset.all_profiles.get(domain_name)
                hits.append(DomainHit(
                    cds_locus_tag=locus_tag,
                    domain=DomainInfo(
                        name=domain_name,
                        acc=_read_hmm_accession(getattr(profile, "hmm_file", "")),
                        description=profile.description if profile else None,
                    ),
                ))
        self._domain_hits = hits
 
        # compound info
        compound: CompoundInfo = compounds.get(self.rule.name)
        self._compound = compound
 
        self._enriched = True
        return self
 
    def _require_enriched(self) -> None:
        if not self._enriched:
            raise RuntimeError(
                f"{self!r} has not been enriched — call enrich(ruleset) first"
            )
 
    def __repr__(self) -> str:
        return (
            f"SubclusterHit(rule={repr(self.rule.name)}, "
            f"core={self.start}-{self.end}, "
            f"cds_count={len(self.cds_results)})"
        )


class SubclusterDetectionResults(DetectionResults): 
    """Results class for the Subcluster detection module """

    schema_version = 1 # increment when the data format in the results changes

    def __init__(
            self, record_id: str, 
            rules_version: str, 
            hits: list[SubclusterHit]
        ) -> None:
        super().__init__(record_id)
        self.rules_version = rules_version
        self.hits = hits

    def get_hits_for_region(self, region: Region) -> list[SubclusterHit]:
        """Return all hits that overlap the given region."""
        return [hit for hit in self.hits if region.overlaps_with(FeatureLocation(hit.start, hit.end))]
 
    def get_hits_outside_regions(self, record: Record) -> list[SubclusterHit]:
        """Return hits that do not overlap any region in the record."""
        return [
            hit for hit in self.hits
            if not any(region.overlaps_with(FeatureLocation(hit.start, hit.end)) for region in record.get_regions())
        ]

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "rules_version": self.rules_version,
            "hits": [hit.to_dict() for hit in self.hits],
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        # check that the previous data version is the same as current, if not, discard the results
        if data["schema_version"] != SubclusterDetectionResults.schema_version:
            return None
        if data["rules_version"] != SubclusterDetectionResults.schema_version:
            logging.warning("Rules version in previous results does not match current version.")

        return cls(
            record_id=data["record_id"], 
            rules_version=data["rules_version"], 
            hits=data["hits"]
            )

    def add_to_record(self, record):
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        # any results would be added here
        # for an example of new features, see antismash.modules.tta
        # for an example of qualifiers, see antismash.modules.t2pks
        # any new feature types or qualifiers would be implemented in antismash.common.secmet,
        #   and would need to be able to be converted to and from biopython's SeqFeature without loss
        raise NotImplementedError()  # remove this when completed