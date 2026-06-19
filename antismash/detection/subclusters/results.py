import logging
import json
from dataclasses import dataclass, field
from  typing import Any, Optional, Self

from antismash.common.hmm_rule_parser.rule_parser import DetectionRule
from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults
from antismash.common.secmet import Record
from antismash.common.secmet.features import SubRegion
from antismash.common.secmet.locations import FeatureLocation
from antismash.common.module_results import DetectionResults


@dataclass(frozen=True)
class RuleInfo:
    """Data for one subcluster detection rule."""
    identifier: str
    description: str
    conditions: str


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


@dataclass(frozen=True)
class SubclusterHitHtmlData:
    """Data for one subcluster hit, as needed by the HTML generator.
    
    Attributes:         
        rule: RuleInfo instance with metadata about the rule.
        domain_hits: List of ``DomainHit`` instances for every domain hit that contributed to this subcluster hit.
        compound: CompoundInfo instance with metadata about the predicted compound.
    """
    rule: RuleInfo
    domain_hits: list[DomainHit]
    compound: CompoundInfo
    
    @property
    def cds_locus_tags(self) -> list[str]:
        """Sorted list of all CDS locus tags that contributed to this hit."""
        return sorted({hit.cds_locus_tag for hit in self.domain_hits})
    

@dataclass
class SubclusterHit:
    """A single detected subcluster.
 
    Attributes:
        rule: The ``DetectionRule`` that fired.
        start: Start of the core location (0-based bp).
        end: End of the core location (0-based bp).
        cds_results: ``CDSResults`` instances for every CDS that contributed 
            to this hit, as returned by the rule-based detection pipeline.
        metadata: Optional ``SubclusterMetadata`` with additional metadata about the hit.
    """
    rule: DetectionRule
    start: int
    end: int
    cds_results: list[CDSResults]
    metadata: Optional[Any] = None  # Placeholder for SubclusterMetadata


class SubclusterDetectionResults(DetectionResults): 
    """Results class for the Subcluster detection module """
    schema_version = 1 # increment when the data format in the results changes

    def __init__(self, record_id: str, rules_version: str, 
                 hits_by_region: dict[int, list[SubclusterHit]]) -> None:
        super().__init__(record_id)
        self.rules_version = rules_version
        self.hits_by_region = hits_by_region

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "rules_version": self.rules_version,
            "hits_by_region": {
                region_number: [hit.to_dict() for hit in hits]
                for region_number, hits in self.hits_by_region.items()
            }
        }

    @classmethod
    def from_json(cls, data: dict[str, Any], record: Record) -> Self:
        # check that the previous data version is the same as current, if not, discard the results
        if data["schema_version"] != SubclusterDetectionResults.schema_version:
            return None
        
        if data["rules_version"] != SubclusterDetectionResults.schema_version:
            logging.warning("Rules version in previous results does not match current version.")

        results = cls(record_id=data["record_id"], rules_version=data["rules_version"], hits_by_region={})
        return results

    def add_to_record(self, record):
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        # any results would be added here
        # for an example of new features, see antismash.modules.tta
        # for an example of qualifiers, see antismash.modules.t2pks
        # any new feature types or qualifiers would be implemented in antismash.common.secmet,
        #   and would need to be able to be converted to and from biopython's SeqFeature without loss
        raise NotImplementedError()  # remove this when completed

# @dataclass(frozen=True)
# class SubclusterCompound:
#     """Compound information for a subcluster.
 
#     Attributes:
#         name: Display name of the compound.
#         smiles: SMILES string for the compound, if available.
#         types: One or more pathway/class labels (e.g. ``["amino acid"]``).
#     """
#     name: str
#     smiles: Optional[str]
#     types: list[str] = field(default_factory=list)

#     @classmethod
#     def from_dict(cls, data: dict) -> "SubclusterCompound":
#         return cls(
#             name=data.get("name"),
#             smiles=data.get("smiles"),
#             type=data.get("type", []),
#         )

# class SubclusterHit:
#     """A single detected subcluster.
 
#     Attributes:
#         rule: The ``DetectionRule`` that fired.
#         core_start: Start of the core location (0-based bp).
#         core_end: End of the core location (0-based bp).
#         cds_results: ``CDSResults`` instances for every CDS that contributed 
#             to this hit, as returned by the rule-based detection pipeline.
#         compound: ``SubclusterCompound`` with metadata about the predicted compound.
#         tool: Written into the ``SubRegion.tool`` qualifier so the feature can
#             be identified in GenBank output and by ``Region.has_subregion_by_tool()``.
#     """
 
#     TOOL = "subcluster-rule-based"
 
#     def __init__(
#         self,
#         rule: DetectionRule,
#         start: int,
#         end: int,
#         cds_results: list[CDSResults],
#         metadata: Optional[SubclusterMetadata] = None,
#     ) -> None:
#         self.rule = rule
#         self.start = start
#         self.end = end
#         self.cds_results = cds_results
#         self.metadata = metadata
 
#     @property
#     def name(self) -> str:
#         """Rule name, e.g. ``"SCG0042"``."""
#         return self.rule.name
    
#     @property
#     def conditions(self) -> str:
#         """String representation of the rule's CONDITIONS block."""
#         condition_text = str(self.rule.conditions)
#         # strip off outer parens if they exist
#         if condition_text[0] == "(" and condition_text[-1] == ')':
#             condition_text = condition_text[1:-1]
#         return condition_text
 
#     @property
#     def all_definition_domains(self) -> set[str]:
#         definition_domains: dict[str, set[str]] = dict()
#         for cds_result in self.cds_results:
#             cds_name = cds_result.get_name()
#             for domain_set in cds_result.definition_domains[self.name]:
#                 for domain in domain_set:
#                     definition_domains[domain].add(cds_name)
#         return definition_domains
 
#     @property
#     def compound_name(self) -> Optional[str]:
#         """Shortcut to the compound name"""
#         return self.metadata.compound.name
#         # Alternatively, could return the rule description, which is currently the same as the compound name:
#         # """Description from the rule file, e.g. the name of the subcluster compound."""
#         # return self.rule.description
    
#     @property
#     def compound_type(self) -> Optional[List[str]]:
#         """Shortcut to the compound type list"""
#         return self.metadata.compound.type
 
#     @property
#     def compound_smiles(self) -> Optional[str]:
#         """Shortcut to the SMILES string for the compound"""
#         return self.metadata.compound.smiles
 
#     def to_subregion(self) -> SubRegion:
#         """Convert to a ``SubRegion`` feature for adding to the antiSMASH record.
 
#         The rule name is stored as ``SubRegion.label`` so that results can be
#         looked up by iterating ``region.subregions`` and matching on ``.label``.
 
#         Returns:
#             A ``SubRegion`` ready for ``record.add_subregion()``.
#         """
#         location = FeatureLocation(self.start, self.end, strand=1)
#         return SubRegion(location, tool=self.TOOL, label=self.name)


#     def to_html_data(self, ruleset: Ruleset) -> SubclusterHitHtmlData:
#         """Build a self-contained HTML data object from this hit.

#         Arguments:
#             ruleset: The ``Ruleset`` used during detection, needed to look up
#                 domain descriptions and accessions from ``all_profiles``.

#         Returns:
#             A ``SubclusterHitHtmlData`` instance ready for the HTML generator.
#         """
#         rule_info = RuleInfo(
#             name=self.rule.name,
#             description=self.rule.description,
#             conditions=self.conditions,  # reuse existing property, strips outer parens
#         )

#         compound_info = None
#         if self.metadata and self.metadata.compound:
#             compound_info = CompoundInfo(
#                 name=self.metadata.compound.name,
#                 smiles=self.metadata.compound.smiles,
#                 classification=list(self.metadata.compound.type),
#             )

#         domain_hits = []
#         for cds_result in self.cds_results:
#             definition_domains = cds_result.definition_domains.get(self.name, set())
#             locus_tag = cds_result.cds.get_name()

#             # build a name->score lookup for this CDS from the SecMetQualifier domains
#             scores = {d.name: d for d in cds_result.domains}

#             for domain_name in sorted(definition_domains):
#                 profile = ruleset.all_profiles.get(domain_name)
#                 hmm_file = getattr(profile, "hmm_file", "")
#                 domain_info = DomainInfo(
#                     name=domain_name,
#                     acc=_read_hmm_accession(hmm_file),
#                     description=profile.description if profile else None,
#                 )
#                 domain_hits.append(DomainHit(
#                     cds_locus_tag=locus_tag,
#                     domain=domain_info,
#                 ))

#         return SubclusterHitHtmlData(
#             rule=rule_info,
#             domain_hits=domain_hits,
#             compound=compound_info,
#         )
    
#     def __repr__(self) -> str:
#         return (
#             f"SubclusterHit(rule={repr(self.name)}, "
#             f"core={self.start}-{self.end}, "
#             f"cds_count={len(self.cds_results)})"
#         )
 

# class SubclusterDetectionResults(DetectionResults): 
#     def __init__(self, record_id: str) -> None:
#         super().__init__(record_id)
#         self.hits_by_region: Dict[int, List[SubclusterHit]] = {}
 
    # # ------------------------------------------------------------------
    # # Accessors for HTML generator
    # # ------------------------------------------------------------------
 
    # def get_hits_for_region(self, region_number: int) -> List[SubclusterHit]:
    #     """Return all hits within the given region.
 
    #     Returns an empty list when the region has no subcluster hits, so
    #     callers do not need to guard against ``KeyError``.
 
    #     Arguments:
    #         region_number: Integer region number from
    #             ``Region.get_region_number()`` / ``RegionLayer.get_region_number()``.
    #     """
    #     return self.hits_by_region.get(region_number, [])
 
    # def get_all_hits(self) -> List[SubclusterHit]:
    #     """Flat list of every ``SubclusterHit`` across all regions."""
    #     hits: List[SubclusterHit] = []
    #     for region_hits in self.hits_by_region.values():
    #         hits.extend(region_hits)
    #     return hits
 
    # # ------------------------------------------------------------------
    # # ModuleResults interface
    # # ------------------------------------------------------------------
    
    # def to_json(self) -> dict:
    #     raise NotImplementedError(
    #         "JSON serialisation for SubclusterDetectionResults is not yet implemented."
    #     )
 
    # @staticmethod
    # def from_json(data: dict, record: Record) -> Optional["SubclusterDetectionResults"]:
    #     raise NotImplementedError(
    #         "JSON deserialisation for SubclusterDetectionResults is not yet implemented."
    #     )
 
    # def add_to_record(self, record: Record) -> None:
    #     """Add detected subclusters as ``SubRegion`` features to the record.
 
    #     Skips silently if features from this tool are already present.
 
    #     Arguments:
    #         record: The ``Record`` to annotate.
    #     """
    #     existing_tools = {sr.tool for sr in record.get_subregions()}
    #     if SubclusterHit.TOOL in existing_tools:
    #         return
    #     for hit in self.get_all_hits():
    #         record.add_subregion(hit.to_subregion())
 
    # def get_predicted_subregions(self) -> List[SubRegion]:
    #     """Return all ``SubRegion`` features predicted by this module.
 
    #     Called by the antiSMASH main loop to collect features before adding
    #     them to the record via ``add_to_record()``.
    #     """
    #     return [hit.to_subregion() for hit in self.get_all_hits()]
 
