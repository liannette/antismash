from dataclasses import dataclass
from typing import Optional

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSection, HTMLSections, Markup
from antismash.common.layers import RecordLayer, RegionLayer
from antismash.common.json import JSONCompatible
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Record, Region

from antismash.config import ConfigType

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
    
    #def to_json():
        
    
def will_handle(products: list[str], categories: set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def generate_html(region_layer: RegionLayer, results: Optional[DetectionResults],
                    record_layer: RecordLayer, options: ConfigType) -> HTMLSections: 
    # build and render template    
    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))

    section = template.render()
    html = HTMLSections(name="subclusters")
    html.add_detail_section("Subclusters", section, class_name="subclusters")
    return html


def generate_javascript_data(record: Record, region: Region,
                                results: DetectionResults) -> JSONCompatible:
    # return dome fake data
    return {}