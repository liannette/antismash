"""Manages HTML construction for the subcluster detection module."""

from types import SimpleNamespace
from typing import Optional

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, Markup
from antismash.common.layers import RecordLayer, RegionLayer
from antismash.common.json import JSONBase
from antismash.common.secmet import Record, Region
from antismash.config import ConfigType

from .results import CompoundInfo, HmmHit, SubclusterDetectionResults, SubclusterPrediction
from .signatures import SubclusterHmmSignature


def will_handle(products: list[str], categories: set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def _get_fake_hits() -> list[SubclusterPrediction]:
    """Hard-coded enriched predictions for template development.

    Remove this function and its call-site in ``generate_html`` once real
    detection data flows through ``SubclusterDetectionResults``.
    """
    def _fake_profile(name: str, accession: Optional[str],
                      description: Optional[str]) -> SubclusterHmmSignature:
        return SubclusterHmmSignature(name=name, description=description or "",
                                      cutoff=0, hmm_path="", accession=accession)

    hit_a = SubclusterPrediction(rule_name="SCG0042", start=27326, end=80190, cds_results=[])
    hit_a._rule = SimpleNamespace(conditions="cds(PDH_N and PDH_C) and Aminotran_1_2 and Glyoxalase and FMN_dh")
    hit_a._domain_hits = [
        HmmHit(profile=_fake_profile("FMH_dh", "PF01070", "FMN-dependent dehydrogenase"), cds_locus_tag="AJAP_32035"),
        HmmHit(profile=_fake_profile("FMH_dh", "PF01070", "FMN-dependent dehydrogenase"), cds_locus_tag="AJAP_32036"),
        HmmHit(profile=_fake_profile("Glyoxylase", "PF00903", "Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily"), cds_locus_tag="AJAP_32040"),
        HmmHit(profile=_fake_profile("Aminotran_1_2", "PF00155", "Aminotransferase class I and II"), cds_locus_tag="AJAP 32060"),
        HmmHit(profile=_fake_profile("PDH_N", "PF02153", "Prephenate dehydrogenase, nucleotide-binding domain"), cds_locus_tag="AJAP_32155"),
        HmmHit(profile=_fake_profile("PDH_C", "PF00903", "Prephenate dehydrogenase, dimerization domain"), cds_locus_tag="AJAP_32155"),
    ]
    hit_a._compound = CompoundInfo(
        name="4-Hydroxyphenylglycine (Hpg)",
        smiles="C1=CC(=CC=C1C(C(=O)O)N)O",
        classification=["amino acid", "precursor"],
    )
    hit_a._enriched = True

    return [hit_a,]


def generate_html(region_layer: RegionLayer, results: Optional[SubclusterDetectionResults],
                  record_layer: RecordLayer, options: ConfigType) -> HTMLSections:
    """Build the detail-panel HTML for subcluster hits in this region."""
    # if results is not None:
    #     hits = results.get_hits_for_region(region_layer.region_feature)
    # else:
    #     hits = []
    hits = _get_fake_hits()

    tooltip = Markup("Subclusters are sets of genes responsible for producing a specific chemical moiety.")

    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    section = template.render(hits=hits, tooltip=tooltip)

    html = HTMLSections(name="subclusters")
    html.add_detail_section("Subclusters", section, class_name="subclusters")
    return html


def generate_javascript_data(record: Record, region: Region,
                                results: SubclusterDetectionResults) -> JSONBase:
    # return some fake data
    return {}