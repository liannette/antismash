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

    def _fake_cds_results(names: list[str]) -> list:
        return [SimpleNamespace(cds=SimpleNamespace(get_name=lambda n=tag: n)) for tag in names]

    hit_a = SubclusterPrediction(rule_name="SCG0006", start=17346, end=21101,
                                 cds_results=_fake_cds_results([
                                     "AJAP_31990", "AJAP_31995", "AJAP_32000",
                                     "AJAP_32005",
                                 ]))
    hit_a._rule = SimpleNamespace(
        conditions="ECH_1 and cds(Chal_sti_synt_C and Chal_sti_synt_N)",
        description="3,5-Dihydroxyphenylglycine (Dhpg)",
    )
    hit_a._domain_hits = [
        HmmHit(profile=_fake_profile("ECH_1", "PF00378", "Enoyl-CoA hydratase/isomerase"), cds_locus_tag="AJAP_31990"),
        HmmHit(profile=_fake_profile("ECH_1", "PF00378", "Enoyl-CoA hydratase/isomerase"), cds_locus_tag="AJAP_31995"),
        HmmHit(profile=_fake_profile("ECH_1", "PF00378", "Enoyl-CoA hydratase/isomerase"), cds_locus_tag="AJAP_32000"),
        HmmHit(profile=_fake_profile("Chal_sti_synt_N", "PF00195", "Chalcone and stilbene synthases, N-terminal domain"), cds_locus_tag="AJAP_32005"),
        HmmHit(profile=_fake_profile("Chal_sti_synt_C", "PF02797", "Chalcone and stilbene synthases, C-terminal domain"), cds_locus_tag="AJAP_32005"),
    ]
    hit_a._compound = CompoundInfo(
        name="3,5-Dihydroxyphenylglycine (Dhpg)",
        smiles="C1=C(O)C=C(O)C=C1[C@@H](C(=O)O)N",
        classification=["amino acid", "precursor"],
    )
    hit_a._enriched = True

    hit_b = SubclusterPrediction(rule_name="SCG0042", start=27326, end=80190,
                                 cds_results=_fake_cds_results([
                                     "AJAP_32035", "AJAP_32036", "AJAP_32040",
                                     "AJAP_32060", "AJAP_32155",
                                 ]))
    hit_b._rule = SimpleNamespace(
        conditions="cds(PDH_N and PDH_C) and Aminotran_1_2 and Glyoxalase and FMN_dh",
        description="4-Hydroxyphenylglycine (Hpg)",
    )
    hit_b._domain_hits = [
        HmmHit(profile=_fake_profile("FMH_dh", "PF01070", "FMN-dependent dehydrogenase"), cds_locus_tag="AJAP_32035"),
        HmmHit(profile=_fake_profile("Glyoxylase", "PF00903", "Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily"), cds_locus_tag="AJAP_32040"),
        HmmHit(profile=_fake_profile("Aminotran_1_2", "PF00155", "Aminotransferase class I and II"), cds_locus_tag="AJAP_32060"),
        HmmHit(profile=_fake_profile("PDH_N", "PF02153", "Prephenate dehydrogenase, nucleotide-binding domain"), cds_locus_tag="AJAP_32155"),
        HmmHit(profile=_fake_profile("PDH_C", "PF00903", "Prephenate dehydrogenase, dimerization domain"), cds_locus_tag="AJAP_32155"),
    ]
    hit_b._compound = CompoundInfo(
        name="4-Hydroxyphenylglycine (Hpg)",
        smiles="C1=CC(=CC=C1C(C(=O)O)N)O",
        classification=["amino acid", "precursor"],
    )
    hit_b._enriched = True

    return [hit_a, hit_b]


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