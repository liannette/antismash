"""Manages HTML construction for the subcluster detection module."""

from types import SimpleNamespace
from typing import Optional

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, Markup
from antismash.common.layers import RecordLayer, RegionLayer
from antismash.common.json import JSONBase
from antismash.common.secmet import Record, Region
from antismash.common.secmet.locations import FeatureLocation
from antismash.config import ConfigType

from .compounds import CompoundInfo
from .results import SubclusterDetectionResults, SubclusterPrediction


def will_handle(products: list[str], categories: set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def _get_fake_hits() -> list[SubclusterPrediction]:
    """Hard-coded enriched predictions for template development.

    Remove this function and its call-site in ``generate_html`` once real
    detection data flows through ``SubclusterDetectionResults``.
    """
    def _fake_cds_results(rule_name: str,
                          cds_domains: dict[str, list[tuple[str, float, float]]]) -> list:
        results = []
        for cds_name, domains in cds_domains.items():
            results.append(SimpleNamespace(
                cds=SimpleNamespace(get_name=lambda n=cds_name: n),
                domains=[
                    SimpleNamespace(name=name, evalue=evalue, bitscore=bitscore)
                    for name, evalue, bitscore in domains
                ],
                definition_domains={rule_name: {name for name, _, _ in domains}},
            ))
        return results

    hit_a = SubclusterPrediction(
        rule_name="SCG0006",
        core_location=FeatureLocation(17346, 21101),
        cds_results=_fake_cds_results("SCG0006", {
            "AJAP_31990": [("ECH_1", 1.2e-18, 65.3)],
            "AJAP_31995": [("ECH_1", 3.4e-21, 72.1)],
            "AJAP_32000": [("ECH_1", 8.7e-20, 68.9)],
            "AJAP_32005": [
                ("Chal_sti_synt_N", 2.1e-45, 152.4),
                ("Chal_sti_synt_C", 5.8e-34, 118.1),
            ],
        }),
        rule=SimpleNamespace(
            conditions="ECH_1 and cds(Chal_sti_synt_C and Chal_sti_synt_N)",
            description="3,5-Dihydroxyphenylglycine (Dhpg)",
        ),
        compound=CompoundInfo(
            name="3,5-Dihydroxyphenylglycine (Dhpg)",
            smiles="C1=C(O)C=C(O)C=C1[C@@H](C(=O)O)N",
            classification=["amino acid", "precursor"],
        ),
    )

    hit_b = SubclusterPrediction(
        rule_name="SCG0042",
        core_location=FeatureLocation(27326, 80190),
        cds_results=_fake_cds_results("SCG0042", {
            "AJAP_32035": [("FMN_dh", 4.1e-29, 98.7)],
            "AJAP_32040": [("Glyoxalase_4", 7.3e-15, 54.2)],
            "AJAP_32060": [("Aminotran_1_2", 9.6e-38, 128.5)],
            "AJAP_32155": [
                ("PDH_N", 1.4e-22, 78.3),
                ("PDH_C", 6.2e-17, 60.1),
            ],
        }),
        rule=SimpleNamespace(
            conditions="cds(PDH_N and PDH_C) and Aminotran_1_2 and Glyoxalase_4 and FMN_dh",
            description="4-Hydroxyphenylglycine (Hpg)",
        ),
        compound=CompoundInfo(
            name="4-Hydroxyphenylglycine (Hpg)",
            smiles="C1=CC(=CC=C1C(C(=O)O)N)O",
            classification=["amino acid", "precursor"],
        ),
    )

    return [hit_a, hit_b]


def generate_html(region_layer: RegionLayer, results: Optional[SubclusterDetectionResults],
                  record_layer: RecordLayer, options: ConfigType) -> HTMLSections:
    """Build the detail-panel HTML for subcluster hits in this region."""
    # if results is not None:
    #     hits = results.get_hits_for_region(region_layer.region_feature)
    # else:
    #     hits = []

    hits = _get_fake_hits()
    enum_hits = list(enumerate(hits, start=1))

    tooltip = Markup("Subclusters are sets of genes responsible for producing a specific chemical moiety.")

    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    section = template.render(enum_hits=enum_hits, tooltip=tooltip, anchor=region_layer.anchor_id)

    html = HTMLSections(name="subclusters")
    html.add_detail_section("Subclusters", section, class_name="subclusters")
    return html


def generate_javascript_data(record: Record, region: Region,
                             results: SubclusterDetectionResults) -> JSONBase:
    region_anchor = f"r{record.record_index}c{region.get_region_number()}"

    predictions = _get_fake_hits()

    javascript_data = []
    for i, prediction in enumerate(predictions, start=1):
        javascript_data.append({
            "identifier": f"subclusters-svg-{region_anchor}-sc{i}",
            "cds_results": [
                {
                    "cds": cds_name,
                    "domains": [
                        {
                            "name": hit.domain_name,
                            "description": hit.domain_description,
                            "accession": hit.domain_accession,
                            "evalue": hit.evalue,
                            "bitscore": hit.bitscore,
                        }
                        for hit in hits
                    ],
                }
                for cds_name, hits in prediction.domain_hits_by_cds.items()
            ],
        })

    return javascript_data
