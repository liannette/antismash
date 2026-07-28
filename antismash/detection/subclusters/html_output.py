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

from .results import SubclusterDetectionResults, SubclusterPrediction


def will_handle(products: list[str], categories: set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def generate_html(region_layer: RegionLayer, results: Optional[SubclusterDetectionResults],
                  record_layer: RecordLayer, options: ConfigType) -> HTMLSections:
    """Build the detail-panel HTML for subcluster predictions in this region."""
    predictions = results.get_predictions_for_region(region_layer.region_feature)

    tooltip = Markup("Subclusters are sets of genes responsible for producing a specific chemical moiety.")

    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    section = template.render(predictions=predictions, tooltip=tooltip, anchor=region_layer.anchor_id)

    html = HTMLSections(name="subclusters")
    html.add_detail_section("Subclusters", section, class_name="subclusters")
    return html


def generate_javascript_data(record: Record, region: Region,
                             results: SubclusterDetectionResults) -> JSONBase:
    region_anchor = f"r{record.record_index}c{region.get_region_number()}"

    predictions = results.get_predictions_for_region(region)

    data = []
    for i, prediction in enumerate(predictions, start=1):
        data.append({
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

    return data
