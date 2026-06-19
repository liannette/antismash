from typing import Optional

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSection, HTMLSections, Markup
from antismash.common.layers import RecordLayer, RegionLayer
from antismash.common.json import JSONBase
from antismash.common.secmet import Record, Region

from antismash.config import ConfigType

from .results import SubclusterDetectionResults
        
    
def will_handle(products: list[str], categories: set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def generate_html(region_layer: RegionLayer, results: Optional[SubclusterDetectionResults],
                    record_layer: RecordLayer, options: ConfigType) -> HTMLSections: 
    # build and render template    
    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))

    section = template.render()
    html = HTMLSections(name="subclusters")
    html.add_detail_section("Subclusters", section, class_name="subclusters")
    return html


def generate_javascript_data(record: Record, region: Region,
                                results: SubclusterDetectionResults) -> JSONBase:
    # return some fake data
    return {}