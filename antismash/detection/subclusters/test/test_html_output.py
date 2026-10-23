# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import re
from types import SimpleNamespace
import unittest

from antismash.common.layers import RegionLayer
from antismash.common.secmet.locations import CompoundLocation, FeatureLocation
from antismash.detection.subclusters import html_output
from antismash.detection.subclusters.compounds import CompoundInfo
from antismash.detection.subclusters.predictions import CDSDomainHit, SubclusterPrediction


def build_hit(cds_name, domain_name, description="desc"):
    return CDSDomainHit(domain_name=domain_name, domain_description=description,
                        domain_accession=f"{domain_name}.1", cds_name=cds_name,
                        evalue=1e-10, bitscore=50.)


def build_prediction(location, hits=None, rule_name="SCG0001"):
    return SubclusterPrediction(
        rule=SimpleNamespace(name=rule_name, conditions="(A and B)"),
        location=location,
        compound=CompoundInfo(name=f"{rule_name} compound", smiles="CCO",
                              classification=["amino acid"]),
        domain_hits=hits or [build_hit("cds1", "A")],
    )


def build_results(predictions):
    """ Builds results returning the given predictions for any region """
    return SimpleNamespace(get_predictions_for_region=lambda region: predictions)


class TestWillHandle(unittest.TestCase):
    def test_always(self):
        assert html_output.will_handle([], set())
        assert html_output.will_handle(["NRPS"], {"NRPS"})


class TestJavascriptData(unittest.TestCase):
    def setUp(self):
        self.record = SimpleNamespace(record_index=2)
        self.region = SimpleNamespace(location=FeatureLocation(100, 500, 1),
                                      get_region_number=lambda: 3)

    def generate(self, predictions):
        return html_output.generate_javascript_data(self.record, self.region,
                                                    build_results(predictions))

    def test_empty(self):
        assert self.generate([]) == []

    def test_identifiers(self):
        data = self.generate([build_prediction(FeatureLocation(150, 250, 1)),
                              build_prediction(FeatureLocation(300, 400, 1))])
        # numbered from one, matching the anchors in the HTML details
        assert [entry["identifier"] for entry in data] == [
            "subclusters-svg-r2c3-sc1",
            "subclusters-svg-r2c3-sc2",
        ]

    def test_hits_by_cds(self):
        hits = [build_hit("cds1", "A"), build_hit("cds2", "A", description=None),
                build_hit("cds1", "B")]
        cds_results = self.generate([build_prediction(FeatureLocation(150, 250, 1), hits)]
                                    )[0]["cds_results"]

        assert [result["cds"] for result in cds_results] == ["cds1", "cds2"]
        assert cds_results[0]["domains"] == [
            {"name": "A", "description": "desc", "accession": "A.1", "evalue": 1e-10,
             "bitscore": 50.},
            {"name": "B", "description": "desc", "accession": "B.1", "evalue": 1e-10,
             "bitscore": 50.},
        ]
        assert [domain["name"] for domain in cds_results[1]["domains"]] == ["A"]

        # each tooltip names the CDS and its domains
        def span(class_name, text):
            return f'<span class="subclusters-tooltip-{class_name}">{text}</span>'

        first, second = (result["tooltip"] for result in cds_results)
        assert span("cds-name", "cds1") in first
        assert span("domain-name", "A") in first
        assert span("domain-name", "B") in first
        assert span("cds-name", "cds2") in second
        assert span("domain-name", "A") in second
        assert span("domain-name", "B") not in second
        assert span("domain-description", "no description available") in second

    def extends(self, location):
        data = self.generate([build_prediction(location)])
        return data[0]["extends_left"], data[0]["extends_right"]

    def test_contained(self):
        assert self.extends(FeatureLocation(150, 250, 1)) == (False, False)

    def test_same_as_region(self):
        # the end of a location is exclusive, so ending with the region doesn't extend it
        assert self.extends(FeatureLocation(100, 500, 1)) == (False, False)

    def test_extends_left(self):
        assert self.extends(FeatureLocation(50, 250, 1)) == (True, False)
        assert self.extends(FeatureLocation(99, 250, 1)) == (True, False)

    def test_extends_right(self):
        assert self.extends(FeatureLocation(150, 550, 1)) == (False, True)
        assert self.extends(FeatureLocation(150, 501, 1)) == (False, True)

    def test_extends_both(self):
        assert self.extends(FeatureLocation(50, 550, 1)) == (True, True)

    def test_cross_origin(self):
        # a region crossing the origin of a circular record of length 1000
        self.region.location = CompoundLocation([FeatureLocation(900, 1000, 1),
                                                 FeatureLocation(0, 100, 1)])
        inside = CompoundLocation([FeatureLocation(950, 1000, 1), FeatureLocation(0, 50, 1)])
        assert self.extends(inside) == (False, False)
        left = CompoundLocation([FeatureLocation(850, 1000, 1), FeatureLocation(0, 50, 1)])
        assert self.extends(left) == (True, False)
        right = CompoundLocation([FeatureLocation(950, 1000, 1), FeatureLocation(0, 150, 1)])
        assert self.extends(right) == (False, True)


class TestAnchorsMatch(unittest.TestCase):
    def test_javascript_identifiers_in_html(self):
        # the javascript draws each subcluster into the HTML element with its
        # identifier, so every identifier must exist in the HTML details
        record = SimpleNamespace(record_index=2)
        region = SimpleNamespace(parent_record=record, get_region_number=lambda: 3,
                                 location=FeatureLocation(100, 500, 1))
        region_layer = SimpleNamespace(region_feature=region,
                                       anchor_id=RegionLayer.build_anchor_id(region))
        results = build_results([build_prediction(FeatureLocation(150, 250, 1)),
                                 build_prediction(FeatureLocation(300, 400, 1))])

        data = html_output.generate_javascript_data(record, region, results)
        html = html_output.generate_html(region_layer, results, None, None)
        html_ids = set(re.findall(r'id="([^"]+)"', str(html.detail_sections[0].content)))

        identifiers = [entry["identifier"] for entry in data]
        assert len(identifiers) == 2
        for identifier in identifiers:
            assert identifier in html_ids, identifier
