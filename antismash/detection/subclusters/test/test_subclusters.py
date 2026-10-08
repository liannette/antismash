# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from types import SimpleNamespace
import unittest
from unittest.mock import patch, sentinel

from antismash.detection import subclusters
from antismash.detection.subclusters.results import SubclusterDetectionResults
from antismash.detection.subclusters.ruleset import get_ruleset


def build_options(strictness="relaxed", mode="clip"):
    return SimpleNamespace(subclusters_strictness=strictness,
                           subclusters_subregion_mode=mode)


def build_results(strictness="relaxed", mode="clip", rule_names=None):
    if rule_names is None:
        rule_names = sorted(get_ruleset(strictness).get_rule_names())
    return {
        "schema_version": SubclusterDetectionResults.schema_version,
        "strictness": strictness,
        "subregion_mode": mode,
        "rule_names": rule_names,
    }


class TestRegenerate(unittest.TestCase):
    def regenerate(self, results, options):
        return subclusters.regenerate_previous_results(results, sentinel.record, options)

    def test_empty(self):
        assert self.regenerate({}, build_options()) is None

    def test_unchanged(self):
        with patch.object(SubclusterDetectionResults, "from_json",
                          return_value=sentinel.regenerated) as patched:
            for mode in ["clip", "extend", "create"]:
                results = build_results(mode=mode)
                assert self.regenerate(results, build_options(mode=mode)) is sentinel.regenerated
                patched.assert_called_with(results, sentinel.record)

    def test_mode_changed(self):
        for previous, current in [("clip", "extend"), ("extend", "create"), ("create", "clip")]:
            with self.assertRaisesRegex(RuntimeError, "region mode changed"):
                self.regenerate(build_results(mode=previous), build_options(mode=current))

    def test_mode_missing(self):
        results = build_results()
        del results["subregion_mode"]
        with self.assertRaisesRegex(RuntimeError, "region mode changed"):
            self.regenerate(results, build_options())

    def test_schema_changed_clip(self):
        results = build_results()
        results["schema_version"] = -1
        assert self.regenerate(results, build_options()) is None

    def test_schema_changed_region_altering(self):
        for mode in ["extend", "create"]:
            results = build_results(mode=mode)
            results["schema_version"] = -1
            with self.assertRaisesRegex(RuntimeError, "results format changed"):
                self.regenerate(results, build_options(mode=mode))

    def test_rule_results_schema_changed_clip(self):
        with patch.object(SubclusterDetectionResults, "from_json", return_value=None):
            assert self.regenerate(build_results(), build_options()) is None

    def test_rule_results_schema_changed_region_altering(self):
        with patch.object(SubclusterDetectionResults, "from_json", return_value=None):
            for mode in ["extend", "create"]:
                with self.assertRaisesRegex(RuntimeError, "rule detection results format changed"):
                    self.regenerate(build_results(mode=mode), build_options(mode=mode))

    def test_strictness_changed_clip(self):
        results = build_results(strictness="strict")
        assert self.regenerate(results, build_options(strictness="relaxed")) is None

    def test_strictness_changed_region_altering(self):
        for mode in ["extend", "create"]:
            results = build_results(strictness="strict", mode=mode)
            with self.assertRaisesRegex(RuntimeError, "strictness or rules changed"):
                self.regenerate(results, build_options(strictness="relaxed", mode=mode))

    def test_rules_changed_clip(self):
        results = build_results(rule_names=["not-a-rule"])
        assert self.regenerate(results, build_options()) is None

    def test_rules_changed_region_altering(self):
        for mode in ["extend", "create"]:
            results = build_results(mode=mode, rule_names=["not-a-rule"])
            with self.assertRaisesRegex(RuntimeError, "strictness or rules changed"):
                self.regenerate(results, build_options(mode=mode))
