import logging
from typing import Any, Optional, Self

from antismash.common.hmm_rule_parser.cluster_prediction import RuleDetectionResults
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Record, Region, SubRegion

from .predictions import SubclusterPrediction, build_predictions
from .subregions import SubRegionMode, build_subregions


class SubclusterDetectionResults(DetectionResults):
    """Results class for the Subcluster detection module """

    schema_version = 1  # increment when the JSON format changes

    def __init__(
            self,
            record_id: str,
            rule_results: RuleDetectionResults,
            rule_names: set[str],
            strictness: str,
            predictions: list[SubclusterPrediction],
            subregions: list[SubRegion],
            subregion_mode: SubRegionMode,
    ) -> None:

        if subregion_mode not in SubRegionMode:
            raise ValueError(f"Unknown subcluster subregion mode: {subregion_mode!r}")

        super().__init__(record_id)
        self.rule_results = rule_results
        self.rule_names = rule_names
        self.strictness = strictness
        self.subregion_mode = subregion_mode
        self.predictions = predictions
        self.subregions = subregions

    def get_predictions_for_region(self, region: Region) -> list[SubclusterPrediction]:
        """Return all predictions overlapping the given region."""
        return [
            prediction for prediction in self.predictions
            if region.overlaps_with(prediction.location)
        ]

    def get_predicted_subregions(self) -> list[SubRegion]:
        """Return sub-region features for the detected subclusters.

        These are built when detection runs, since they depend on the areas
        found by the other detection modules, and are added to the record during
        the region-formation step of the main pipeline.
        """
        return self.subregions

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "strictness": self.strictness,
            "subregion_mode": self.subregion_mode,
            "rule_names": sorted(self.rule_names),
            "rule_results": self.rule_results.to_json(),
        }

    @classmethod
    def from_json(cls, data: dict[str, Any], record: Record) -> Optional[Self]:
        if data.get("schema_version") != cls.schema_version:
            logging.debug(
                "Discarding subcluster results: schema version %s != %s",
                data.get("schema_version"), cls.schema_version,
            )
            return None
        
        rule_results = RuleDetectionResults.from_json(data["rule_results"], record)
        if rule_results is None:
            logging.debug("Discarding subcluster results: rule detection result schema changed")
            return None
        
        # neither the predictions nor the subregions are stored, as they're
        # derived from the rule results and, for the subregions, the areas found
        # by the other detection modules, so both are rebuilt here
        subregion_mode = SubRegionMode(data["subregion_mode"])
        return cls(
            record_id=data["record_id"],
            rule_results=rule_results,
            rule_names=set(data["rule_names"]),
            strictness=data["strictness"],
            subregion_mode=subregion_mode,
            predictions=build_predictions(rule_results, data["strictness"]),
            subregions=build_subregions(record, rule_results.protoclusters,
                                        tool=rule_results.tool, mode=subregion_mode),
        )

    def add_to_record(self, record: Record) -> None:
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        self.rule_results.annotate_cds_features()
