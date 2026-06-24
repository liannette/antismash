# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Detection of subclusters
"""

from typing import Any, Optional

from antismash.common.hmm_rule_parser.cluster_prediction import Ruleset
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from .results import SubclusterDetectionResults, SubclusterPrediction
from .signatures import get_subcluster_profiles
from .html_output import generate_html, will_handle, generate_javascript_data

NAME = "subclusters"
SHORT_DESCRIPTION = "Detects subclusters for specific chemical moieties, e.g. precursor molecules."
DETECTION_STAGE = DetectionStage.AREA_FORMATION


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module

        Returns:
            an empty or populated ModuleArgs instance
    """
    args = ModuleArgs("Subcluster detection options", "subclusters")

    args.add_analysis_toggle('--subclusters',
                             dest='subclusters',
                             default=False,
                             action='store_true',
                             help="Run Subcluster detection.")
    return args


def check_options(options: ConfigType) -> list[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    return []


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return options.subclusters


def check_prereqs(_options: ConfigType) -> list[str]:
    """ Checks that prerequisites are satisfied.
    """
    return []


def regenerate_previous_results(results: dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[SubclusterDetectionResults]:
    """ Regenerate previous results. 
    """
    if not results:
        return None

    return SubclusterDetectionResults.from_json(results)


def run_on_record(record: Record, previous_results: Optional[SubclusterDetectionResults],
                  options: ConfigType) -> Optional[SubclusterDetectionResults]:
    """ Finds the external annoations for the given record """
    if previous_results:
        return previous_results
    
    # # TODO: build ruleset from rule/HMM files
    # ruleset: Ruleset = _build_ruleset(options)

    # hits: list[SubclusterPrediction] = [
    #     hit.enrich(ruleset) for hit in _detect_hits(record, ruleset)
    # ]

    # return SubclusterDetectionResults(
    #     record_id=record.id,
    #     rules_version=ruleset.tool,
    #     hits=hits,
    # )

    return SubclusterDetectionResults(record.id, "ruleset_version_placeholder", [])  # placeholder until real detection is implemented



def _build_ruleset(options: ConfigType) -> Ruleset:
    """Construct the Ruleset from rule/HMM files.  Placeholder."""
    raise NotImplementedError
 
 
def _detect_hits(record: Record, ruleset: Ruleset) -> list[SubclusterPrediction]:
    """Run rule-based detection for one region.  Placeholder."""
    raise NotImplementedError

