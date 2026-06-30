# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Detection of subclusters
"""
import logging
from typing import Any, Optional

from antismash.common.hmm_rule_parser.cluster_prediction import detect_protoclusters_and_signatures
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from .results import SubclusterDetectionResults
from .ruleset import get_ruleset
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
                                options: ConfigType) -> Optional[SubclusterDetectionResults]:
    """Regenerate previous results."""
    return None
    if not results:
        return None
    previous = SubclusterDetectionResults.from_json(results, record)
    if previous is None:
        return None

    # TODO: replace "strict" with options.subclusters_strictness once that option exists
    current_strictness = "strict"
    if previous.strictness != current_strictness:
        logging.debug("Subcluster strictness changed from %r to %r; forcing re-detection.",
                      previous.strictness, current_strictness)
        return None

    current_rule_names = get_ruleset(current_strictness).get_rule_names()
    if previous.rule_names != current_rule_names:
        logging.debug("Subcluster rules changed; forcing re-detection.")
        return None

    return previous


def run_on_record(record: Record, previous_results: Optional[SubclusterDetectionResults],
                  options: ConfigType) -> SubclusterDetectionResults:
    """Run subcluster detection on a single record."""
    if previous_results:
        return previous_results

    # TODO: replace "strict" with options.subclusters_strictness once that option exists
    current_strictness = "strict"
    ruleset = get_ruleset(current_strictness)
    rule_results = detect_protoclusters_and_signatures(record, ruleset)

    return SubclusterDetectionResults(
        record_id=record.id,
        rule_results=rule_results,
        rule_names=ruleset.get_rule_names(),
        strictness=current_strictness,
    )

