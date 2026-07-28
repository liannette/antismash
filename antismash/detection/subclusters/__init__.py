# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Detection of subclusters
"""
import logging
import os
from typing import Any, Optional

from antismash.common import hmmer, path
from antismash.common.hmm_rule_parser.cluster_prediction import detect_protoclusters_and_signatures
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from .results import SubclusterDetectionResults
from .ruleset import get_ruleset, _STRICTNESS_LEVELS
from .signatures import get_signature_profiles
from .html_output import generate_html, will_handle, generate_javascript_data

NAME = "subclusters"
SHORT_DESCRIPTION = "Detects subclusters for specific chemical moieties, e.g. precursor molecules."
DETECTION_STAGE = DetectionStage.AREA_FORMATION


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module

        Returns:
            an empty or populated set of module arguments
    """
    args = ModuleArgs("Subcluster detection options", "subclusters")

    args.add_analysis_toggle('--subclusters',
                             dest='subclusters',
                             default=False,
                             action='store_true',
                             help="Run Subcluster detection.")
    args.add_option('strictness',
                    dest='strictness',
                    type=str,
                    choices=list(_STRICTNESS_LEVELS),
                    default="relaxed",
                    help=("Defines which level of strictness to use for "
                          "subcluster detection. Levels are cumulative, so "
                          "looser levels also include all stricter rules "
                          "(default: %(default)s)."))
    args.add_option('as-subregions',
                    dest='as_subregions',
                    action='store_true',
                    default=False,
                    help=("Treat detected subclusters as subregions so they "
                          "seed and extend regions during region formation, "
                          "(default: %(default)s)."))
    args.add_option('require-overlap',
                    dest='require_overlap',
                    action='store_true',
                    default=False,
                    help=("Only relevant with --subclusters-as-subregions: keep "
                          "a subcluster subregion only if it overlaps a cluster "
                          "found by another detection module, so subclusters can "
                          "extend existing regions but never create new regions "
                          "on their own or merge with each other "
                          "(default: %(default)s)."))
    return args


def check_options(options: ConfigType) -> list[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    if options.subclusters_strictness not in _STRICTNESS_LEVELS:
        return [f"Unknown subcluster strictness level: {options.subclusters_strictness}"]

    if options.subclusters_require_overlap and not options.subclusters_as_subregions:
        logging.warning("--subclusters-require-overlap has no effect without "
                        "--subclusters-as-subregions; subclusters will not be "
                        "treated as subregions")

    return []


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return options.subclusters


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared.

        Aggregates the individual subcluster HMM signatures into a single
        combined profile database and presses it with hmmpress, regenerating
        whenever it is missing or out of date. 

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failure_messages: list[str] = []

    # Check that hmmdetails.txt is readable and well-formatted
    try:
        signatures = get_signature_profiles()
    except ValueError as err:
        if not logging_only:
            raise
        return [str(err)]

    aggregate_hmm = path.get_full_path(__file__, "data", "subclusters.hmm")
    hmm_files = [signature.hmm_file for signature in signatures]

    description_file = path.get_full_path(__file__, "data", "hmmdetails.txt")
    force_replace = not (path.locate_file(aggregate_hmm)
                         and os.path.getmtime(description_file) < os.path.getmtime(aggregate_hmm))

    failure_messages.extend(hmmer.aggregate_profiles(aggregate_hmm, hmm_files, force_replace=force_replace,
                                                     return_not_raise=logging_only))

    return failure_messages


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check that prereqs are satisfied. hmmpress is only required if the
        databases have not yet been generated.
    """
    failure_messages = []
    for binary_name in ["hmmsearch", "hmmpress"]:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable for {binary_name!r}")

    # no point checking the data if we can't use it
    if failure_messages:
        return failure_messages

    failure_messages.extend(prepare_data(logging_only=True))

    return failure_messages


def _get_strictness(options: ConfigType) -> str:
    """ Returns the subcluster detection strictness to use for the given options. """
    return options.subclusters_strictness


def regenerate_previous_results(results: dict[str, Any], record: Record,
                                options: ConfigType) -> Optional[SubclusterDetectionResults]:
    """Regenerate previous results."""
    return None
    if not results:
        return None
    previous = SubclusterDetectionResults.from_json(results, record)
    if previous is None:
        return None

    current_strictness = _get_strictness(options)
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

    current_strictness = _get_strictness(options)
    ruleset = get_ruleset(current_strictness)
    rule_results = detect_protoclusters_and_signatures(record, ruleset)

    # The shared rule-based detection pipeline hardcodes the protocluster
    # "aStool" qualifier to "rule-based-clusters". These protoclusters were
    # produced by this module, so relabel them with this module's tool name
    # to avoid the misleading qualifier in the output.
    for protocluster in rule_results.protoclusters:
        protocluster.tool = rule_results.tool

    return SubclusterDetectionResults(
        record_id=record.id,
        rule_results=rule_results,
        rule_names=ruleset.get_rule_names(),
        strictness=current_strictness,
        as_subregions=options.subclusters_as_subregions,
        require_overlap=options.subclusters_require_overlap,
        record=record,
    )

