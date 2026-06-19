# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of subclustersers
"""

import argparse
import os
from typing import Any, AnyStr, Dict, List, Optional

from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs, MultipleFullPathAction, SplitCommaAction
from antismash.detection import DetectionStage
from .results import SubclusterDetectionResults
from .html_output import generate_html, will_handle, generate_javascript_data

NAME = "subclusters"
SHORT_DESCRIPTION = ""
DETECTION_STAGE = DetectionStage.AREA_FORMATION


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module
    """
    args = ModuleArgs("Subcluster options", "subclusters")
    args.add_option("subclusters",
                    dest="subclusters",
                    type=bool,
                    default=False,
                    help="help me")
    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    return []


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return options.subclusters


def regenerate_previous_results(results: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[SubclusterDetectionResults]:
    """ Regenerate previous results. """
    if not results:
        return None

    return SubclusterDetectionResults.from_json(results, record)


def run_on_record(record: Record, previous_results: Optional[SubclusterDetectionResults],
                  options: ConfigType) -> Optional[DetectionResults]:
    """ Finds the external annoations for the given record """
    if previous_results:
        return previous_results
    import logging; logging.critical("%s", SubclusterDetectionResults(record.id).to_json())
    return SubclusterDetectionResults(record.id)

def check_prereqs(_options: ConfigType) -> List[str]:
    """ Checks that prerequisites are satisfied.
    """
    return []
