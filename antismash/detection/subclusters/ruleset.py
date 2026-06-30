"""Ruleset construction for the subcluster detection module."""

from antismash.common import path
from antismash.common.hmm_rule_parser import rule_parser
from antismash.common.hmm_rule_parser.cluster_prediction import Ruleset

from .signatures import get_subcluster_profiles

# Strictness levels are cumulative: each level includes all levels above it.
# Add new levels here and provide a matching rule file in subcluster_rules/.
_STRICTNESS_LEVELS = ("strict",)

# Categories referenced across all rule files, should always be "Subcluster"
_CATEGORIES = {"Subcluster"}

_HMM_FILE = path.get_full_path(__file__, "data", "subclusters.hmm")

_CACHE: dict[str, Ruleset] = {}


def get_ruleset(strictness: str = "strict") -> Ruleset:
    """Return a Ruleset for the given strictness level.

    Results are cached per strictness level.

    Arguments:
        strictness: one of the levels defined in ``_STRICTNESS_LEVELS``

    Returns:
        a configured ``Ruleset`` instance
    """
    if strictness not in _STRICTNESS_LEVELS:
        raise ValueError(f"Unknown strictness level {strictness!r}. "
                         f"Valid levels: {', '.join(_STRICTNESS_LEVELS)}")

    if strictness in _CACHE:
        return _CACHE[strictness]

    ruleset = _build_ruleset(strictness)
    _CACHE[strictness] = ruleset
    return ruleset


def _rule_files_for_strictness(strictness: str) -> list[str]:
    """Return the ordered list of rule files to load for a given strictness level.

    Files are loaded in order from most strict to the requested level, so that
    rules from stricter levels are always included.
    """
    levels = _STRICTNESS_LEVELS[: _STRICTNESS_LEVELS.index(strictness) + 1]
    return [path.get_full_path(__file__, "subcluster_rules", f"{level}.txt") for level in levels]


def _build_ruleset(strictness: str) -> Ruleset:
    profiles = get_subcluster_profiles()

    rules: list[rule_parser.DetectionRule] = []
    aliases: dict[str, list[rule_parser.Token]] = {}
    for rule_file in _rule_files_for_strictness(strictness):
        with open(rule_file, encoding="utf-8") as f:
            rules = rule_parser.Parser(
                f.read(), set(profiles), _CATEGORIES, rules, aliases,
            ).rules

    return Ruleset(
        tuple(rules),
        dict(profiles),
        _HMM_FILE,
        _CATEGORIES,
        "subclusters",
        equivalence_groups=[],
    )