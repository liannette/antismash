# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Construction of the view objects describing each detected subcluster.

Predictions are a plain derivation of the rule detection results and the ruleset
detection ran with, so they aren't serialised with the results and are instead
rebuilt whenever those results are reconstructed. Every compound and signature
lookup happens here, so the results themselves never have to reach back into the
packaged details files.
"""
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

from antismash.common.hmm_rule_parser.rule_parser import DetectionRule
from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults, RuleDetectionResults
from antismash.common.secmet.locations import FeatureLocation

from .compounds import CompoundInfo, get_compound
from .ruleset import get_ruleset
from .signatures import get_signatures


@dataclass(frozen=True)
class CDSDomainHit:
    """A single signature match found in a specific CDS."""
    domain_name: str
    domain_description: Optional[str]
    domain_accession: Optional[str]
    cds_name: str
    evalue: float
    bitscore: float


class SubclusterPrediction:
    """A single predicted subcluster, used as a view object for HTML rendering.

    Attributes:
        rule: The detection rule whose conditions were met
        location: location of the protocluster that produced this prediction,
            i.e. the matching core along with any neighbourhood the rule defines
        cds_results: the per-CDS detection results for every CDS that contributed
            to this prediction, as returned by the detection pipeline
        compound: the compound associated with the detection rule
        domain_hits: a flat list of every domain hit, each paired with its CDS name
    """

    def __init__(
            self,
            *,
            rule: DetectionRule,
            location: FeatureLocation,
            cds_results: list[CDSResults],
            compound: CompoundInfo,
            domain_hits: list[CDSDomainHit],
    ) -> None:
        self.rule = rule
        self.location = location
        self.cds_results = cds_results
        self.compound = compound
        self.domain_hits = domain_hits

    @property
    def rule_name(self) -> str:
        """Name of the matching detection rule."""
        return self.rule.name

    @property
    def conditions_str(self) -> str:
        text = str(self.rule.conditions)
        if text.startswith("(") and text.endswith(")"):
            return text[1:-1]
        return text

    @property
    def domain_hits_by_cds(self) -> dict[str, list[CDSDomainHit]]:
        hits_by_cds: dict[str, list[CDSDomainHit]] = defaultdict(list)
        for hit in self.domain_hits:
            hits_by_cds[hit.cds_name].append(hit)
        return dict(hits_by_cds)

    @property
    def domain_hits_by_domain(self) -> dict[str, list[CDSDomainHit]]:
        hits_by_domain: dict[str, list[CDSDomainHit]] = defaultdict(list)
        for hit in self.domain_hits:
            hits_by_domain[hit.domain_name].append(hit)
        return dict(hits_by_domain)

    def __repr__(self) -> str:
        return (
            f"SubclusterPrediction(rule_name={self.rule_name!r}, "
            f"location={self.location.start}-{self.location.end}, "
            f"cds_count={len(self.cds_results)})"
        )


def _build_domain_hits(rule_name: str, cds_results: list[CDSResults]) -> list[CDSDomainHit]:
    """The domain hits that fired the given rule, flattened across the given CDS results.

    Where a domain matched a CDS more than once, only the strongest hit is kept.
    """
    hits: list[CDSDomainHit] = []
    for cds_result in cds_results:
        cds_name = cds_result.cds.get_name()
        fired = cds_result.definition_domains.get(rule_name, set())
        for domain_name in sorted(fired):
            matching = [domain for domain in cds_result.domains if domain.name == domain_name]
            if not matching:
                continue
            best = max(matching, key=lambda domain: domain.bitscore)
            signature = get_signatures()[domain_name]
            hits.append(CDSDomainHit(
                domain_name=signature.name,
                domain_description=signature.description,
                domain_accession=signature.accession,
                cds_name=cds_name,
                evalue=best.evalue,
                bitscore=best.bitscore,
            ))
    return hits


def build_predictions(rule_results: RuleDetectionResults, strictness: str,
                      ) -> list[SubclusterPrediction]:
    """Build the view objects for each detected subcluster.

    Arguments:
        rule_results: the detection results to build the predictions from
        strictness: the strictness level detection was run at, which determines
            the ruleset the matching rules are taken from

    Returns:
        a list of predictions, one per detected subcluster
    """
    ruleset = get_ruleset(strictness)
    predictions: list[SubclusterPrediction] = []
    for protocluster, cds_results in rule_results.cds_by_cluster.items():
        rule = ruleset.get_rule_by_name(protocluster.product)
        # only the CDS features actually contributing to this rule are relevant,
        # since a CDS within the cutoff may only define a neighbouring subcluster
        contributing = [result for result in cds_results
                        if result.definition_domains.get(protocluster.product)]
        predictions.append(SubclusterPrediction(
            rule=rule,
            location=protocluster.location,
            cds_results=contributing,
            # every rule must have a matching entry in the compound details file
            compound=get_compound(rule.name),
            domain_hits=_build_domain_hits(rule.name, contributing),
        ))
    return predictions
