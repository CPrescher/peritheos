"""Private audit manifest for the pre-0.7 convenience-record catalog."""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, Mapping

LegacyRecordRelationship = Literal["equivalent", "different", "absent"]


@dataclass(frozen=True)
class LegacyRecordAudit:
    """Relationship between one historical record and the canonical catalog."""

    relationship: LegacyRecordRelationship
    canonical_identifier: str | None


def _entry(
    relationship: LegacyRecordRelationship, canonical_identifier: str | None = None
) -> LegacyRecordAudit:
    return LegacyRecordAudit(relationship, canonical_identifier)


# This explicit manifest is intentionally exhaustive. ``test_catalog`` checks
# it against the historical compatibility mapping and canonical identifiers so
# adding, deleting, or redirecting a convenience record cannot silently bypass
# the scientific migration audit.
LEGACY_RECORD_AUDIT: Mapping[str, LegacyRecordAudit] = MappingProxyType(
    {
        "mgo_b1_tange_2009_vinet": _entry("absent"),
        "mgo_b1_sokolova_2013": _entry("different", "mgo_sokolova_2013_holzapfel_4"),
        "diamond_sokolova_2013": _entry(
            "different", "diamond_sokolova_2013_holzapfel_3"
        ),
        "al_fcc_sokolova_2013": _entry(
            "different", "aluminum_sokolova_2013_holzapfel_2"
        ),
        "cu_fcc_sokolova_2013": _entry("different", "copper_sokolova_2013_holzapfel_2"),
        "ag_fcc_sokolova_2013": _entry("different", "silver_sokolova_2013_holzapfel_2"),
        "au_fcc_sokolova_2013": _entry("different", "gold_sokolova_2013_holzapfel_4"),
        "pt_fcc_sokolova_2013": _entry(
            "different", "platinum_sokolova_2013_holzapfel_3"
        ),
        "nb_bcc_sokolova_2013": _entry(
            "different", "niobium_sokolova_2013_holzapfel_2"
        ),
        "ta_bcc_sokolova_2013": _entry(
            "different", "tantalum_sokolova_2013_holzapfel_3"
        ),
        "mo_bcc_sokolova_2013": _entry(
            "different", "molybdenum_sokolova_2013_holzapfel_2"
        ),
        "w_bcc_sokolova_2013": _entry(
            "different", "tungsten_sokolova_2013_holzapfel_4"
        ),
        "au_fcc_fei_2007": _entry("equivalent", "gold_fei_2007_vinet_2"),
        "pt_fcc_fei_2007": _entry("absent"),
        "nacl_b2_fei_2007": _entry("absent"),
        "ne_fcc_fei_2007": _entry("equivalent", "neon_fcc_fei_2007_vinet_2"),
        "au_fcc_dorfman_2012": _entry("absent"),
        "pt_fcc_dorfman_2012": _entry("absent"),
        "mo_bcc_dorfman_2012": _entry("absent"),
        "nacl_b2_dorfman_2012": _entry("absent"),
        "ne_fcc_dorfman_2012": _entry("absent"),
        "lif_b1_dewaele_2019": _entry("different", "lif_b1_dewaele_2019_vinet_1"),
        "nacl_b1_dewaele_2019": _entry("absent"),
        "nacl_b2_dewaele_2019": _entry("absent"),
        "kcl_b1_dewaele_2012": _entry("equivalent", "kcl_b1_dewaele_2012_vinet_1"),
        "kcl_b2_dewaele_2012": _entry("equivalent", "kcl_b2_dewaele_2012_vinet_3"),
        "kbr_b1_dewaele_2012": _entry("equivalent", "kbr_b1_dewaele_2012_vinet_1"),
        "kbr_b2_dewaele_2012": _entry("different", "kbr_b2_dewaele_2012_vinet_1"),
        "cbn_zincblende_datchi_2007": _entry(
            "different", "boron_nitride_datchi_2007_vinet_1"
        ),
        "diamond_correa_2008": _entry(
            "equivalent", "diamond_correa_2008_double_debye_log_moment_5"
        ),
        "diamond_benedict_2014": _entry(
            "equivalent", "diamond_benedict_2014_double_debye_4"
        ),
        "diamond_datchi_dewaele_2008": _entry("absent"),
        "diamond_correa_2008_dewaele_anchored": _entry(
            "equivalent", "diamond_correa_2008_dewaele_anchored"
        ),
        "diamond_benedict_2014_dewaele_anchored": _entry(
            "equivalent", "diamond_benedict_2014_dewaele_anchored"
        ),
        "ni_fcc_dewaele_2008": _entry("different", "nickel_dewaele_2008_vinet_1"),
        "ag_fcc_dewaele_2008": _entry("different", "silver_dewaele_2008_vinet_1"),
        "re_hcp_anzellini_2014": _entry("different", "rhenium_anzellini_2014_vinet_1"),
    }
)


__all__ = ["LEGACY_RECORD_AUDIT", "LegacyRecordAudit", "LegacyRecordRelationship"]
