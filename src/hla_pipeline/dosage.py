"""Readers for per-sample HLA allele dosages.

Two layouts are supported:

* PLINK ``--recode A`` (``.raw``): one row per sample, whitespace-delimited,
  header ``FID IID PAT MAT SEX PHENOTYPE <marker>_<counted allele> ...``.
  PLINK counts the A1 allele, which is not always the one you want.
* SNP2HLA ``.dosage``: one row per marker, ``marker allele1 allele2`` then one
  value per sample in ``.fam`` order, equal to ``2*P(11) + P(12)``, the
  expected count of allele1.

SNP2HLA codes HLA allele markers (``HLA_DRB1_0401`` and so on) as binary
presence (``P``) / absence (``A``). Both readers return the expected count of
the *present* allele, flipping ``x -> 2 - x`` when the counted allele is ``A``.
That orientation is only applied to ``HLA_`` markers, because ``A`` means
adenine for ordinary SNPs.
"""

from __future__ import annotations

import math
from pathlib import Path

PRESENT = "P"
ABSENT = "A"


def _orient(marker: str, counted_allele: str, value: float) -> float:
    if marker.startswith("HLA_") and counted_allele == ABSENT:
        return 2.0 - value
    return value


def _parse_dosage_value(raw: str, marker: str, sample_id: str) -> float:
    val = float(raw)
    if not math.isfinite(val):
        raise ValueError(
            f"{marker}: non-finite dosage value {raw!r} for sample {sample_id!r}"
        )
    return val


def find_raw_column(header: list[str], marker: str) -> tuple[int, str]:
    """Locate ``marker`` in a .raw header.

    Accepts the bare marker name or PLINK's ``<marker>_<allele>`` form and
    returns ``(column index, counted allele)``; the allele is "" for a bare
    name.
    """
    for i, col in enumerate(header):
        if col == marker:
            return i, ""
        if col.startswith(marker + "_"):
            allele = col[len(marker) + 1 :]
            if allele.isalpha():
                return i, allele
    raise KeyError(f"marker {marker!r} not found in header")


def read_plink_raw(path: str | Path, marker: str) -> dict[str, float]:
    """Return ``{IID: dosage of the present allele}`` from a PLINK .raw file.

    Samples with a missing value (``NA``) are omitted.
    """
    with open(path) as fh:
        header = fh.readline().split()
        col, allele = find_raw_column(header, marker)
        iid_col = header.index("IID") if "IID" in header else 1
        out: dict[str, float] = {}
        for line in fh:
            parts = line.split()
            if len(parts) <= max(col, iid_col) or parts[col].upper() == "NA":
                continue
            val = _parse_dosage_value(parts[col], marker, parts[iid_col])
            out[parts[iid_col]] = _orient(marker, allele, val)
    return out


def read_fam_iids(fam_path: str | Path) -> list[str]:
    """Return the IIDs (column 2) of a PLINK .fam file, in order."""
    iids: list[str] = []
    with open(fam_path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(
                    f"{fam_path}:{line_no}: expected at least 2 columns in .fam"
                )
            iids.append(parts[1])
    return iids


def read_snp2hla_dosage(
    dosage_path: str | Path, fam_path: str | Path, marker: str
) -> dict[str, float]:
    """Return ``{IID: dosage of the present allele}`` from a SNP2HLA .dosage.

    ``fam_path`` must be the .fam written by the same SNP2HLA run, because the
    dosage columns carry no sample IDs. Samples with a missing value (``NA``)
    are omitted.
    """
    iids = read_fam_iids(fam_path)
    with open(dosage_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3 or parts[0] != marker:
                continue
            allele1 = parts[1]
            values = parts[3:]
            if len(values) != len(iids):
                raise ValueError(
                    f"{marker}: {len(values)} dosage values but {len(iids)} "
                    "samples in the .fam file"
                )
            out: dict[str, float] = {}
            for iid, raw in zip(iids, values, strict=True):
                if raw.upper() == "NA":
                    continue
                val = _parse_dosage_value(raw, marker, iid)
                out[iid] = _orient(marker, allele1, val)
            return out
    raise KeyError(f"marker {marker!r} not found in {dosage_path}")


def read_r2(path: str | Path) -> dict[str, float]:
    """Parse a Beagle ``.bgl.r2`` file into ``{marker: r2}``.

    Non-numeric or non-finite values (for example ``NaN`` for monomorphic
    markers in older Beagle versions) are skipped.
    """
    scores: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                value = float(parts[1])
            except ValueError:
                continue
            if math.isfinite(value):
                scores[parts[0]] = value
    return scores
