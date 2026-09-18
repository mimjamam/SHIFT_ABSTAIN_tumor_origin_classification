"""
COSMIC-style SBS-96 trinucleotide substitution categories.

Convention: every substitution is expressed with a pyrimidine (C or T)
reference base. If the reported reference base is a purine (A or G), the
mutation and its flanking bases are reverse-complemented before counting.
This is the standard convention used by COSMIC / SigProfiler mutational
signatures.
"""

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

BASES = ["A", "C", "G", "T"]
PYR_SUBS = {
    "C": ["A", "G", "T"],
    "T": ["A", "C", "G"],
}

# Canonical ordering: substitution type (C>A, C>G, C>T, T>A, T>C, T>G),
# then 5' base, then 3' base, each ascending over A,C,G,T.
CATEGORIES = []
for ref in ["C", "T"]:
    for alt in PYR_SUBS[ref]:
        for five in BASES:
            for three in BASES:
                CATEGORIES.append(f"{five}[{ref}>{alt}]{three}")

assert len(CATEGORIES) == 96
CATEGORY_INDEX = {cat: i for i, cat in enumerate(CATEGORIES)}


def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def sbs96_category(context3: str, ref: str, alt: str):
    """context3: the 3 bases (5', ref, 3') around the mutated site, ref/alt:
    the reported single-base substitution. Returns the SBS-96 category
    string, or None if inputs aren't a clean single-base substitution with
    valid ACGT bases."""
    if len(context3) != 3 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        return None
    if any(b not in "ACGT" for b in context3):
        return None
    if context3[1] != ref:
        return None  # reference mismatch between MAF and genome lookup
    if ref in ("A", "G"):
        context3 = revcomp(context3)
        ref = COMPLEMENT[ref]
        alt = COMPLEMENT[alt]
    five, _, three = context3
    return f"{five}[{ref}>{alt}]{three}"
