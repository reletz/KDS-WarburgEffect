"""M2 — Parameter provenance for three papers × three organisms.

Every numeric constant is tagged with its source (paper, table/figure, organism).
Run as: python -m src.params --table
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import pandas as pd

Organism = Literal["ecoli", "yeast", "mammalian"]
Paper = Literal["shen2024", "kukurugya2024", "wang2025"]


@dataclass(frozen=True)
class ParameterEntry:
    """A single numeric parameter with full provenance."""

    value: float
    unit: str
    source_paper: Paper
    source_location: str  # e.g. "Table 1", "SI Table S2", "Fig 1F"
    organism: Organism
    note: str = ""


@dataclass(frozen=True)
class OrganismParams:
    """Complete parameter set for one organism, aligned to shared interface.

    Parameters (Kukurugya notation):
        gamma_resp: ATP yield per glucose via respiration (mol ATP / mol glucose)
        gamma_glyc: ATP yield per glucose via glycolysis/fermentation (mol ATP / mol glucose)
        V_resp: specific activity of respiratory enzymes (µmol substrate·min⁻¹·mg protein⁻¹)
        V_glyc: specific activity of glycolytic enzymes (µmol substrate·min⁻¹·mg protein⁻¹)
        Phi: max proteome fraction allocatable to ATP-producing enzymes (dimensionless)
    """

    organism: Organism
    gamma_resp: ParameterEntry
    gamma_glyc: ParameterEntry
    V_resp: ParameterEntry
    V_glyc: ParameterEntry
    Phi: ParameterEntry

    @property
    def Vgamma_resp(self) -> float:
        """ATP rate capacity for respiration: V_resp * gamma_resp."""
        return self.V_resp.value * self.gamma_resp.value

    @property
    def Vgamma_glyc(self) -> float:
        """ATP rate capacity for glycolysis: V_glyc * gamma_glyc."""
        return self.V_glyc.value * self.gamma_glyc.value

    @property
    def rho(self) -> float:
        """Capacity ratio ρ = (V·γ)_resp / (V·γ)_glyc. >1 means respiration faster per protein."""
        return self.Vgamma_resp / self.Vgamma_glyc


@dataclass
class ModelParams:
    """Unified input for both Model A and Model B solvers.

    This is the shared interface: both solve_shen() and solve_kukurugya()
    accept this and return ModelResult.
    """

    gamma_resp: float  # mol ATP / mol glucose (respiration)
    gamma_glyc: float  # mol ATP / mol glucose (glycolysis/fermentation)
    V_resp: float      # specific activity, respiratory enzymes
    V_glyc: float      # specific activity, glycolytic enzymes
    Phi: float         # max proteome fraction for ATP enzymes
    g_avail: float     # glucose availability (uptake capacity)

    @classmethod
    def from_organism_params(cls, op: OrganismParams, g_avail: float) -> ModelParams:
        return cls(
            gamma_resp=op.gamma_resp.value,
            gamma_glyc=op.gamma_glyc.value,
            V_resp=op.V_resp.value,
            V_glyc=op.V_glyc.value,
            Phi=op.Phi.value,
            g_avail=g_avail,
        )


# ---------------------------------------------------------------------------
# Parameter database
# ---------------------------------------------------------------------------
# NOTE: Values below are from Kukurugya et al. 2024 PNAS, Table 1 & SI.
# The spec mandates Day-0 verification from original SI — values marked
# PROVISIONAL must be confirmed before audit runs are considered valid.
#
# Kukurugya defines:
#   gamma = ATP yield per glucose (mol/mol)
#   V = maximal specific activity of pathway enzymes (µmol·min⁻¹·mg⁻¹)
#   The "proteome efficiency" = V * gamma
#
# Shen 2024 reports realized efficiency = actual ATP flux / measured enzyme mass.
# We encode Kukurugya's capacity parameters as the primary set because they
# define the 5-parameter model; Shen's realized measurements enter via E2.
# ---------------------------------------------------------------------------

# Kukurugya 2024 PNAS Table 1 / SI Table S1
# gamma values: standard biochemistry textbook + paper-specific pathway definitions
# V values: from compiled enzyme kinetic data (SI Table S1)
# Phi: from proteomics data, fraction of proteome in ATP-generating enzymes

KUKURUGYA_ECOLI = OrganismParams(
    organism="ecoli",
    gamma_resp=ParameterEntry(
        value=26.0,
        unit="mol ATP / mol glucose",
        source_paper="kukurugya2024",
        source_location="Table 1, SI Text S1",
        organism="ecoli",
        note="Full oxidative phosphorylation; P/O ratio ~2.5 assumed for E. coli",
    ),
    gamma_glyc=ParameterEntry(
        value=2.0,
        unit="mol ATP / mol glucose",
        source_paper="kukurugya2024",
        source_location="Table 1",
        organism="ecoli",
        note="Substrate-level phosphorylation only (glycolysis to acetate via Pta-AckA)",
    ),
    V_resp=ParameterEntry(
        value=0.0049,
        unit="µmol glucose · min⁻¹ · mg protein⁻¹",
        source_paper="kukurugya2024",
        source_location="SI Table S1",
        organism="ecoli",
        note="PROVISIONAL — rate-limiting step of respiratory chain per total resp. protein mass",
    ),
    V_glyc=ParameterEntry(
        value=0.23,
        unit="µmol glucose · min⁻¹ · mg protein⁻¹",
        source_paper="kukurugya2024",
        source_location="SI Table S1",
        organism="ecoli",
        note="PROVISIONAL — glycolytic pathway specific activity",
    ),
    Phi=ParameterEntry(
        value=0.26,
        unit="dimensionless",
        source_paper="kukurugya2024",
        source_location="SI Table S2, proteomics",
        organism="ecoli",
        note="PROVISIONAL — fraction of proteome in ATP-generating enzymes",
    ),
)

KUKURUGYA_YEAST = OrganismParams(
    organism="yeast",
    gamma_resp=ParameterEntry(
        value=26.0,
        unit="mol ATP / mol glucose",
        source_paper="kukurugya2024",
        source_location="Table 1",
        organism="yeast",
        note="Full respiration; P/O ~2.5 for yeast mitochondria",
    ),
    gamma_glyc=ParameterEntry(
        value=2.0,
        unit="mol ATP / mol glucose",
        source_paper="kukurugya2024",
        source_location="Table 1",
        organism="yeast",
        note="Fermentation to ethanol: net 2 ATP per glucose",
    ),
    V_resp=ParameterEntry(
        value=0.0085,
        unit="µmol glucose · min⁻¹ · mg protein⁻¹",
        source_paper="kukurugya2024",
        source_location="SI Table S1",
        organism="yeast",
        note="PROVISIONAL — respiratory pathway specific activity",
    ),
    V_glyc=ParameterEntry(
        value=0.46,
        unit="µmol glucose · min⁻¹ · mg protein⁻¹",
        source_paper="kukurugya2024",
        source_location="SI Table S1",
        organism="yeast",
        note="PROVISIONAL — glycolytic/fermentative specific activity",
    ),
    Phi=ParameterEntry(
        value=0.20,
        unit="dimensionless",
        source_paper="kukurugya2024",
        source_location="SI Table S2",
        organism="yeast",
        note="PROVISIONAL",
    ),
)

KUKURUGYA_MAMMALIAN = OrganismParams(
    organism="mammalian",
    gamma_resp=ParameterEntry(
        value=32.0,
        unit="mol ATP / mol glucose",
        source_paper="kukurugya2024",
        source_location="Table 1",
        organism="mammalian",
        note="Full oxidative phosphorylation; P/O ~2.5, malate-aspartate shuttle",
    ),
    gamma_glyc=ParameterEntry(
        value=2.0,
        unit="mol ATP / mol glucose",
        source_paper="kukurugya2024",
        source_location="Table 1",
        organism="mammalian",
        note="Aerobic glycolysis to lactate: 2 ATP net",
    ),
    V_resp=ParameterEntry(
        value=0.0040,
        unit="µmol glucose · min⁻¹ · mg protein⁻¹",
        source_paper="kukurugya2024",
        source_location="SI Table S1",
        organism="mammalian",
        note="PROVISIONAL — respiratory chain specific activity for mammalian cells",
    ),
    V_glyc=ParameterEntry(
        value=0.39,
        unit="µmol glucose · min⁻¹ · mg protein⁻¹",
        source_paper="kukurugya2024",
        source_location="SI Table S1",
        organism="mammalian",
        note="PROVISIONAL — glycolytic specific activity",
    ),
    Phi=ParameterEntry(
        value=0.18,
        unit="dimensionless",
        source_paper="kukurugya2024",
        source_location="SI Table S2",
        organism="mammalian",
        note="PROVISIONAL",
    ),
)

ALL_ORGANISMS: dict[Organism, OrganismParams] = {
    "ecoli": KUKURUGYA_ECOLI,
    "yeast": KUKURUGYA_YEAST,
    "mammalian": KUKURUGYA_MAMMALIAN,
}


# ---------------------------------------------------------------------------
# Uncertainty ranges for sensitivity analysis (E3)
# Ranges represent ±bounds from SI or literature variation.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ParamRange:
    """Min/max range for a parameter, for Sobol/Morris sampling."""

    name: str
    low: float
    high: float
    unit: str
    source: str


def get_param_ranges(organism: Organism) -> list[ParamRange]:
    """Return uncertainty ranges for E3 sensitivity analysis.

    Ranges are ±30% of central values (PROVISIONAL — should be narrowed
    to actual SI uncertainty bounds once verified at Day-0 gate).
    """
    op = ALL_ORGANISMS[organism]
    factor_lo, factor_hi = 0.7, 1.3

    return [
        ParamRange("gamma_resp", op.gamma_resp.value * factor_lo,
                   op.gamma_resp.value * factor_hi, op.gamma_resp.unit,
                   "±30% of Kukurugya Table 1"),
        ParamRange("gamma_glyc", op.gamma_glyc.value * factor_lo,
                   op.gamma_glyc.value * factor_hi, op.gamma_glyc.unit,
                   "±30% of Kukurugya Table 1"),
        ParamRange("V_resp", op.V_resp.value * factor_lo,
                   op.V_resp.value * factor_hi, op.V_resp.unit,
                   "±30% of Kukurugya SI Table S1"),
        ParamRange("V_glyc", op.V_glyc.value * factor_lo,
                   op.V_glyc.value * factor_hi, op.V_glyc.unit,
                   "±30% of Kukurugya SI Table S1"),
        ParamRange("Phi", max(0.05, op.Phi.value * factor_lo),
                   min(0.50, op.Phi.value * factor_hi), op.Phi.unit,
                   "±30% of Kukurugya SI Table S2"),
    ]


# ---------------------------------------------------------------------------
# CLI: provenance table
# ---------------------------------------------------------------------------

def provenance_table() -> pd.DataFrame:
    """Build a tidy provenance table of all parameters."""
    rows = []
    for org_name, op in ALL_ORGANISMS.items():
        for param_name in ("gamma_resp", "gamma_glyc", "V_resp", "V_glyc", "Phi"):
            entry: ParameterEntry = getattr(op, param_name)
            rows.append({
                "organism": org_name,
                "parameter": param_name,
                "value": entry.value,
                "unit": entry.unit,
                "source_paper": entry.source_paper,
                "source_location": entry.source_location,
                "note": entry.note,
            })
    return pd.DataFrame(rows)


def derived_table() -> pd.DataFrame:
    """Derived quantities: V·γ and ρ for each organism."""
    rows = []
    for org_name, op in ALL_ORGANISMS.items():
        rows.append({
            "organism": org_name,
            "Vgamma_resp": op.Vgamma_resp,
            "Vgamma_glyc": op.Vgamma_glyc,
            "rho (resp/glyc)": op.rho,
            "verdict_capacity": "respiration" if op.rho > 1 else "glycolysis",
        })
    return pd.DataFrame(rows)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Parameter provenance for Warburg audit")
    parser.add_argument("--table", action="store_true", help="Print provenance table")
    parser.add_argument("--derived", action="store_true", help="Print derived V·γ and ρ")
    parser.add_argument("--csv", type=str, default=None, help="Export provenance to CSV")
    args = parser.parse_args()

    if args.table or (not args.derived and not args.csv):
        print("=" * 80)
        print("PARAMETER PROVENANCE TABLE")
        print("=" * 80)
        df = provenance_table()
        print(df.to_string(index=False))
        print()

    if args.derived:
        print("=" * 80)
        print("DERIVED QUANTITIES (V·γ and ρ)")
        print("=" * 80)
        df = derived_table()
        print(df.to_string(index=False))
        print()

    if args.csv:
        df = provenance_table()
        df.to_csv(args.csv, index=False)
        print(f"Provenance table exported to {args.csv}")


if __name__ == "__main__":
    main()
