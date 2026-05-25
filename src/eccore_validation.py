"""M6 — ecCore enzyme-constrained FBA validation (STRETCH).

Loads the published E. coli core model via COBRApy, adds a single
enzyme-mass constraint (Σ mᵢ·|vᵢ| ≤ Φ), sweeps Φ and glucose uptake,
and confirms the respiration→fermentation transition predicted by the
toy models (Model A & B).

The key observable: acetate overflow onset as proteome budget (Φ) tightens
or glucose uptake increases — qualitatively matching the audit direction.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

try:
    import cobra
    from cobra.io import load_model
    HAS_COBRA = True
except ImportError:
    HAS_COBRA = False

RESULTS_DIR = Path("results")


# Approximate molecular weights (kDa) for E. coli core model enzymes.
# These are rough estimates for the enzyme-mass constraint demonstration.
# In a full GECKO model, these come from UniProt; here we use uniform weights
# to show the *qualitative* effect of proteome limitation.
DEFAULT_ENZYME_WEIGHT = 50.0  # kDa, uniform approximation


def _get_metabolic_reaction_ids(model: "cobra.Model") -> list[str]:
    """Get IDs of metabolic reactions (exclude exchange, demand, biomass)."""
    excluded_prefixes = ("EX_", "DM_", "SK_")
    excluded_ids = {"BIOMASS_Ecoli_core_w_GAM", "ATPM"}
    rxn_ids = []
    for rxn in model.reactions:
        if rxn.id in excluded_ids:
            continue
        if any(rxn.id.startswith(p) for p in excluded_prefixes):
            continue
        rxn_ids.append(rxn.id)
    return rxn_ids


def run_eccore_phi_sweep(
    phi_values: np.ndarray | None = None,
    glucose_uptake: float = 10.0,
    n_phi: int = 30,
) -> pd.DataFrame:
    """Sweep proteome budget Φ and observe metabolic phenotype.

    Adds constraint: Σ (enzyme_weight_i * |v_i|) ≤ Φ_total
    for all metabolic reactions.

    Returns DataFrame with Φ, growth rate, acetate secretion, O2 uptake, etc.
    """
    if not HAS_COBRA:
        raise ImportError("COBRApy is required for M6. Install via: pip install cobra")

    model = load_model("textbook")  # e_coli_core

    # Set glucose uptake
    model.reactions.get_by_id("EX_glc__D_e").lower_bound = -glucose_uptake

    # Get metabolic reactions for enzyme constraint
    met_rxn_ids = _get_metabolic_reaction_ids(model)

    if phi_values is None:
        # Range: from very tight (forces fermentation) to loose (unconstrained)
        # Units: kDa · mmol/gDW/h (enzyme_weight × flux)
        # With ~70 reactions at weight 50 kDa and fluxes ~0-20, max ~70000
        phi_values = np.linspace(500, 50000, n_phi)

    rows = []
    for phi in phi_values:
        with model:
            # Add enzyme capacity constraint using indicator constraints
            # For reversible reactions, constrain the absolute flux via two constraints
            enzyme_constraint_fwd = model.problem.Constraint(
                sum(
                    DEFAULT_ENZYME_WEIGHT * model.reactions.get_by_id(rid).flux_expression
                    for rid in met_rxn_ids
                    if not model.reactions.get_by_id(rid).reversibility
                ),
                ub=phi,
                name="enzyme_cap_irrev",
            )

            # Simpler approach: constrain sum of forward fluxes only
            # (this is an approximation but captures the qualitative effect)
            constraint_expr = sum(
                DEFAULT_ENZYME_WEIGHT * model.reactions.get_by_id(rid).forward_variable
                for rid in met_rxn_ids
            )
            enzyme_constraint = model.problem.Constraint(
                constraint_expr,
                ub=phi,
                name="enzyme_capacity",
            )
            model.add_cons_vars(enzyme_constraint)

            sol = model.optimize()

            if sol.status == "optimal":
                # Key fluxes
                growth = sol.objective_value
                acetate_flux = sol.fluxes.get("EX_ac_e", 0.0)
                o2_flux = sol.fluxes.get("EX_o2_e", 0.0)  # negative = uptake
                ethanol_flux = sol.fluxes.get("EX_etoh_e", 0.0)
                co2_flux = sol.fluxes.get("EX_co2_e", 0.0)

                # Compute fermentative fraction (acetate + ethanol secretion vs total carbon out)
                ferm_carbon = max(0, acetate_flux) + max(0, ethanol_flux)
                total_carbon_out = ferm_carbon + max(0, co2_flux)
                frac_ferm = ferm_carbon / total_carbon_out if total_carbon_out > 1e-6 else 0.0

                rows.append({
                    "Phi": float(phi),
                    "growth_rate": growth,
                    "acetate_secretion": acetate_flux,
                    "o2_uptake": -o2_flux,  # make positive
                    "ethanol_secretion": ethanol_flux,
                    "co2_secretion": co2_flux,
                    "frac_fermentative": frac_ferm,
                    "glucose_uptake": glucose_uptake,
                    "status": "optimal",
                })
            else:
                rows.append({
                    "Phi": float(phi),
                    "growth_rate": 0.0,
                    "acetate_secretion": 0.0,
                    "o2_uptake": 0.0,
                    "ethanol_secretion": 0.0,
                    "co2_secretion": 0.0,
                    "frac_fermentative": 0.0,
                    "glucose_uptake": glucose_uptake,
                    "status": sol.status,
                })

    return pd.DataFrame(rows)


def run_eccore_glucose_sweep(
    glucose_values: np.ndarray | None = None,
    phi_tight: float = 5000.0,
    n_glc: int = 30,
) -> pd.DataFrame:
    """Sweep glucose uptake at fixed tight proteome budget.

    This mirrors the toy-model glucose sweep (E1) but in a real GEM.
    """
    if not HAS_COBRA:
        raise ImportError("COBRApy is required for M6. Install via: pip install cobra")

    model = load_model("textbook")

    met_rxn_ids = _get_metabolic_reaction_ids(model)

    if glucose_values is None:
        glucose_values = np.linspace(1.0, 20.0, n_glc)

    rows = []
    for glc in glucose_values:
        with model:
            model.reactions.get_by_id("EX_glc__D_e").lower_bound = -float(glc)

            constraint_expr = sum(
                DEFAULT_ENZYME_WEIGHT * model.reactions.get_by_id(rid).forward_variable
                for rid in met_rxn_ids
            )
            enzyme_constraint = model.problem.Constraint(
                constraint_expr,
                ub=phi_tight,
                name="enzyme_capacity",
            )
            model.add_cons_vars(enzyme_constraint)

            sol = model.optimize()

            if sol.status == "optimal":
                growth = sol.objective_value
                acetate_flux = sol.fluxes.get("EX_ac_e", 0.0)
                o2_flux = sol.fluxes.get("EX_o2_e", 0.0)
                ethanol_flux = sol.fluxes.get("EX_etoh_e", 0.0)
                co2_flux = sol.fluxes.get("EX_co2_e", 0.0)

                ferm_carbon = max(0, acetate_flux) + max(0, ethanol_flux)
                total_carbon_out = ferm_carbon + max(0, co2_flux)
                frac_ferm = ferm_carbon / total_carbon_out if total_carbon_out > 1e-6 else 0.0

                rows.append({
                    "glucose_uptake": float(glc),
                    "Phi": phi_tight,
                    "growth_rate": growth,
                    "acetate_secretion": acetate_flux,
                    "o2_uptake": -o2_flux,
                    "ethanol_secretion": ethanol_flux,
                    "co2_secretion": co2_flux,
                    "frac_fermentative": frac_ferm,
                    "status": "optimal",
                })
            else:
                rows.append({
                    "glucose_uptake": float(glc),
                    "Phi": phi_tight,
                    "growth_rate": 0.0,
                    "acetate_secretion": 0.0,
                    "o2_uptake": 0.0,
                    "ethanol_secretion": 0.0,
                    "co2_secretion": 0.0,
                    "frac_fermentative": 0.0,
                    "status": sol.status,
                })

    return pd.DataFrame(rows)


def run_m6() -> dict[str, pd.DataFrame]:
    """Run full M6 validation: Φ-sweep + glucose-sweep."""
    RESULTS_DIR.mkdir(exist_ok=True)

    print("[M6] Running ecCore Φ-sweep (enzyme-constrained FBA)...")
    df_phi = run_eccore_phi_sweep()
    phi_path = RESULTS_DIR / "m6_phi_sweep.csv"
    df_phi.to_csv(phi_path, index=False)
    print(f"  → {phi_path} ({len(df_phi)} rows)")

    # Check for acetate overflow onset
    overflow_rows = df_phi[df_phi["acetate_secretion"] > 0.1]
    if not overflow_rows.empty:
        onset_phi = overflow_rows["Phi"].max()  # tightest Phi where overflow appears
        print(f"  Acetate overflow onset at Φ ≤ {onset_phi:.0f} kDa·mmol/gDW/h")
    else:
        print("  No acetate overflow detected in Φ range")

    print("[M6] Running ecCore glucose-sweep at tight Φ...")
    df_glc = run_eccore_glucose_sweep()
    glc_path = RESULTS_DIR / "m6_glucose_sweep.csv"
    df_glc.to_csv(glc_path, index=False)
    print(f"  → {glc_path} ({len(df_glc)} rows)")

    # Check for fermentation onset with rising glucose
    ferm_rows = df_glc[df_glc["frac_fermentative"] > 0.05]
    if not ferm_rows.empty:
        onset_glc = ferm_rows["glucose_uptake"].min()
        print(f"  Fermentation onset at glucose uptake ≥ {onset_glc:.1f} mmol/gDW/h")
    else:
        print("  No fermentation onset detected in glucose range")

    return {"phi_sweep": df_phi, "glucose_sweep": df_glc}


if __name__ == "__main__":
    run_m6()
