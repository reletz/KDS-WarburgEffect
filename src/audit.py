"""M3 — Audit experiments E1 (shared-axis overlay) and E2 (capacity vs. realized).

E1: Run Model A and Model B on identical input vectors, collect outputs on shared
    observable axes (x = glucose uptake rate, y = fermentative ATP fraction / total ATP).
E2: Sweep glycolytic utilization fraction u_G to find where the capacity vs. realized
    verdict flips. Secondary: vary enzyme accounting attribution.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src.models import (
    ModelResult,
    sweep_glucose,
    compute_efficiency_with_utilization,
    find_verdict_flip_uG,
)
from src.params import ALL_ORGANISMS, ModelParams, Organism

RESULTS_DIR = Path("results")


def _results_to_df(results: list[ModelResult], model_name: str) -> pd.DataFrame:
    rows = []
    for r in results:
        rows.append({
            "model": model_name,
            "g_avail": r.g_avail,
            "J_ATP": r.J_ATP,
            "J_ATP_glyc": r.J_ATP_glyc,
            "J_ATP_resp": r.J_ATP_resp,
            "frac_glyc": r.frac_glyc,
            "phi_G": r.phi_G,
            "phi_R": r.phi_R,
            "phenotype": r.phenotype,
            "eps_ratio": r.eps_ratio,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# E1 — Shared-axis overlay
# ---------------------------------------------------------------------------

def run_e1(organism: Organism = "yeast", n_points: int = 300) -> pd.DataFrame:
    """Run both models on identical glucose grid, return combined dataframe."""
    op = ALL_ORGANISMS[organism]
    base = ModelParams.from_organism_params(op, g_avail=0.0)

    g_max = 2.5 * op.V_glyc.value * op.Phi.value
    g_range = np.linspace(1e-6, g_max, n_points)

    results_a = sweep_glucose(base, g_range=g_range, model="shen")
    results_b = sweep_glucose(base, g_range=g_range, model="kukurugya")

    df_a = _results_to_df(results_a, "Model A (Shen LP)")
    df_b = _results_to_df(results_b, "Model B (Kukurugya)")
    df = pd.concat([df_a, df_b], ignore_index=True)
    df["organism"] = organism
    return df


def run_e1_all_organisms(n_points: int = 300) -> pd.DataFrame:
    """Run E1 for all three organisms."""
    dfs = []
    for org in ALL_ORGANISMS:
        dfs.append(run_e1(organism=org, n_points=n_points))
    return pd.concat(dfs, ignore_index=True)


# ---------------------------------------------------------------------------
# E2 — Capacity vs. Realized (sweep u_G)
# ---------------------------------------------------------------------------

def run_e2(
    organism: Organism = "yeast",
    u_R: float = 1.0,
    n_points: int = 500,
) -> tuple[pd.DataFrame, dict]:
    """Sweep glycolytic utilization u_G, compute capacity & realized efficiency.

    Returns:
        (dataframe of sweep results, summary dict with flip point info)
    """
    op = ALL_ORGANISMS[organism]
    base = ModelParams.from_organism_params(op, g_avail=0.0)

    u_G_values = np.linspace(1.0, 0.01, n_points)
    rows = []
    for u_G in u_G_values:
        r = compute_efficiency_with_utilization(base, u_G=float(u_G), u_R=u_R)
        r["organism"] = organism
        rows.append(r)
    df = pd.DataFrame(rows)

    flip_uG = find_verdict_flip_uG(base, u_R=u_R, n_points=n_points)

    summary = {
        "organism": organism,
        "u_R": u_R,
        "flip_uG": flip_uG,
        "Vgamma_G_capacity": base.V_glyc * base.gamma_glyc,
        "Vgamma_R_capacity": base.V_resp * base.gamma_resp,
        "rho_capacity": (base.V_resp * base.gamma_resp) / (base.V_glyc * base.gamma_glyc),
        "verdict_at_full_capacity": "glycolysis" if base.V_glyc * base.gamma_glyc > base.V_resp * base.gamma_resp else "respiration",
    }
    if flip_uG is not None:
        r_at_flip = compute_efficiency_with_utilization(base, u_G=flip_uG, u_R=u_R)
        summary["margin_at_flip"] = r_at_flip["margin"]
        summary["verdict_below_flip"] = "respiration"

    return df, summary


def run_e2_all_organisms(u_R: float = 1.0) -> tuple[pd.DataFrame, list[dict]]:
    """Run E2 for all organisms."""
    dfs = []
    summaries = []
    for org in ALL_ORGANISMS:
        df, summary = run_e2(organism=org, u_R=u_R)
        dfs.append(df)
        summaries.append(summary)
    return pd.concat(dfs, ignore_index=True), summaries


# ---------------------------------------------------------------------------
# E2b — Secondary: enzyme accounting attribution variation
# ---------------------------------------------------------------------------

def run_e2b_accounting(
    organism: Organism = "yeast",
    membrane_fractions: np.ndarray | None = None,
) -> pd.DataFrame:
    """Vary how much membrane/maintenance protein is attributed to each sector.

    Simulates the "which enzymes count" axis of the disagreement.
    membrane_frac: fraction of Phi attributed to overhead (neither G nor R).
    As overhead increases, effective Phi for ATP enzymes shrinks.
    """
    op = ALL_ORGANISMS[organism]
    if membrane_fractions is None:
        membrane_fractions = np.linspace(0.0, 0.5, 50)

    rows = []
    for mf in membrane_fractions:
        effective_Phi = op.Phi.value * (1.0 - mf)
        if effective_Phi < 1e-6:
            continue
        p = ModelParams(
            gamma_resp=op.gamma_resp.value,
            gamma_glyc=op.gamma_glyc.value,
            V_resp=op.V_resp.value,
            V_glyc=op.V_glyc.value,
            Phi=effective_Phi,
            g_avail=0.0,
        )
        flip = find_verdict_flip_uG(p, u_R=1.0, n_points=500)
        rows.append({
            "organism": organism,
            "membrane_frac": float(mf),
            "effective_Phi": effective_Phi,
            "flip_uG": flip,
            "Vgamma_R_cap": p.V_resp * p.gamma_resp,
            "Vgamma_G_cap": p.V_glyc * p.gamma_glyc,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Full audit pipeline
# ---------------------------------------------------------------------------

def run_full_audit() -> dict[str, pd.DataFrame | list]:
    """Run E1 + E2 + E2b, save results, return all dataframes."""
    RESULTS_DIR.mkdir(exist_ok=True)

    print("[M3] Running E1 — shared-axis overlay for all organisms...")
    e1_df = run_e1_all_organisms()
    e1_path = RESULTS_DIR / "e1_overlay.csv"
    e1_df.to_csv(e1_path, index=False)
    print(f"  → {e1_path} ({len(e1_df)} rows)")

    print("[M3] Running E2 — capacity vs. realized sweep...")
    e2_df, e2_summaries = run_e2_all_organisms()
    e2_path = RESULTS_DIR / "e2_capacity_vs_realized.csv"
    e2_df.to_csv(e2_path, index=False)
    print(f"  → {e2_path} ({len(e2_df)} rows)")
    for s in e2_summaries:
        flip = s["flip_uG"]
        label = f"u_G = {flip:.4f}" if flip else "NO FLIP"
        print(f"  {s['organism']:>10}: capacity ρ = {s['rho_capacity']:.4f}, "
              f"verdict flips at {label}")

    print("[M3] Running E2b — enzyme accounting variation...")
    e2b_dfs = []
    for org in ALL_ORGANISMS:
        e2b_dfs.append(run_e2b_accounting(organism=org))
    e2b_df = pd.concat(e2b_dfs, ignore_index=True)
    e2b_path = RESULTS_DIR / "e2b_accounting.csv"
    e2b_df.to_csv(e2b_path, index=False)
    print(f"  → {e2b_path} ({len(e2b_df)} rows)")

    return {
        "e1": e1_df,
        "e2": e2_df,
        "e2_summaries": e2_summaries,
        "e2b": e2b_df,
    }


if __name__ == "__main__":
    run_full_audit()
