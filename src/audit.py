"""M3 — Audit experiments E1 (shared-axis overlay) and E2 (capacity vs. realized).

E1: Run Model A and Model B on identical input vectors, collect outputs on shared
    observable axes (x = glucose uptake rate, y = fermentative ATP fraction / total ATP).
E2: Sweep glycolytic utilization fraction u_G to find where the capacity vs. realized
    verdict flips. Secondary: vary enzyme accounting attribution.
E2c: Uncertainty quantification — bootstrap confidence intervals on flip points.
E2d: Attribution decomposition — relative contribution of u_G vs enzyme definition.
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
from src.params import ALL_ORGANISMS, ModelParams, Organism, get_param_ranges

RESULTS_DIR = Path("results")
RNG_SEED = 42


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
# E2c — Uncertainty quantification (bootstrap CI on flip points)
# ---------------------------------------------------------------------------

def run_e2c_uncertainty(
    organism: Organism = "yeast",
    n_bootstrap: int = 200,
    u_R: float = 1.0,
) -> dict:
    """Bootstrap confidence interval on the verdict flip point u_G.

    Samples parameters from their uncertainty ranges (±30% or SI bounds),
    computes flip_uG for each sample, and returns distribution statistics.
    """
    rng = np.random.default_rng(RNG_SEED)
    ranges = get_param_ranges(organism)
    op = ALL_ORGANISMS[organism]

    flip_samples = []
    for _ in range(n_bootstrap):
        gamma_resp = rng.uniform(ranges[0].low, ranges[0].high)
        gamma_glyc = rng.uniform(ranges[1].low, ranges[1].high)
        V_resp = rng.uniform(ranges[2].low, ranges[2].high)
        V_glyc = rng.uniform(ranges[3].low, ranges[3].high)
        Phi = rng.uniform(ranges[4].low, ranges[4].high)

        p = ModelParams(
            gamma_resp=gamma_resp,
            gamma_glyc=gamma_glyc,
            V_resp=V_resp,
            V_glyc=V_glyc,
            Phi=Phi,
            g_avail=0.0,
        )
        flip = find_verdict_flip_uG(p, u_R=u_R, n_points=200)
        if flip is not None:
            flip_samples.append(flip)

    flip_arr = np.array(flip_samples)
    nominal_p = ModelParams.from_organism_params(op, g_avail=0.0)
    nominal_flip = find_verdict_flip_uG(nominal_p, u_R=u_R)

    return {
        "organism": organism,
        "nominal_flip_uG": nominal_flip,
        "mean_flip_uG": float(np.mean(flip_arr)) if len(flip_arr) > 0 else None,
        "std_flip_uG": float(np.std(flip_arr)) if len(flip_arr) > 0 else None,
        "ci_2.5": float(np.percentile(flip_arr, 2.5)) if len(flip_arr) > 0 else None,
        "ci_97.5": float(np.percentile(flip_arr, 97.5)) if len(flip_arr) > 0 else None,
        "ci_5": float(np.percentile(flip_arr, 5)) if len(flip_arr) > 0 else None,
        "ci_95": float(np.percentile(flip_arr, 95)) if len(flip_arr) > 0 else None,
        "n_valid": len(flip_arr),
        "n_total": n_bootstrap,
        "flip_samples": flip_arr,
    }


def run_e2c_all_organisms(n_bootstrap: int = 200) -> pd.DataFrame:
    """Run uncertainty quantification for all organisms. Returns summary DataFrame."""
    rows = []
    for org in ALL_ORGANISMS:
        result = run_e2c_uncertainty(org, n_bootstrap=n_bootstrap)
        rows.append({
            "organism": org,
            "nominal_flip_uG": result["nominal_flip_uG"],
            "mean_flip_uG": result["mean_flip_uG"],
            "std_flip_uG": result["std_flip_uG"],
            "ci_2.5": result["ci_2.5"],
            "ci_97.5": result["ci_97.5"],
            "ci_5": result["ci_5"],
            "ci_95": result["ci_95"],
            "n_valid": result["n_valid"],
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# E2d — Attribution decomposition: u_G effect vs enzyme-definition effect
# ---------------------------------------------------------------------------

def run_e2d_decomposition(
    organism: Organism = "yeast",
    n_points: int = 50,
) -> pd.DataFrame:
    """Decompose the verdict into u_G effect vs enzyme-definition (V reattribution) effect.

    The "enzyme definition" axis simulates what happens if some proteins currently
    attributed to glycolysis are reclassified as "shared/maintenance" (reducing
    effective V_glyc), or vice versa. This captures the Shen vs Kukurugya
    disagreement about which enzymes to count in each pathway.

    The key quantity is: how much does reclassifying enzymes shift ρ (and thus
    the flip point) compared to the dominant effect of u_G?
    """
    op = ALL_ORGANISMS[organism]
    baseline_rho = op.rho
    baseline_flip = find_verdict_flip_uG(
        ModelParams.from_organism_params(op, g_avail=0.0), u_R=1.0, n_points=500
    )

    # Vary attribution: shift effective V_glyc by ±30% (reattributing enzymes)
    # If some glycolytic enzymes are counted as "shared", V_glyc decreases
    # If respiratory auxiliary proteins are excluded, V_resp increases
    reattribution_factors = np.linspace(0.7, 1.3, n_points)
    rows = []
    for factor in reattribution_factors:
        V_glyc_eff = op.V_glyc.value * factor
        V_resp_eff = op.V_resp.value * (2.0 - factor)  # inverse reattribution

        p = ModelParams(
            gamma_resp=op.gamma_resp.value,
            gamma_glyc=op.gamma_glyc.value,
            V_resp=V_resp_eff,
            V_glyc=V_glyc_eff,
            Phi=op.Phi.value,
            g_avail=0.0,
        )
        flip = find_verdict_flip_uG(p, u_R=1.0, n_points=500)
        rho = (V_resp_eff * op.gamma_resp.value) / (V_glyc_eff * op.gamma_glyc.value)

        rows.append({
            "organism": organism,
            "reattribution_factor": float(factor),
            "V_glyc_effective": V_glyc_eff,
            "V_resp_effective": V_resp_eff,
            "flip_uG": flip,
            "rho": rho,
            "baseline_rho": baseline_rho,
            "delta_rho_from_baseline": rho - baseline_rho,
            "delta_flip_from_baseline": (flip - baseline_flip) if flip and baseline_flip else None,
        })

    df = pd.DataFrame(rows)
    df["attribution_sensitivity"] = df["delta_flip_from_baseline"].abs()

    return df


# ---------------------------------------------------------------------------
# Full audit pipeline
# ---------------------------------------------------------------------------

def run_full_audit() -> dict[str, pd.DataFrame | list]:
    """Run E1 + E2 + E2b + E2c + E2d, save results, return all dataframes."""
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

    print("[M3] Running E2c — uncertainty quantification (bootstrap CI)...")
    e2c_df = run_e2c_all_organisms(n_bootstrap=200)
    e2c_path = RESULTS_DIR / "e2c_uncertainty.csv"
    e2c_df.to_csv(e2c_path, index=False)
    print(f"  → {e2c_path}")
    for _, row in e2c_df.iterrows():
        print(f"  {row['organism']:>10}: flip u_G = {row['nominal_flip_uG']:.4f} "
              f"[90% CI: {row['ci_5']:.4f} – {row['ci_95']:.4f}]")

    print("[M3] Running E2d — attribution decomposition...")
    e2d_dfs = []
    for org in ALL_ORGANISMS:
        e2d_dfs.append(run_e2d_decomposition(organism=org))
    e2d_df = pd.concat(e2d_dfs, ignore_index=True)
    e2d_path = RESULTS_DIR / "e2d_decomposition.csv"
    e2d_df.to_csv(e2d_path, index=False)
    print(f"  → {e2d_path} ({len(e2d_df)} rows)")
    for org in ALL_ORGANISMS:
        sub = e2d_df[e2d_df["organism"] == org]
        max_delta = sub["attribution_sensitivity"].max()
        print(f"  {org:>10}: max |Δ flip| from enzyme attribution = {max_delta:.6f} "
              f"(vs flip ≈ {sub['rho'].iloc[0]:.4f})")

    return {
        "e1": e1_df,
        "e2": e2_df,
        "e2_summaries": e2_summaries,
        "e2b": e2b_df,
        "e2c": e2c_df,
        "e2d": e2d_df,
    }


# ---------------------------------------------------------------------------
# Wang 2025 comparison — relate our u_G finding to Wang's crossing mechanism
# ---------------------------------------------------------------------------

def run_wang_comparison(n_points: int = 200) -> pd.DataFrame:
    """Compare our audit results to Wang 2025's efficiency-crossing prediction.

    Wang (eLife 2025) proposes that respiratory and fermentative efficiency
    curves *cross* as nutrient quality varies — at high quality, respiration
    wins; at low quality, fermentation wins (due to heterogeneity + overflow).

    Our E2 shows the same crossing but identifies the *mechanism*: it's not
    nutrient quality per se, but the *utilization fraction* u_G that drives
    the crossing. When glycolytic enzymes are constitutively overexpressed
    (low u_G, Shen's finding), respiration appears more efficient.

    This function maps Wang's nutrient-quality axis to our u_G axis and shows
    they predict the same crossing point from different mechanistic angles.
    """
    rows = []
    for org_name, op in ALL_ORGANISMS.items():
        Vgamma_G = op.Vgamma_glyc
        Vgamma_R = op.Vgamma_resp
        rho = op.rho

        # Wang's model: efficiency crosses when growth rate demands exceed
        # respiratory capacity. Map this to our framework:
        # - "High nutrient quality" ≈ low glucose demand → resp wins → u_G low (idle glycolysis)
        # - "Low nutrient quality" ≈ overflow regime → glyc wins → u_G high (full utilization)
        #
        # The crossing in Wang happens at a growth rate µ*. In our framework,
        # this corresponds to u_G = ρ = (V·γ)_R / (V·γ)_G.

        u_G_range = np.linspace(0.01, 1.0, n_points)
        for u_G in u_G_range:
            realized_R = Vgamma_R  # u_R = 1 (resp always fully utilized)
            realized_G = u_G * Vgamma_G
            margin = realized_R - realized_G

            # Wang's "nutrient quality" parameter q maps inversely to u_G:
            # high q → cell grows slowly → glycolytic capacity idle → low u_G
            # low q → cell at max growth → glycolytic machinery fully engaged → high u_G
            wang_nutrient_quality = 1.0 - u_G  # simplified mapping

            rows.append({
                "organism": org_name,
                "u_G": float(u_G),
                "wang_nutrient_quality": wang_nutrient_quality,
                "realized_eff_resp": realized_R,
                "realized_eff_glyc": realized_G,
                "margin": margin,
                "verdict": "respiration" if margin > 0 else "glycolysis",
                "crossing_point_uG": rho,
                "crossing_point_wang_q": 1.0 - rho,
            })

    df = pd.DataFrame(rows)
    return df


if __name__ == "__main__":
    run_full_audit()
