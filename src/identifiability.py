"""M4 — Identifiability / value-of-information analysis (E3).

Runs Sobol and Morris global sensitivity analysis (via SALib) with CONTINUOUS
output: the efficiency margin m = u_R·(V·γ)_R − u_G·(V·γ)_G, or the critical
glucose threshold g* where the metabolic switch occurs.

Ranks parameters by total-order Sobol index → top parameter = "the measurement
that, if pinned down, most reduces verdict uncertainty" = dispute-settling experiment.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from SALib.analyze import morris as morris_analyze, sobol as sobol_analyze
from SALib.sample import morris as morris_sample, sobol as sobol_sample

from src.models import ModelParams, compute_efficiency_with_utilization, sweep_glucose, solve_shen
from src.params import ALL_ORGANISMS, Organism, get_param_ranges

RESULTS_DIR = Path("results")
RNG_SEED = 42


def _build_salib_problem(organism: Organism, include_uG: bool = True) -> dict:
    """Build SALib problem dict from parameter ranges."""
    ranges = get_param_ranges(organism)
    names = [r.name for r in ranges]
    bounds = [[r.low, r.high] for r in ranges]

    if include_uG:
        names.append("u_G")
        bounds.append([0.1, 1.0])  # utilization fraction range

    return {
        "num_vars": len(names),
        "names": names,
        "bounds": bounds,
    }


def _evaluate_margin(X: np.ndarray, u_R: float = 1.0) -> np.ndarray:
    """Evaluate efficiency margin for each parameter sample row.

    Output: m = u_R·(V·γ)_R − u_G·(V·γ)_G
    Positive → respiration wins on realized basis.
    """
    n = X.shape[0]
    Y = np.zeros(n)

    for i in range(n):
        row = X[i]
        gamma_resp = row[0]
        gamma_glyc = row[1]
        V_resp = row[2]
        V_glyc = row[3]
        _Phi = row[4]  # doesn't affect margin directly but keeps problem structure
        u_G = row[5] if X.shape[1] > 5 else 1.0

        Vgamma_R = V_resp * gamma_resp
        Vgamma_G = V_glyc * gamma_glyc
        margin = u_R * Vgamma_R - u_G * Vgamma_G
        Y[i] = margin

    return Y


def _evaluate_critical_glucose(X: np.ndarray) -> np.ndarray:
    """Evaluate critical glucose g* where metabolic switch occurs.

    g* = V_R * Phi (the point where respiration alone saturates the proteome
    and overflow begins). This is a continuous, interpretable output.
    """
    n = X.shape[0]
    Y = np.zeros(n)

    for i in range(n):
        V_resp = X[i, 2]
        Phi = X[i, 4]
        Y[i] = V_resp * Phi

    return Y


# ---------------------------------------------------------------------------
# Sobol analysis
# ---------------------------------------------------------------------------

def run_sobol(
    organism: Organism = "yeast",
    N: int = 1024,
    output: str = "margin",
) -> tuple[pd.DataFrame, dict]:
    """Run Sobol sensitivity analysis.

    Args:
        organism: Which organism's parameter ranges to use.
        N: Base sample size (total evaluations = N * (2D + 2)).
        output: "margin" or "g_star" — which continuous output to analyze.

    Returns:
        (DataFrame of Sobol indices, problem dict for reference)
    """
    include_uG = (output == "margin")
    problem = _build_salib_problem(organism, include_uG=include_uG)

    X = sobol_sample.sample(problem, N, calc_second_order=False, seed=RNG_SEED)

    if output == "margin":
        Y = _evaluate_margin(X)
    else:
        Y = _evaluate_critical_glucose(X)

    Si = sobol_analyze.analyze(problem, Y, calc_second_order=False, seed=RNG_SEED)

    df = pd.DataFrame({
        "parameter": problem["names"],
        "S1": Si["S1"],
        "S1_conf": Si["S1_conf"],
        "ST": Si["ST"],
        "ST_conf": Si["ST_conf"],
    })
    df = df.sort_values("ST", ascending=False).reset_index(drop=True)
    df["organism"] = organism
    df["output"] = output

    return df, problem


# ---------------------------------------------------------------------------
# Morris analysis (cross-check)
# ---------------------------------------------------------------------------

def run_morris(
    organism: Organism = "yeast",
    N: int = 512,
    output: str = "margin",
) -> pd.DataFrame:
    """Run Morris elementary effects screening.

    Returns DataFrame with mu_star (importance) and sigma (interaction/nonlinearity).
    """
    include_uG = (output == "margin")
    problem = _build_salib_problem(organism, include_uG=include_uG)

    X = morris_sample.sample(problem, N, seed=RNG_SEED)

    if output == "margin":
        Y = _evaluate_margin(X)
    else:
        Y = _evaluate_critical_glucose(X)

    Si = morris_analyze.analyze(problem, X, Y)

    df = pd.DataFrame({
        "parameter": problem["names"],
        "mu_star": Si["mu_star"],
        "sigma": Si["sigma"],
        "mu_star_conf": Si["mu_star_conf"],
    })
    df = df.sort_values("mu_star", ascending=False).reset_index(drop=True)
    df["organism"] = organism
    df["output"] = output

    return df


# ---------------------------------------------------------------------------
# E3b — Sobol on LP phenotype output (nonlinear, non-tautological)
# ---------------------------------------------------------------------------

def _evaluate_lp_frac_glyc(X: np.ndarray) -> np.ndarray:
    """Evaluate Model A (Shen LP) fermentative fraction at a diagnostic glucose level.

    Unlike the margin formula (linear in u_G), the LP solution is genuinely
    nonlinear: it has corner-point transitions and constraint switching.
    The glucose level is set to the proteome-saturating threshold g* = V_R * Phi
    for each sample — the regime where the phenotype decision actually happens.
    """
    n = X.shape[0]
    Y = np.zeros(n)

    for i in range(n):
        gamma_resp = X[i, 0]
        gamma_glyc = X[i, 1]
        V_resp = X[i, 2]
        V_glyc = X[i, 3]
        Phi = X[i, 4]
        u_G = X[i, 5] if X.shape[1] > 5 else 1.0

        g_star = V_resp * Phi
        g_avail = g_star * 1.5

        effective_V_glyc = V_glyc * u_G

        p = ModelParams(
            gamma_resp=gamma_resp,
            gamma_glyc=gamma_glyc,
            V_resp=V_resp,
            V_glyc=effective_V_glyc,
            Phi=Phi,
            g_avail=g_avail,
        )
        result = solve_shen(p)
        Y[i] = result.frac_glyc

    return Y


def run_sobol_lp(
    organism: Organism = "yeast",
    N: int = 1024,
) -> tuple[pd.DataFrame, dict]:
    """Run Sobol on LP model's fermentative fraction (nonlinear output).

    This addresses the tautology in E3: the margin m is linear in u_G,
    so Sobol trivially ranks u_G first. Here, the output is the LP's
    optimal frac_glyc — a genuinely nonlinear function with corner-point
    transitions — giving Sobol something meaningful to decompose.
    """
    problem = _build_salib_problem(organism, include_uG=True)

    X = sobol_sample.sample(problem, N, calc_second_order=False, seed=RNG_SEED)
    Y = _evaluate_lp_frac_glyc(X)

    Si = sobol_analyze.analyze(problem, Y, calc_second_order=False, seed=RNG_SEED)

    df = pd.DataFrame({
        "parameter": problem["names"],
        "S1": Si["S1"],
        "S1_conf": Si["S1_conf"],
        "ST": Si["ST"],
        "ST_conf": Si["ST_conf"],
    })
    df = df.sort_values("ST", ascending=False).reset_index(drop=True)
    df["organism"] = organism
    df["output"] = "lp_frac_glyc"

    return df, problem


def run_e3(
    organisms: list[Organism] | None = None,
    N_sobol: int = 1024,
    N_morris: int = 512,
) -> dict[str, pd.DataFrame]:
    """Run complete E3 identifiability analysis.

    Returns dict of DataFrames: sobol_margin, sobol_gstar, sobol_lp,
    morris_margin, morris_gstar.
    """
    if organisms is None:
        organisms = list(ALL_ORGANISMS.keys())

    RESULTS_DIR.mkdir(exist_ok=True)

    sobol_margin_dfs = []
    sobol_gstar_dfs = []
    sobol_lp_dfs = []
    morris_margin_dfs = []
    morris_gstar_dfs = []

    for org in organisms:
        print(f"[M4] Running Sobol (margin) for {org}...")
        df_sm, _ = run_sobol(org, N=N_sobol, output="margin")
        sobol_margin_dfs.append(df_sm)

        print(f"[M4] Running Sobol (g*) for {org}...")
        df_sg, _ = run_sobol(org, N=N_sobol, output="g_star")
        sobol_gstar_dfs.append(df_sg)

        print(f"[M4] Running Sobol (LP frac_glyc) for {org}...")
        df_sl, _ = run_sobol_lp(org, N=N_sobol)
        sobol_lp_dfs.append(df_sl)

        print(f"[M4] Running Morris (margin) for {org}...")
        df_mm = run_morris(org, N=N_morris, output="margin")
        morris_margin_dfs.append(df_mm)

        print(f"[M4] Running Morris (g*) for {org}...")
        df_mg = run_morris(org, N=N_morris, output="g_star")
        morris_gstar_dfs.append(df_mg)

    sobol_margin = pd.concat(sobol_margin_dfs, ignore_index=True)
    sobol_gstar = pd.concat(sobol_gstar_dfs, ignore_index=True)
    sobol_lp = pd.concat(sobol_lp_dfs, ignore_index=True)
    morris_margin = pd.concat(morris_margin_dfs, ignore_index=True)
    morris_gstar = pd.concat(morris_gstar_dfs, ignore_index=True)

    sobol_margin.to_csv(RESULTS_DIR / "sobol_margin.csv", index=False)
    sobol_gstar.to_csv(RESULTS_DIR / "sobol_gstar.csv", index=False)
    sobol_lp.to_csv(RESULTS_DIR / "sobol_lp.csv", index=False)
    morris_margin.to_csv(RESULTS_DIR / "morris_margin.csv", index=False)
    morris_gstar.to_csv(RESULTS_DIR / "morris_gstar.csv", index=False)

    print("\n[M4] Top parameters by Sobol total-order (margin output):")
    for org in organisms:
        subset = sobol_margin[sobol_margin["organism"] == org].head(3)
        print(f"  {org}:")
        for _, row in subset.iterrows():
            print(f"    {row['parameter']:>12}: ST = {row['ST']:.4f} ± {row['ST_conf']:.4f}")

    print("\n[M4] Top parameters by Sobol total-order (LP frac_glyc — nonlinear):")
    for org in organisms:
        subset = sobol_lp[sobol_lp["organism"] == org].head(3)
        print(f"  {org}:")
        for _, row in subset.iterrows():
            print(f"    {row['parameter']:>12}: ST = {row['ST']:.4f} ± {row['ST_conf']:.4f}")

    print("\n[M4] Cross-check: top parameters by Morris µ* (margin output):")
    for org in organisms:
        subset = morris_margin[morris_margin["organism"] == org].head(3)
        print(f"  {org}:")
        for _, row in subset.iterrows():
            print(f"    {row['parameter']:>12}: µ* = {row['mu_star']:.6f}")

    return {
        "sobol_margin": sobol_margin,
        "sobol_gstar": sobol_gstar,
        "sobol_lp": sobol_lp,
        "morris_margin": morris_margin,
        "morris_gstar": morris_gstar,
    }


if __name__ == "__main__":
    run_e3()
