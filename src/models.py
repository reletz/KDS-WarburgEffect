"""M1 — Two solvers (Model A: Shen-style LP, Model B: Kukurugya analytical), one interface.

Both accept ModelParams and return ModelResult. This is the only place model
mathematics lives; called everywhere else.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.optimize import linprog

from src.params import ModelParams


Phenotype = Literal["respiration", "glycolysis", "mixed"]


@dataclass(frozen=True)
class ModelResult:
    """Unified output from either solver."""

    phenotype: Phenotype
    phi_G: float          # proteome fraction allocated to glycolysis
    phi_R: float          # proteome fraction allocated to respiration
    J_ATP: float          # total ATP production rate
    J_ATP_glyc: float     # ATP rate from glycolysis
    J_ATP_resp: float     # ATP rate from respiration
    frac_glyc: float      # fraction of ATP from glycolysis (0 = pure resp, 1 = pure glyc)
    eps_ratio: float      # realized ρ = J_ATP_resp/phi_R ÷ J_ATP_glyc/phi_G (if both >0)
    g_avail: float        # input glucose availability (for tracing)


# ---------------------------------------------------------------------------
# Model A — 2-sector proteome LP (Shen/Basan style)
# ---------------------------------------------------------------------------
# Decision variables: [J_glc_G, J_glc_R, phi_G, phi_R]
#
# Maximize: gamma_G * J_glc_G + gamma_R * J_glc_R
#
# Subject to:
#   J_glc_G - V_G * phi_G <= 0          (glucose uptake limited by enzyme capacity)
#   J_glc_R - V_R * phi_R <= 0
#   phi_G + phi_R <= Phi                 (proteome budget)
#   J_glc_G + J_glc_R <= g_avail        (glucose availability)
#   all variables >= 0
# ---------------------------------------------------------------------------

def solve_shen(params: ModelParams) -> ModelResult:
    """Solve Model A: 2-sector LP maximizing ATP rate.

    Uses scipy.optimize.linprog with HiGHS backend.
    Variables: x = [J_glc_G, J_glc_R, phi_G, phi_R]
    """
    gamma_G = params.gamma_glyc
    gamma_R = params.gamma_resp
    V_G = params.V_glyc
    V_R = params.V_resp
    Phi = params.Phi
    g = params.g_avail

    # linprog minimizes c^T x, so negate for maximization
    c = [-gamma_G, -gamma_R, 0.0, 0.0]

    # Inequality constraints: A_ub @ x <= b_ub
    # (1) J_glc_G - V_G * phi_G <= 0
    # (2) J_glc_R - V_R * phi_R <= 0
    # (3) phi_G + phi_R <= Phi
    # (4) J_glc_G + J_glc_R <= g_avail
    A_ub = [
        [1.0, 0.0, -V_G, 0.0],    # J_glc_G <= V_G * phi_G
        [0.0, 1.0, 0.0, -V_R],    # J_glc_R <= V_R * phi_R
        [0.0, 0.0, 1.0, 1.0],     # phi_G + phi_R <= Phi
        [1.0, 1.0, 0.0, 0.0],     # J_glc_G + J_glc_R <= g_avail
    ]
    b_ub = [0.0, 0.0, Phi, g]

    # Bounds: all >= 0
    bounds = [(0, None), (0, None), (0, None), (0, None)]

    result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")

    if not result.success:
        raise RuntimeError(f"LP infeasible or unbounded: {result.message}")

    J_glc_G, J_glc_R, phi_G, phi_R = result.x
    J_ATP_G = gamma_G * J_glc_G
    J_ATP_R = gamma_R * J_glc_R
    J_ATP = J_ATP_G + J_ATP_R

    frac_glyc = J_ATP_G / J_ATP if J_ATP > 1e-12 else 0.0

    # Classify phenotype
    if frac_glyc < 0.01:
        phenotype: Phenotype = "respiration"
    elif frac_glyc > 0.99:
        phenotype = "glycolysis"
    else:
        phenotype = "mixed"

    # Realized efficiency ratio (ρ): (J_ATP_R / phi_R) / (J_ATP_G / phi_G)
    if phi_G > 1e-12 and phi_R > 1e-12:
        eps_ratio = (J_ATP_R / phi_R) / (J_ATP_G / phi_G)
    else:
        eps_ratio = float("inf") if phi_G < 1e-12 and phi_R > 1e-12 else 0.0

    return ModelResult(
        phenotype=phenotype,
        phi_G=phi_G,
        phi_R=phi_R,
        J_ATP=J_ATP,
        J_ATP_glyc=J_ATP_G,
        J_ATP_resp=J_ATP_R,
        frac_glyc=frac_glyc,
        eps_ratio=eps_ratio,
        g_avail=g,
    )


# ---------------------------------------------------------------------------
# Model B — Kukurugya 5-parameter analytical model
# ---------------------------------------------------------------------------
# The model maximizes ATP rate J_ATP = gamma_G * V_G * phi_G + gamma_R * V_R * phi_R
# subject to:
#   phi_G + phi_R <= Phi
#   V_G * phi_G + V_R * phi_R <= g_avail   (total glucose uptake <= availability)
#   phi_G, phi_R >= 0
#
# This is still an LP but small enough for analytical corner-point solution:
#
# Three regimes depending on which constraint binds:
#   (i)  Glucose-limited (g_avail small): allocate to maximize yield per glucose
#        → pure respiration if gamma_R > gamma_G (always true biologically)
#   (ii) Proteome-limited (g_avail large): allocate to maximize ATP per protein
#        → pure pathway with highest V*gamma
#   (iii) Transition: both constraints bind simultaneously → mixed
#
# Critical glucose threshold g* where transition occurs:
#   At pure respiration: J_glc = V_R * Phi, need g >= V_R * Phi
#   At transition: both constraints bind → solve 2×2 system
# ---------------------------------------------------------------------------

def solve_kukurugya(params: ModelParams) -> ModelResult:
    """Solve Model B: Kukurugya 5-parameter analytical model.

    Returns optimal allocation at given glucose availability.
    """
    gamma_G = params.gamma_glyc
    gamma_R = params.gamma_resp
    V_G = params.V_glyc
    V_R = params.V_resp
    Phi = params.Phi
    g = params.g_avail

    # Corner solutions:
    # Pure respiration: phi_R = min(Phi, g/V_R), phi_G = 0
    # Pure glycolysis:  phi_G = min(Phi, g/V_G), phi_R = 0
    # Mixed: both constraints bind → phi_G + phi_R = Phi AND V_G*phi_G + V_R*phi_R = g

    # ATP rate at pure respiration
    phi_R_pure = min(Phi, g / V_R) if V_R > 0 else Phi
    J_ATP_resp_pure = gamma_R * V_R * phi_R_pure

    # ATP rate at pure glycolysis
    phi_G_pure = min(Phi, g / V_G) if V_G > 0 else Phi
    J_ATP_glyc_pure = gamma_G * V_G * phi_G_pure

    # Check mixed solution: both constraints bind
    # phi_G + phi_R = Phi
    # V_G * phi_G + V_R * phi_R = g
    # → phi_G = (g - V_R * Phi) / (V_G - V_R)
    # → phi_R = Phi - phi_G
    mixed_feasible = False
    phi_G_mix, phi_R_mix = 0.0, 0.0
    J_ATP_mix = -float("inf")

    if abs(V_G - V_R) > 1e-15:
        phi_G_mix = (g - V_R * Phi) / (V_G - V_R)
        phi_R_mix = Phi - phi_G_mix
        if phi_G_mix >= -1e-12 and phi_R_mix >= -1e-12:
            phi_G_mix = max(0.0, phi_G_mix)
            phi_R_mix = max(0.0, phi_R_mix)
            # Check glucose constraint is satisfied
            g_used = V_G * phi_G_mix + V_R * phi_R_mix
            if g_used <= g + 1e-12:
                mixed_feasible = True
                J_ATP_mix = gamma_G * V_G * phi_G_mix + gamma_R * V_R * phi_R_mix

    # Pick the best feasible solution
    candidates = [
        (J_ATP_resp_pure, 0.0, phi_R_pure, "resp_pure"),
        (J_ATP_glyc_pure, phi_G_pure, 0.0, "glyc_pure"),
    ]
    if mixed_feasible:
        candidates.append((J_ATP_mix, phi_G_mix, phi_R_mix, "mixed"))

    best_J, best_phi_G, best_phi_R, best_label = max(candidates, key=lambda t: t[0])

    # Compute outputs
    J_ATP_G = gamma_G * V_G * best_phi_G
    J_ATP_R = gamma_R * V_R * best_phi_R
    J_ATP = J_ATP_G + J_ATP_R
    frac_glyc = J_ATP_G / J_ATP if J_ATP > 1e-12 else 0.0

    if best_label == "resp_pure" or frac_glyc < 0.01:
        phenotype: Phenotype = "respiration"
    elif best_label == "glyc_pure" or frac_glyc > 0.99:
        phenotype = "glycolysis"
    else:
        phenotype = "mixed"

    if best_phi_G > 1e-12 and best_phi_R > 1e-12:
        eps_ratio = (J_ATP_R / best_phi_R) / (J_ATP_G / best_phi_G)
    else:
        eps_ratio = float("inf") if best_phi_G < 1e-12 and best_phi_R > 1e-12 else 0.0

    return ModelResult(
        phenotype=phenotype,
        phi_G=best_phi_G,
        phi_R=best_phi_R,
        J_ATP=J_ATP,
        J_ATP_glyc=J_ATP_G,
        J_ATP_resp=J_ATP_R,
        frac_glyc=frac_glyc,
        eps_ratio=eps_ratio,
        g_avail=g,
    )


# ---------------------------------------------------------------------------
# Sweep utility — run a model across a glucose availability grid
# ---------------------------------------------------------------------------

def sweep_glucose(
    params: ModelParams,
    g_range: np.ndarray | None = None,
    model: Literal["shen", "kukurugya"] = "shen",
) -> list[ModelResult]:
    """Run one model across a glucose availability sweep.

    Args:
        params: Base parameters (g_avail will be overwritten per point).
        g_range: Array of glucose availability values. If None, auto-generates
                 a range from 0 to 2× the proteome-saturating level.
        model: Which solver to use.

    Returns:
        List of ModelResult, one per g_avail value.
    """
    solver = solve_shen if model == "shen" else solve_kukurugya

    if g_range is None:
        # Auto range: 0 to well past the proteome-saturating glucose level
        g_max = 2.0 * params.V_glyc * params.Phi
        g_range = np.linspace(1e-6, g_max, 200)

    results = []
    for g in g_range:
        p = ModelParams(
            gamma_resp=params.gamma_resp,
            gamma_glyc=params.gamma_glyc,
            V_resp=params.V_resp,
            V_glyc=params.V_glyc,
            Phi=params.Phi,
            g_avail=float(g),
        )
        results.append(solver(p))

    return results


# ---------------------------------------------------------------------------
# E2 utility — capacity vs. realized efficiency
# ---------------------------------------------------------------------------

def compute_efficiency_with_utilization(
    params: ModelParams,
    u_G: float = 1.0,
    u_R: float = 1.0,
) -> dict[str, float]:
    """Compute capacity and realized ATP efficiency per protein for both sectors.

    Args:
        u_G: utilization fraction for glycolysis (1 = at Vmax, <1 = idle capacity)
        u_R: utilization fraction for respiration

    Returns:
        Dict with capacity and realized efficiencies, plus verdict.
    """
    Vgamma_G_cap = params.V_glyc * params.gamma_glyc   # capacity
    Vgamma_R_cap = params.V_resp * params.gamma_resp   # capacity
    Vgamma_G_real = u_G * Vgamma_G_cap                 # realized
    Vgamma_R_real = u_R * Vgamma_R_cap                 # realized

    margin = Vgamma_R_real - Vgamma_G_real  # >0 means respiration wins

    return {
        "Vgamma_G_capacity": Vgamma_G_cap,
        "Vgamma_R_capacity": Vgamma_R_cap,
        "Vgamma_G_realized": Vgamma_G_real,
        "Vgamma_R_realized": Vgamma_R_real,
        "rho_capacity": Vgamma_R_cap / Vgamma_G_cap if Vgamma_G_cap > 0 else float("inf"),
        "rho_realized": Vgamma_R_real / Vgamma_G_real if Vgamma_G_real > 0 else float("inf"),
        "margin": margin,
        "verdict_capacity": "respiration" if Vgamma_R_cap > Vgamma_G_cap else "glycolysis",
        "verdict_realized": "respiration" if margin > 0 else "glycolysis",
        "u_G": u_G,
        "u_R": u_R,
    }


def find_verdict_flip_uG(
    params: ModelParams,
    u_R: float = 1.0,
    n_points: int = 1000,
) -> float | None:
    """Find the glycolytic utilization fraction u_G where verdict flips.

    Sweeps u_G from 1.0 down to near 0 and finds where
    realized verdict changes from 'respiration wins' to 'glycolysis wins'
    (or vice versa).

    Returns:
        u_G value at flip point, or None if no flip occurs.
    """
    u_G_values = np.linspace(1.0, 0.01, n_points)
    prev_margin = None

    for u_G in u_G_values:
        result = compute_efficiency_with_utilization(params, u_G=float(u_G), u_R=u_R)
        margin = result["margin"]

        if prev_margin is not None and prev_margin * margin < 0:
            # Sign change — interpolate
            return float(u_G)
        prev_margin = margin

    return None
