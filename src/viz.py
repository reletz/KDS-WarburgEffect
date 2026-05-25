"""M5 — Visualization: publication-quality figures for the audit report.

F1: E1 overlay (two models, shared glucose axis) with disagreement shading
F2: E2 capacity vs. realized (verdict flip curve) with confidence intervals
F3: E3 Sobol/Morris tornado (dispute-settling measurement)
F3c: E3b Sobol LP comparison (nonlinear vs linear)
F5: Wang 2025 comparison (crossing-point mechanism)
F6: Attribution decomposition (u_G vs enzyme definition)
F7: E4 phase diagram (2D verdict map)

All figures regenerated from cached results in <30s via `python -m src.viz`.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures")

# Color palette: teal = respiration, amber = glycolysis (consistent across deck & paper)
C_RESP = "#0d7377"
C_GLYC = "#e6a817"
C_MODEL_A = "#2c3e50"
C_MODEL_B = "#c0392b"
C_WANG = "#7b2d8b"
C_CI = "#aaaaaa"


def _setup_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


# ---------------------------------------------------------------------------
# F1 — E1 Overlay: two models on shared glucose axis
# ---------------------------------------------------------------------------

def plot_f1_overlay(df: pd.DataFrame | None = None) -> plt.Figure:
    """Plot E1: both models' fermentative fraction vs. glucose availability.

    Enhanced with shaded disagreement regions where models predict different phenotypes.
    """
    _setup_style()

    if df is None:
        df = pd.read_csv(RESULTS_DIR / "e1_overlay.csv")

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4.5), squeeze=False)

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        sub = df[df["organism"] == org]

        m_a = sub[sub["model"] == "Model A (Shen LP)"].sort_values("g_avail")
        m_b = sub[sub["model"] == "Model B (Kukurugya)"].sort_values("g_avail")

        ax.plot(m_a["g_avail"], m_a["frac_glyc"], color=C_MODEL_A, ls="-", lw=2, label="Model A (Shen LP)")
        ax.plot(m_b["g_avail"], m_b["frac_glyc"], color=C_MODEL_B, ls="--", lw=2, label="Model B (Kukurugya)")

        # Shade disagreement region (where models differ by >0.1)
        if len(m_a) == len(m_b):
            g_vals = m_a["g_avail"].values
            diff = np.abs(m_a["frac_glyc"].values - m_b["frac_glyc"].values)
            disagree = diff > 0.1
            ax.fill_between(
                g_vals, 0, 1, where=disagree,
                alpha=0.08, color="red", label="Disagreement zone"
            )

        ax.set_xlabel("Glucose availability (g_avail)")
        ax.set_ylabel("Fermentative ATP fraction")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(0.5, color="gray", ls=":", lw=0.8, alpha=0.5)
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle("E1: Model Overlay — Fermentative Fraction vs. Glucose Availability", fontsize=13, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# F2 — E2 Capacity vs. Realized (verdict flip)
# ---------------------------------------------------------------------------

def plot_f2_verdict_flip(df: pd.DataFrame | None = None, ci_df: pd.DataFrame | None = None) -> plt.Figure:
    """Plot E2: realized efficiency margin as function of u_G.

    Enhanced with confidence interval bands on the flip point from E2c bootstrap.
    """
    _setup_style()

    if df is None:
        df = pd.read_csv(RESULTS_DIR / "e2_capacity_vs_realized.csv")
    if ci_df is None:
        ci_path = RESULTS_DIR / "e2c_uncertainty.csv"
        ci_df = pd.read_csv(ci_path) if ci_path.exists() else None

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4.5), squeeze=False)

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        sub = df[df["organism"] == org]

        ax.plot(sub["u_G"], sub["Vgamma_R_realized"], color=C_RESP, lw=2, label="Realized (V·γ)_R")
        ax.plot(sub["u_G"], sub["Vgamma_G_realized"], color=C_GLYC, lw=2, label="Realized (V·γ)_G")
        ax.axhline(sub["Vgamma_R_capacity"].iloc[0], color=C_RESP, ls=":", lw=1, alpha=0.6, label="Capacity (V·γ)_R")
        ax.axhline(sub["Vgamma_G_capacity"].iloc[0], color=C_GLYC, ls=":", lw=1, alpha=0.6, label="Capacity (V·γ)_G")

        # Find and annotate flip point
        margin = sub["Vgamma_R_realized"].values - sub["Vgamma_G_realized"].values
        sign_changes = np.where(np.diff(np.sign(margin)))[0]
        for sc in sign_changes:
            flip_u = sub["u_G"].iloc[sc]
            ax.axvline(flip_u, color="red", ls="--", lw=1.5, alpha=0.8)

            # Add CI band if available
            if ci_df is not None:
                ci_row = ci_df[ci_df["organism"] == org]
                if not ci_row.empty:
                    ci_lo = ci_row["ci_5"].iloc[0]
                    ci_hi = ci_row["ci_95"].iloc[0]
                    y_min, y_max = ax.get_ylim() if ax.get_ylim()[1] != 1.0 else (sub["Vgamma_R_realized"].min(), sub["Vgamma_G_realized"].max())
                    ax.axvspan(ci_lo, ci_hi, alpha=0.12, color=C_CI, label=f"90% CI [{ci_lo:.3f}, {ci_hi:.3f}]")
                    ax.annotate(
                        f"flip = {flip_u:.3f}\n90% CI: [{ci_lo:.3f}, {ci_hi:.3f}]",
                        xy=(flip_u, sub["Vgamma_G_realized"].max() * 0.85),
                        fontsize=7, color="red", ha="center",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                    )
                else:
                    ax.annotate(f"flip @ u_G={flip_u:.3f}", xy=(flip_u, sub["Vgamma_G_realized"].max() * 0.9),
                                fontsize=8, color="red", ha="center")
            else:
                ax.annotate(f"flip @ u_G={flip_u:.3f}", xy=(flip_u, sub["Vgamma_G_realized"].max() * 0.9),
                            fontsize=8, color="red", ha="center")

        # Shade verdict regions
        u_vals = sub["u_G"].values
        ax.fill_between(u_vals, sub["Vgamma_R_realized"], sub["Vgamma_G_realized"],
                        where=sub["Vgamma_R_realized"] > sub["Vgamma_G_realized"],
                        alpha=0.06, color=C_RESP)
        ax.fill_between(u_vals, sub["Vgamma_R_realized"], sub["Vgamma_G_realized"],
                        where=sub["Vgamma_R_realized"] < sub["Vgamma_G_realized"],
                        alpha=0.06, color=C_GLYC)

        ax.set_xlabel("Glycolytic utilization fraction (u_G)")
        ax.set_ylabel("ATP rate per protein (V·γ × u)")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))
        ax.legend(loc="upper right", fontsize=7)
        ax.grid(alpha=0.3)

    fig.suptitle("E2: Capacity vs. Realized — Verdict Flips When Glycolysis Is Underutilized", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# F3 — E3 Sobol tornado (dispute-settling measurement)
# ---------------------------------------------------------------------------

def plot_f3_sobol_tornado(df: pd.DataFrame | None = None) -> plt.Figure:
    """Plot E3: horizontal bar chart of Sobol total-order indices."""
    _setup_style()

    if df is None:
        df = pd.read_csv(RESULTS_DIR / "sobol_margin.csv")

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4), squeeze=False)

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        sub = df[df["organism"] == org].sort_values("ST", ascending=True)

        colors = [C_GLYC if "glyc" in p or p == "u_G" else C_RESP if "resp" in p else "#555555"
                  for p in sub["parameter"]]

        bars = ax.barh(sub["parameter"], sub["ST"], xerr=sub["ST_conf"],
                       color=colors, edgecolor="white", linewidth=0.5)

        ax.set_xlabel("Sobol Total-Order Index (ST)")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))
        ax.set_xlim(0, 1.0)
        ax.grid(axis="x", alpha=0.3)

    fig.suptitle("E3: Parameter Importance — Which Measurement Settles the Dispute?", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


def plot_f3_morris_scatter(df: pd.DataFrame | None = None) -> plt.Figure:
    """Plot Morris µ* vs σ scatter for cross-check."""
    _setup_style()

    if df is None:
        df = pd.read_csv(RESULTS_DIR / "morris_margin.csv")

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4), squeeze=False)

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        sub = df[df["organism"] == org]

        for _, row in sub.iterrows():
            color = C_GLYC if "glyc" in row["parameter"] or row["parameter"] == "u_G" else C_RESP if "resp" in row["parameter"] else "#555555"
            ax.scatter(row["mu_star"], row["sigma"], color=color, s=80, zorder=3)
            ax.annotate(row["parameter"], (row["mu_star"], row["sigma"]),
                        fontsize=8, ha="left", va="bottom")

        ax.set_xlabel("µ* (importance)")
        ax.set_ylabel("σ (nonlinearity/interaction)")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))
        ax.grid(alpha=0.3)

    fig.suptitle("E3 Cross-Check: Morris Elementary Effects", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# F4 — M6 ecCore validation (stretch)
# ---------------------------------------------------------------------------

def plot_f4_eccore(
    df_phi: pd.DataFrame | None = None,
    df_glc: pd.DataFrame | None = None,
) -> plt.Figure:
    """Plot M6: enzyme-constrained FBA showing overflow metabolism onset."""
    _setup_style()

    if df_phi is None:
        df_phi = pd.read_csv(RESULTS_DIR / "m6_phi_sweep.csv")
    if df_glc is None:
        df_glc = pd.read_csv(RESULTS_DIR / "m6_glucose_sweep.csv")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: Φ sweep — acetate overflow vs proteome budget
    opt_phi = df_phi[df_phi["status"] == "optimal"]
    ax1.plot(opt_phi["Phi"], opt_phi["growth_rate"], color=C_MODEL_A, lw=2, label="Growth rate")
    ax1_twin = ax1.twinx()
    ax1_twin.plot(opt_phi["Phi"], opt_phi["acetate_secretion"], color=C_GLYC, lw=2, ls="--", label="Acetate secretion")
    ax1_twin.fill_between(opt_phi["Phi"], 0, opt_phi["acetate_secretion"],
                          where=opt_phi["acetate_secretion"] > 0.1,
                          alpha=0.15, color=C_GLYC)
    ax1.set_xlabel("Proteome budget Φ (kDa·mmol/gDW/h)")
    ax1.set_ylabel("Growth rate (1/h)", color=C_MODEL_A)
    ax1_twin.set_ylabel("Acetate secretion (mmol/gDW/h)", color=C_GLYC)
    ax1.set_title("Φ-Sweep: Overflow Onset Under Proteome Constraint")
    ax1.grid(alpha=0.3)
    lines1 = ax1.get_lines() + ax1_twin.get_lines()
    labels1 = [l.get_label() for l in lines1]
    ax1.legend(lines1, labels1, loc="center right", fontsize=8)

    # Right: Glucose sweep at tight Φ — fermentation fraction vs glucose
    opt_glc = df_glc[df_glc["status"] == "optimal"]
    ax2.plot(opt_glc["glucose_uptake"], opt_glc["frac_fermentative"], color=C_GLYC, lw=2, label="Ferm. fraction")
    ax2.fill_between(opt_glc["glucose_uptake"], 0, opt_glc["frac_fermentative"], alpha=0.15, color=C_GLYC)
    ax2_twin = ax2.twinx()
    ax2_twin.plot(opt_glc["glucose_uptake"], opt_glc["growth_rate"], color=C_MODEL_A, lw=2, ls="--", label="Growth rate")
    ax2.set_xlabel("Glucose uptake (mmol/gDW/h)")
    ax2.set_ylabel("Fermentative carbon fraction", color=C_GLYC)
    ax2_twin.set_ylabel("Growth rate (1/h)", color=C_MODEL_A)
    ax2.set_title(f"Glucose Sweep at Tight Φ = {opt_glc['Phi'].iloc[0]:.0f}")
    ax2.set_ylim(-0.05, 1.05)
    ax2.grid(alpha=0.3)
    lines2 = ax2.get_lines() + ax2_twin.get_lines()
    labels2 = [l.get_label() for l in lines2]
    ax2.legend(lines2, labels2, loc="upper left", fontsize=8)

    fig.suptitle("M6 Validation: ecCore Enzyme-Constrained FBA Confirms Overflow Transition", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# F5 — Wang 2025 comparison (crossing-point mechanism)
# ---------------------------------------------------------------------------

def plot_f5_wang_comparison(df: pd.DataFrame | None = None) -> plt.Figure:
    """Plot comparison between our u_G mechanism and Wang 2025's crossing prediction.

    Shows that both frameworks predict the same efficiency crossing but from
    different mechanistic angles: Wang uses nutrient quality, we use utilization fraction.
    """
    _setup_style()

    if df is None:
        df_path = RESULTS_DIR / "wang_comparison.csv"
        if df_path.exists():
            df = pd.read_csv(df_path)
        else:
            from src.audit import run_wang_comparison
            df = run_wang_comparison()

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(2, len(organisms), figsize=(5 * len(organisms), 8), squeeze=False)

    for idx, org in enumerate(organisms):
        sub = df[df["organism"] == org]
        crossing = sub["crossing_point_uG"].iloc[0]

        # Top row: Our framework (u_G axis)
        ax_top = axes[0, idx]
        ax_top.plot(sub["u_G"], sub["realized_eff_resp"], color=C_RESP, lw=2, label="ε_resp (realized)")
        ax_top.plot(sub["u_G"], sub["realized_eff_glyc"], color=C_GLYC, lw=2, label="ε_glyc (realized)")
        ax_top.axvline(crossing, color="red", ls="--", lw=1.5, alpha=0.8)
        ax_top.fill_between(sub["u_G"], sub["realized_eff_resp"], sub["realized_eff_glyc"],
                            where=sub["realized_eff_resp"] > sub["realized_eff_glyc"],
                            alpha=0.08, color=C_RESP)
        ax_top.fill_between(sub["u_G"], sub["realized_eff_resp"], sub["realized_eff_glyc"],
                            where=sub["realized_eff_resp"] < sub["realized_eff_glyc"],
                            alpha=0.08, color=C_GLYC)
        ax_top.annotate(f"crossing @ ρ = {crossing:.3f}", xy=(crossing, ax_top.get_ylim()[0]),
                        xytext=(crossing + 0.05, sub["realized_eff_resp"].iloc[0] * 0.5),
                        fontsize=8, color="red", arrowprops=dict(arrowstyle="->", color="red"))
        ax_top.set_xlabel("Glycolytic utilization fraction (u_G)")
        ax_top.set_ylabel("Efficiency (V·γ × u)")
        ax_top.set_title(f"{org} — Our Framework")
        ax_top.legend(loc="upper left", fontsize=8)
        ax_top.grid(alpha=0.3)

        # Bottom row: Wang's framework (nutrient quality axis)
        ax_bot = axes[1, idx]
        wang_crossing = sub["crossing_point_wang_q"].iloc[0]
        ax_bot.plot(sub["wang_nutrient_quality"], sub["realized_eff_resp"], color=C_RESP, lw=2, label="ε_resp")
        ax_bot.plot(sub["wang_nutrient_quality"], sub["realized_eff_glyc"], color=C_GLYC, lw=2, label="ε_glyc")
        ax_bot.axvline(wang_crossing, color=C_WANG, ls="--", lw=1.5, alpha=0.8)
        ax_bot.fill_between(sub["wang_nutrient_quality"], sub["realized_eff_resp"], sub["realized_eff_glyc"],
                            where=sub["realized_eff_resp"] > sub["realized_eff_glyc"],
                            alpha=0.08, color=C_RESP)
        ax_bot.fill_between(sub["wang_nutrient_quality"], sub["realized_eff_resp"], sub["realized_eff_glyc"],
                            where=sub["realized_eff_resp"] < sub["realized_eff_glyc"],
                            alpha=0.08, color=C_GLYC)
        ax_bot.annotate(f"Wang crossing @ q = {wang_crossing:.3f}", xy=(wang_crossing, ax_bot.get_ylim()[0]),
                        xytext=(wang_crossing + 0.05, sub["realized_eff_resp"].iloc[0] * 0.5),
                        fontsize=8, color=C_WANG, arrowprops=dict(arrowstyle="->", color=C_WANG))
        ax_bot.set_xlabel("Wang nutrient quality (q = 1 − u_G)")
        ax_bot.set_ylabel("Efficiency (V·γ × u)")
        ax_bot.set_title(f"{org} — Wang 2025 Framework")
        ax_bot.legend(loc="upper right", fontsize=8)
        ax_bot.grid(alpha=0.3)

    fig.suptitle("Comparison: Our Audit vs. Wang 2025 — Same Crossing, Different Mechanism", fontsize=12, y=1.01)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# F6 — Attribution decomposition
# ---------------------------------------------------------------------------

def plot_f6_decomposition(df: pd.DataFrame | None = None) -> plt.Figure:
    """Plot E2d: how enzyme reattribution (varying which proteins are counted) affects ρ and flip.

    Key insight: even ±30% reattribution of enzyme mass between pathways shifts
    the flip point proportionally (flip ≈ ρ always), but the Sobol analysis shows
    u_G dominates (ST ≈ 0.76-0.78) while V_glyc/V_resp contribute only ST ≈ 0.12.
    """
    _setup_style()

    if df is None:
        df_path = RESULTS_DIR / "e2d_decomposition.csv"
        if df_path.exists():
            df = pd.read_csv(df_path)
        else:
            return plt.figure()

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4.5), squeeze=False)

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        sub = df[df["organism"] == org].dropna(subset=["flip_uG"])

        if len(sub) == 0:
            continue

        # Plot flip_uG vs reattribution factor
        ax.plot(sub["reattribution_factor"], sub["flip_uG"], color=C_MODEL_A, lw=2.5, label="flip u_G")
        ax.plot(sub["reattribution_factor"], sub["rho"], color="red", ls="--", lw=1.5, alpha=0.7, label="ρ (tracks flip)")

        # Mark baseline
        baseline_rho = sub["baseline_rho"].iloc[0]
        ax.axhline(baseline_rho, color="gray", ls=":", lw=1, alpha=0.5)
        ax.axvline(1.0, color="gray", ls=":", lw=1, alpha=0.5)
        ax.annotate(f"baseline ρ = {baseline_rho:.3f}",
                    xy=(1.0, baseline_rho), xytext=(1.1, baseline_rho * 1.4),
                    fontsize=8, color="gray",
                    arrowprops=dict(arrowstyle="->", color="gray"))

        # Show the range of effect
        flip_range = sub["flip_uG"].max() - sub["flip_uG"].min()
        ax.annotate(
            f"±30% reattribution\nshifts flip by {flip_range:.3f}\n"
            f"(but u_G effect spans\nfull range [0,1])",
            xy=(0.85, sub["flip_uG"].min()),
            fontsize=7, ha="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.9),
        )

        ax.set_xlabel("Enzyme reattribution factor\n(<1 = fewer glyc. enzymes counted; >1 = more)")
        ax.set_ylabel("Verdict flip point (u_G = ρ)")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle("E2d: Enzyme Definition Shifts ρ, But u_G Still Dominates Verdict (ST > 0.75)", fontsize=11, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Generate all figures
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# F3c — E3b Sobol LP comparison (nonlinear vs linear)
# ---------------------------------------------------------------------------

def plot_f3c_sobol_lp(
    df_lp: pd.DataFrame | None = None,
    df_margin: pd.DataFrame | None = None,
) -> plt.Figure:
    """Compare Sobol indices from linear margin vs nonlinear LP output.

    Shows side-by-side bars: left = margin (tautological), right = LP frac_glyc
    (genuinely nonlinear). If u_G still dominates the LP output, the finding
    is robust and not an artifact of the linear formula.
    """
    _setup_style()

    if df_lp is None:
        df_lp = pd.read_csv(RESULTS_DIR / "sobol_lp.csv")
    if df_margin is None:
        df_margin = pd.read_csv(RESULTS_DIR / "sobol_margin.csv")

    organisms = df_lp["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4.5), squeeze=False)

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        lp = df_lp[df_lp["organism"] == org].sort_values("parameter")
        mg = df_margin[df_margin["organism"] == org].sort_values("parameter")

        params = lp["parameter"].values
        x = np.arange(len(params))
        w = 0.35

        colors_lp = [C_GLYC if "glyc" in p or p == "u_G" else C_RESP if "resp" in p else "#555555"
                     for p in params]

        ax.barh(x - w/2, lp["ST"].values, w, xerr=lp["ST_conf"].values,
                color=colors_lp, alpha=0.9, label="LP frac_glyc (nonlinear)")
        ax.barh(x + w/2, mg["ST"].values, w, xerr=mg["ST_conf"].values,
                color=[c + "66" if len(c) == 7 else c for c in colors_lp],
                alpha=0.5, edgecolor="gray", linewidth=0.5,
                label="Margin (linear)")

        ax.set_yticks(x)
        ax.set_yticklabels(params)
        ax.set_xlabel("Sobol Total-Order Index (ST)")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))
        ax.set_xlim(0, 1.0)
        ax.legend(loc="lower right", fontsize=7)
        ax.grid(axis="x", alpha=0.3)

    fig.suptitle("E3b: Sobol on LP Output (Nonlinear) vs. Margin (Linear)", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# F7 — E4 Phase diagram (2D verdict map)
# ---------------------------------------------------------------------------

def plot_f7_phase_diagram(df: pd.DataFrame | None = None) -> plt.Figure:
    """Plot E4: 2D heatmap of fermentative fraction over (u_G, g_avail).

    Shows regime boundaries between respiration-dominant, mixed, and
    glycolysis-dominant regions. This phase diagram is a genuinely new
    visualization not present in Shen, Kukurugya, or Wang.
    """
    _setup_style()

    if df is None:
        df = pd.read_csv(RESULTS_DIR / "e4_phase_diagram.csv")

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4.5), squeeze=False)

    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "verdict", [C_RESP, "#f0f0f0", C_GLYC], N=256
    )

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        sub = df[df["organism"] == org]

        u_vals = sorted(sub["u_G"].unique())
        g_vals = sorted(sub["g_avail"].unique())

        Z = np.zeros((len(u_vals), len(g_vals)))
        for i, u in enumerate(u_vals):
            for j, g in enumerate(g_vals):
                row = sub[(sub["u_G"] == u) & (sub["g_avail"] == g)]
                if len(row) > 0:
                    Z[i, j] = row["frac_glyc"].iloc[0]

        im = ax.pcolormesh(
            g_vals, u_vals, Z,
            cmap=cmap, vmin=0, vmax=1, shading="auto",
        )

        # Draw regime boundary contours
        ax.contour(
            g_vals, u_vals, Z,
            levels=[0.05, 0.5, 0.95],
            colors=["white", "black", "white"],
            linewidths=[1, 2, 1],
            linestyles=["--", "-", "--"],
        )

        ax.set_xlabel("Glucose availability (g_avail)")
        ax.set_ylabel("Glycolytic utilization (u_G)")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))

    cbar = fig.colorbar(im, ax=axes[0, :], shrink=0.8, pad=0.02)
    cbar.set_label("Fermentative ATP fraction")

    fig.suptitle("E4: Phase Diagram — Metabolic Regime Map (u_G vs. Glucose)", fontsize=12, y=1.02)
    return fig


def generate_all_figures(include_m6: bool = False) -> None:
    """Generate and save all report/deck figures from cached results."""
    FIGURES_DIR.mkdir(exist_ok=True)

    print("[M5] Generating F1 (E1 overlay with disagreement zones)...")
    fig1 = plot_f1_overlay()
    fig1.savefig(FIGURES_DIR / "F1_overlay.png")
    fig1.savefig(FIGURES_DIR / "F1_overlay.svg")
    plt.close(fig1)

    print("[M5] Generating F2 (E2 verdict flip with CI)...")
    fig2 = plot_f2_verdict_flip()
    fig2.savefig(FIGURES_DIR / "F2_verdict_flip.png")
    fig2.savefig(FIGURES_DIR / "F2_verdict_flip.svg")
    plt.close(fig2)

    print("[M5] Generating F3 (E3 Sobol tornado)...")
    fig3 = plot_f3_sobol_tornado()
    fig3.savefig(FIGURES_DIR / "F3_sobol_tornado.png")
    fig3.savefig(FIGURES_DIR / "F3_sobol_tornado.svg")
    plt.close(fig3)

    print("[M5] Generating F3b (Morris cross-check)...")
    fig3b = plot_f3_morris_scatter()
    fig3b.savefig(FIGURES_DIR / "F3b_morris_scatter.png")
    fig3b.savefig(FIGURES_DIR / "F3b_morris_scatter.svg")
    plt.close(fig3b)

    if include_m6:
        print("[M5] Generating F4 (M6 ecCore validation)...")
        fig4 = plot_f4_eccore()
        fig4.savefig(FIGURES_DIR / "F4_eccore.png")
        fig4.savefig(FIGURES_DIR / "F4_eccore.svg")
        plt.close(fig4)

    print("[M5] Generating F5 (Wang 2025 comparison)...")
    fig5 = plot_f5_wang_comparison()
    fig5.savefig(FIGURES_DIR / "F5_wang_comparison.png")
    fig5.savefig(FIGURES_DIR / "F5_wang_comparison.svg")
    plt.close(fig5)

    print("[M5] Generating F6 (attribution decomposition)...")
    fig6 = plot_f6_decomposition()
    fig6.savefig(FIGURES_DIR / "F6_decomposition.png")
    fig6.savefig(FIGURES_DIR / "F6_decomposition.svg")
    plt.close(fig6)

    if (RESULTS_DIR / "sobol_lp.csv").exists():
        print("[M5] Generating F3c (Sobol LP vs margin comparison)...")
        fig3c = plot_f3c_sobol_lp()
        fig3c.savefig(FIGURES_DIR / "F3c_sobol_lp.png")
        fig3c.savefig(FIGURES_DIR / "F3c_sobol_lp.svg")
        plt.close(fig3c)

    if (RESULTS_DIR / "e4_phase_diagram.csv").exists():
        print("[M5] Generating F7 (phase diagram)...")
        fig7 = plot_f7_phase_diagram()
        fig7.savefig(FIGURES_DIR / "F7_phase_diagram.png")
        fig7.savefig(FIGURES_DIR / "F7_phase_diagram.svg")
        plt.close(fig7)

    print(f"[M5] All figures saved to {FIGURES_DIR}/")


if __name__ == "__main__":
    generate_all_figures()
