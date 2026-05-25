"""M5 — Visualization: publication-quality figures for the audit report.

F1: E1 overlay (two models, shared glucose axis)
F2: E2 capacity vs. realized (verdict flip curve)
F3: E3 Sobol/Morris tornado (dispute-settling measurement)

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
    """Plot E1: both models' fermentative fraction vs. glucose availability."""
    _setup_style()

    if df is None:
        df = pd.read_csv(RESULTS_DIR / "e1_overlay.csv")

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4), squeeze=False)

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        sub = df[df["organism"] == org]

        for model_name, color, ls in [
            ("Model A (Shen LP)", C_MODEL_A, "-"),
            ("Model B (Kukurugya)", C_MODEL_B, "--"),
        ]:
            m = sub[sub["model"] == model_name]
            ax.plot(m["g_avail"], m["frac_glyc"], color=color, ls=ls, lw=2, label=model_name)

        ax.set_xlabel("Glucose availability (g_avail)")
        ax.set_ylabel("Fermentative ATP fraction")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(0.5, color="gray", ls=":", lw=0.8, alpha=0.5)
        ax.legend(loc="lower right")
        ax.grid(alpha=0.3)

    fig.suptitle("E1: Model Overlay — Fermentative Fraction vs. Glucose Availability", fontsize=13, y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# F2 — E2 Capacity vs. Realized (verdict flip)
# ---------------------------------------------------------------------------

def plot_f2_verdict_flip(df: pd.DataFrame | None = None) -> plt.Figure:
    """Plot E2: realized efficiency margin as function of u_G."""
    _setup_style()

    if df is None:
        df = pd.read_csv(RESULTS_DIR / "e2_capacity_vs_realized.csv")

    organisms = df["organism"].unique()
    fig, axes = plt.subplots(1, len(organisms), figsize=(5 * len(organisms), 4), squeeze=False)

    for idx, org in enumerate(organisms):
        ax = axes[0, idx]
        sub = df[df["organism"] == org]

        ax.plot(sub["u_G"], sub["Vgamma_R_realized"], color=C_RESP, lw=2, label="Realized (V·γ)_R")
        ax.plot(sub["u_G"], sub["Vgamma_G_realized"], color=C_GLYC, lw=2, label="Realized (V·γ)_G")
        ax.axhline(sub["Vgamma_R_capacity"].iloc[0], color=C_RESP, ls=":", lw=1, alpha=0.6, label="Capacity (V·γ)_R")
        ax.axhline(sub["Vgamma_G_capacity"].iloc[0], color=C_GLYC, ls=":", lw=1, alpha=0.6, label="Capacity (V·γ)_G")

        margin = sub["Vgamma_R_realized"].values - sub["Vgamma_G_realized"].values
        sign_changes = np.where(np.diff(np.sign(margin)))[0]
        for sc in sign_changes:
            flip_u = sub["u_G"].iloc[sc]
            ax.axvline(flip_u, color="red", ls="--", lw=1.5, alpha=0.8)
            ax.annotate(f"flip @ u_G={flip_u:.3f}", xy=(flip_u, ax.get_ylim()[1] * 0.9),
                        fontsize=8, color="red", ha="center")

        ax.set_xlabel("Glycolytic utilization fraction (u_G)")
        ax.set_ylabel("ATP rate per protein (V·γ × u)")
        ax.set_title(org.replace("ecoli", "E. coli").replace("yeast", "S. cerevisiae").replace("mammalian", "Mammalian"))
        ax.legend(loc="upper right", fontsize=8)
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
# Generate all figures
# ---------------------------------------------------------------------------

def generate_all_figures() -> None:
    """Generate and save all report/deck figures from cached results."""
    FIGURES_DIR.mkdir(exist_ok=True)

    print("[M5] Generating F1 (E1 overlay)...")
    fig1 = plot_f1_overlay()
    fig1.savefig(FIGURES_DIR / "F1_overlay.png")
    fig1.savefig(FIGURES_DIR / "F1_overlay.svg")
    plt.close(fig1)

    print("[M5] Generating F2 (E2 verdict flip)...")
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

    print(f"[M5] All figures saved to {FIGURES_DIR}/")


if __name__ == "__main__":
    generate_all_figures()
