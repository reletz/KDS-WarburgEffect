"""System tests — end-to-end pipeline integration.

These tests verify that the full pipeline runs correctly as a system:
- M2→M1→M3→M4→M5 produce expected artifacts
- Outputs are internally consistent across modules
- Results match known biological constraints
- Pipeline is reproducible (deterministic with fixed seed)
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).parent.parent


# ===========================================================================
# SYSTEM TEST — Full pipeline execution
# ===========================================================================

class TestPipelineExecution:
    """Test that `python main.py` runs end-to-end and produces all artifacts."""

    def test_main_runs_without_error(self):
        """Full pipeline (no stretch) completes with exit code 0."""
        result = subprocess.run(
            [sys.executable, "main.py"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, f"Pipeline failed:\n{result.stderr}"

    def test_main_produces_results(self):
        """Pipeline generates all expected result CSV files."""
        subprocess.run(
            [sys.executable, "main.py"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            timeout=60,
        )
        results_dir = PROJECT_ROOT / "results"
        expected_files = [
            "provenance.csv",
            "e1_overlay.csv",
            "e2_capacity_vs_realized.csv",
            "e2b_accounting.csv",
            "e2c_uncertainty.csv",
            "e2d_decomposition.csv",
            "wang_comparison.csv",
            "sobol_margin.csv",
            "sobol_gstar.csv",
            "morris_margin.csv",
            "morris_gstar.csv",
        ]
        for f in expected_files:
            assert (results_dir / f).exists(), f"Missing result: {f}"

    def test_main_produces_figures(self):
        """Pipeline generates all expected figure files."""
        subprocess.run(
            [sys.executable, "main.py"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            timeout=60,
        )
        figures_dir = PROJECT_ROOT / "figures"
        expected_figures = [
            "F1_overlay.png",
            "F2_verdict_flip.png",
            "F3_sobol_tornado.png",
            "F3b_morris_scatter.png",
            "F5_wang_comparison.png",
            "F6_decomposition.png",
        ]
        for f in expected_figures:
            assert (figures_dir / f).exists(), f"Missing figure: {f}"
            assert (figures_dir / f).stat().st_size > 1000, f"Figure too small: {f}"

    def test_pipeline_stretch_runs(self):
        """Pipeline with --stretch completes (requires cobra)."""
        result = subprocess.run(
            [sys.executable, "main.py", "--stretch"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode == 0, f"Stretch pipeline failed:\n{result.stderr}"
        assert (PROJECT_ROOT / "results" / "m6_phi_sweep.csv").exists()
        assert (PROJECT_ROOT / "figures" / "F4_eccore.png").exists()


# ===========================================================================
# SYSTEM TEST — Cross-module data consistency
# ===========================================================================

class TestCrossModuleConsistency:
    """Verify that outputs from different modules are consistent with each other."""

    @pytest.fixture(autouse=True)
    def run_pipeline(self):
        """Ensure pipeline has been run before these tests."""
        subprocess.run(
            [sys.executable, "main.py"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            timeout=60,
        )

    def test_e1_organisms_match_provenance(self):
        """E1 overlay should contain the same organisms as provenance table."""
        prov = pd.read_csv(PROJECT_ROOT / "results" / "provenance.csv")
        e1 = pd.read_csv(PROJECT_ROOT / "results" / "e1_overlay.csv")
        prov_organisms = set(prov["organism"].unique())
        e1_organisms = set(e1["organism"].unique())
        assert prov_organisms == e1_organisms

    def test_e2_flip_consistent_with_derived_rho(self):
        """E2 flip points should approximately equal the capacity ρ values."""
        from src.params import ALL_ORGANISMS, ModelParams
        from src.models import find_verdict_flip_uG

        for org_name, op in ALL_ORGANISMS.items():
            p = ModelParams.from_organism_params(op, g_avail=0.1)
            flip = find_verdict_flip_uG(p, u_R=1.0)
            rho = op.rho
            assert flip is not None
            assert abs(flip - rho) < 0.02, (
                f"{org_name}: flip u_G ({flip:.4f}) should ≈ ρ ({rho:.4f})"
            )

    def test_e2c_uncertainty_valid(self):
        """E2c uncertainty CSV should have valid CI bounds for all organisms."""
        ci = pd.read_csv(PROJECT_ROOT / "results" / "e2c_uncertainty.csv")
        assert len(ci) == 3
        for _, row in ci.iterrows():
            assert row["ci_5"] < row["ci_95"]
            assert row["n_valid"] > 100  # most bootstraps should succeed

    def test_e2d_decomposition_consistent(self):
        """E2d decomposition should show flip ≈ ρ across all reattribution levels."""
        df = pd.read_csv(PROJECT_ROOT / "results" / "e2d_decomposition.csv")
        valid = df.dropna(subset=["flip_uG"])
        assert len(valid) > 0
        for _, row in valid.iterrows():
            assert abs(row["flip_uG"] - row["rho"]) < 0.02

    def test_sobol_morris_top_param_agreement(self):
        """Sobol and Morris should agree on the top parameter."""
        sobol = pd.read_csv(PROJECT_ROOT / "results" / "sobol_margin.csv")
        morris = pd.read_csv(PROJECT_ROOT / "results" / "morris_margin.csv")

        for org in sobol["organism"].unique():
            sobol_top = sobol[sobol["organism"] == org].sort_values("ST", ascending=False).iloc[0]["parameter"]
            morris_top = morris[morris["organism"] == org].sort_values("mu_star", ascending=False).iloc[0]["parameter"]
            assert sobol_top == morris_top, (
                f"{org}: Sobol top={sobol_top}, Morris top={morris_top} — should agree"
            )

    def test_e1_both_models_same_row_count(self):
        """Each model in E1 should have equal number of data points per organism."""
        e1 = pd.read_csv(PROJECT_ROOT / "results" / "e1_overlay.csv")
        for org in e1["organism"].unique():
            sub = e1[e1["organism"] == org]
            counts = sub.groupby("model").size()
            assert counts.nunique() == 1, (
                f"{org}: models have different row counts: {counts.to_dict()}"
            )


# ===========================================================================
# SYSTEM TEST — Biological invariants (domain constraints)
# ===========================================================================

class TestBiologicalInvariants:
    """Verify that results respect known biological constraints."""

    @pytest.fixture(autouse=True)
    def run_pipeline(self):
        subprocess.run(
            [sys.executable, "main.py"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            timeout=60,
        )

    def test_low_glucose_always_respiration(self):
        """At lowest glucose in E1, both models should show near-pure respiration."""
        e1 = pd.read_csv(PROJECT_ROOT / "results" / "e1_overlay.csv")
        for org in e1["organism"].unique():
            sub = e1[e1["organism"] == org]
            for model in sub["model"].unique():
                m = sub[sub["model"] == model].sort_values("g_avail")
                lowest = m.iloc[0]
                assert lowest["frac_glyc"] < 0.05, (
                    f"{org}/{model}: frac_glyc={lowest['frac_glyc']:.4f} at lowest glucose — "
                    "should be near 0 (respiration)"
                )

    def test_frac_glyc_bounded_in_e1(self):
        """Fermentative fraction must always be in [0, 1]."""
        e1 = pd.read_csv(PROJECT_ROOT / "results" / "e1_overlay.csv")
        assert (e1["frac_glyc"] >= -1e-6).all()
        assert (e1["frac_glyc"] <= 1.0 + 1e-6).all()

    def test_atp_nonnegative_in_e1(self):
        """ATP production must be non-negative everywhere."""
        e1 = pd.read_csv(PROJECT_ROOT / "results" / "e1_overlay.csv")
        assert (e1["J_ATP"] >= -1e-9).all()
        assert (e1["J_ATP_glyc"] >= -1e-9).all()
        assert (e1["J_ATP_resp"] >= -1e-9).all()

    def test_proteome_fractions_valid(self):
        """phi_G and phi_R must be non-negative and sum ��� Phi."""
        e1 = pd.read_csv(PROJECT_ROOT / "results" / "e1_overlay.csv")
        assert (e1["phi_G"] >= -1e-9).all()
        assert (e1["phi_R"] >= -1e-9).all()
        # Phi is at most 0.26 (ecoli), so phi_G + phi_R ≤ 0.30 with tolerance
        assert ((e1["phi_G"] + e1["phi_R"]) <= 0.30).all()

    def test_e2_margin_monotone_with_uG(self):
        """In E2, margin should be monotonically decreasing as u_G increases."""
        e2 = pd.read_csv(PROJECT_ROOT / "results" / "e2_capacity_vs_realized.csv")
        for org in e2["organism"].unique():
            sub = e2[e2["organism"] == org].sort_values("u_G")
            margins = sub["margin"].values
            # Margin = V·γ_R - u_G * V·γ_G → strictly decreasing with u_G
            diffs = np.diff(margins)
            # Allow tiny numerical noise
            assert (diffs <= 1e-9).all(), (
                f"{org}: margin not monotonically decreasing with u_G"
            )


# ===========================================================================
# SYSTEM TEST — Reproducibility
# ===========================================================================

class TestReproducibility:
    """Verify that pipeline produces identical results on repeated runs."""

    def test_deterministic_e3(self):
        """Sobol results should be identical across two runs (fixed seed)."""
        from src.identifiability import run_sobol

        df1, _ = run_sobol(organism="yeast", N=64, output="margin")
        df2, _ = run_sobol(organism="yeast", N=64, output="margin")

        # Sort by parameter name for stable comparison (indices with near-equal
        # values may sort differently)
        df1_s = df1.sort_values("parameter").reset_index(drop=True)
        df2_s = df2.sort_values("parameter").reset_index(drop=True)
        pd.testing.assert_frame_equal(df1_s, df2_s)

    def test_deterministic_e1(self):
        """E1 results should be identical across two runs."""
        from src.audit import run_e1

        df1 = run_e1(organism="yeast", n_points=30)
        df2 = run_e1(organism="yeast", n_points=30)

        pd.testing.assert_frame_equal(df1, df2)

    def test_deterministic_e2(self):
        """E2 results should be identical across two runs."""
        from src.audit import run_e2

        df1, s1 = run_e2(organism="yeast", n_points=50)
        df2, s2 = run_e2(organism="yeast", n_points=50)

        pd.testing.assert_frame_equal(df1, df2)
        assert s1["flip_uG"] == s2["flip_uG"]


# ===========================================================================
# SYSTEM TEST — M6 stretch (ecCore)
# ===========================================================================

class TestM6System:
    """System tests for M6 ecCore validation."""

    @pytest.fixture(autouse=True)
    def run_stretch(self):
        subprocess.run(
            [sys.executable, "main.py", "--stretch"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            timeout=120,
        )

    def test_phi_sweep_shows_overflow(self):
        """Tightening Φ should eventually produce acetate overflow."""
        df = pd.read_csv(PROJECT_ROOT / "results" / "m6_phi_sweep.csv")
        optimal = df[df["status"] == "optimal"]
        assert not optimal.empty
        assert optimal["acetate_secretion"].max() > 0.1, (
            "No acetate overflow detected — expected overflow at tight Φ"
        )

    def test_glucose_sweep_shows_fermentation(self):
        """Increasing glucose at tight Φ should trigger fermentation."""
        df = pd.read_csv(PROJECT_ROOT / "results" / "m6_glucose_sweep.csv")
        optimal = df[df["status"] == "optimal"]
        assert not optimal.empty
        assert optimal["frac_fermentative"].max() > 0.05, (
            "No fermentation onset detected in glucose sweep"
        )

    def test_growth_rate_positive(self):
        """Growth rate should be positive for optimal solutions."""
        df = pd.read_csv(PROJECT_ROOT / "results" / "m6_phi_sweep.csv")
        optimal = df[df["status"] == "optimal"]
        assert (optimal["growth_rate"] > 0).all()

    def test_f4_figure_generated(self):
        """F4 ecCore figure should exist after stretch run."""
        assert (PROJECT_ROOT / "figures" / "F4_eccore.png").exists()
        assert (PROJECT_ROOT / "figures" / "F4_eccore.svg").exists()
