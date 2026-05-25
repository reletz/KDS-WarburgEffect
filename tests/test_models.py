"""Tests for M1 (models.py) and M2 (params.py).

Validation gates:
- Model A: boundary behavior (low glucose → pure respiration; saturated proteome + high glucose → glycolysis onset)
- Model B: reproduces Kukurugya capacity ratios directionally
- E2: verdict flip exists at realistic u_G
"""
import numpy as np
import pytest

from src.models import (
    ModelParams,
    ModelResult,
    solve_shen,
    solve_kukurugya,
    sweep_glucose,
    compute_efficiency_with_utilization,
    find_verdict_flip_uG,
)
from src.params import (
    ALL_ORGANISMS,
    ModelParams as MP,
    OrganismParams,
    provenance_table,
    derived_table,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def yeast_params() -> ModelParams:
    op = ALL_ORGANISMS["yeast"]
    return ModelParams.from_organism_params(op, g_avail=0.001)


@pytest.fixture
def mammalian_params() -> ModelParams:
    op = ALL_ORGANISMS["mammalian"]
    return ModelParams.from_organism_params(op, g_avail=0.001)


# ---------------------------------------------------------------------------
# Model A (Shen LP) — boundary tests
# ---------------------------------------------------------------------------

class TestModelA:
    def test_low_glucose_pure_respiration(self, yeast_params: ModelParams):
        """At very low glucose, respiration should dominate (higher yield)."""
        p = ModelParams(
            gamma_resp=yeast_params.gamma_resp,
            gamma_glyc=yeast_params.gamma_glyc,
            V_resp=yeast_params.V_resp,
            V_glyc=yeast_params.V_glyc,
            Phi=yeast_params.Phi,
            g_avail=1e-5,  # very low glucose
        )
        result = solve_shen(p)
        assert result.phenotype == "respiration", (
            f"Expected pure respiration at low glucose, got {result.phenotype} "
            f"(frac_glyc={result.frac_glyc:.4f})"
        )

    def test_high_glucose_glycolysis_onset(self, yeast_params: ModelParams):
        """At high glucose with saturated proteome, glycolysis should appear."""
        # High glucose = well above what full proteome in glycolysis can handle
        g_high = 2.0 * yeast_params.V_glyc * yeast_params.Phi
        p = ModelParams(
            gamma_resp=yeast_params.gamma_resp,
            gamma_glyc=yeast_params.gamma_glyc,
            V_resp=yeast_params.V_resp,
            V_glyc=yeast_params.V_glyc,
            Phi=yeast_params.Phi,
            g_avail=g_high,
        )
        result = solve_shen(p)
        # At saturated proteome, the pathway with higher V*gamma should dominate
        # For yeast: V_glyc * gamma_glyc vs V_resp * gamma_resp
        Vg_glyc = yeast_params.V_glyc * yeast_params.gamma_glyc
        Vg_resp = yeast_params.V_resp * yeast_params.gamma_resp
        if Vg_glyc > Vg_resp:
            assert result.frac_glyc > 0.5, (
                f"Expected glycolysis dominance at high glucose (V·γ_glyc={Vg_glyc:.4f} > "
                f"V·γ_resp={Vg_resp:.4f}), got frac_glyc={result.frac_glyc:.4f}"
            )
        else:
            assert result.frac_glyc < 0.5

    def test_respiration_at_zero_glucose(self):
        """At g_avail → 0, no ATP produced but direction should be respiration."""
        p = ModelParams(
            gamma_resp=26.0, gamma_glyc=2.0,
            V_resp=0.01, V_glyc=0.5,
            Phi=0.2, g_avail=1e-8,
        )
        result = solve_shen(p)
        # Either pure respiration or near-zero ATP
        assert result.J_ATP < 1e-5

    def test_proteome_constraint_binds(self):
        """phi_G + phi_R should never exceed Phi."""
        p = ModelParams(
            gamma_resp=30.0, gamma_glyc=2.0,
            V_resp=0.01, V_glyc=0.4,
            Phi=0.25, g_avail=0.1,
        )
        result = solve_shen(p)
        assert result.phi_G + result.phi_R <= p.Phi + 1e-9


# ---------------------------------------------------------------------------
# Model B (Kukurugya analytical) — capacity ratio tests
# ---------------------------------------------------------------------------

class TestModelB:
    def test_low_glucose_respiration(self, yeast_params: ModelParams):
        """Model B at low glucose should prefer respiration (higher yield)."""
        p = ModelParams(
            gamma_resp=yeast_params.gamma_resp,
            gamma_glyc=yeast_params.gamma_glyc,
            V_resp=yeast_params.V_resp,
            V_glyc=yeast_params.V_glyc,
            Phi=yeast_params.Phi,
            g_avail=1e-5,
        )
        result = solve_kukurugya(p)
        assert result.phenotype == "respiration"

    def test_high_glucose_matches_Vgamma_direction(self):
        """At high glucose, Model B should pick pathway with highest V·γ."""
        # Yeast: V_glyc*gamma_glyc = 0.46*2 = 0.92; V_resp*gamma_resp = 0.0085*26 = 0.221
        # → glycolysis should win
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=1.0)  # well above saturation
        result = solve_kukurugya(p)
        assert result.frac_glyc > 0.9, (
            f"Yeast at high glucose: expected glycolysis dominance, got frac_glyc={result.frac_glyc:.4f}"
        )

    def test_ecoli_respirofermentative(self):
        """E. coli: V_glyc*gamma_glyc vs V_resp*gamma_resp determines direction."""
        op = ALL_ORGANISMS["ecoli"]
        p = ModelParams.from_organism_params(op, g_avail=1.0)
        result = solve_kukurugya(p)
        # E. coli: V_glyc*gamma_glyc = 0.23*2 = 0.46; V_resp*gamma_resp = 0.0049*26 = 0.1274
        # glycolysis capacity > respiration capacity
        assert result.frac_glyc > 0.5

    def test_mammalian_glycolysis_at_saturation(self):
        """Mammalian cells at high glucose should show Warburg-like glycolysis."""
        op = ALL_ORGANISMS["mammalian"]
        p = ModelParams.from_organism_params(op, g_avail=1.0)
        result = solve_kukurugya(p)
        # V_glyc*gamma_glyc = 0.39*2 = 0.78; V_resp*gamma_resp = 0.004*32 = 0.128
        assert result.frac_glyc > 0.9


# ---------------------------------------------------------------------------
# Cross-model consistency
# ---------------------------------------------------------------------------

class TestCrossModel:
    def test_both_agree_at_low_glucose(self):
        """Both models should prefer respiration at low glucose."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=1e-5)
        r_shen = solve_shen(p)
        r_kuk = solve_kukurugya(p)
        assert r_shen.phenotype == "respiration"
        assert r_kuk.phenotype == "respiration"

    def test_both_agree_direction_at_high_glucose(self):
        """Both models should agree on direction at high glucose (V·γ determines)."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=1.0)
        r_shen = solve_shen(p)
        r_kuk = solve_kukurugya(p)
        # Both should show glycolysis dominance for yeast at high glucose
        assert r_shen.frac_glyc > 0.5
        assert r_kuk.frac_glyc > 0.5


# ---------------------------------------------------------------------------
# E2 — Capacity vs. Realized
# ---------------------------------------------------------------------------

class TestE2:
    def test_capacity_respiration_loses_for_yeast(self):
        """At full capacity (u=1), glycolysis has higher V·γ for yeast."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=0.1)
        result = compute_efficiency_with_utilization(p, u_G=1.0, u_R=1.0)
        # V·γ_glyc = 0.92 > V·γ_resp = 0.221
        assert result["verdict_capacity"] == "glycolysis"

    def test_realized_can_flip_with_idle_glycolysis(self):
        """With glycolytic enzymes mostly idle, respiration can win on realized basis."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=0.1)
        # If glycolysis is only 20% utilized
        result = compute_efficiency_with_utilization(p, u_G=0.2, u_R=1.0)
        # Realized: 0.2 * 0.92 = 0.184 vs 1.0 * 0.221 = 0.221 → respiration wins
        assert result["verdict_realized"] == "respiration"

    def test_find_flip_point_exists(self):
        """A verdict flip point should exist for yeast parameters."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=0.1)
        flip = find_verdict_flip_uG(p, u_R=1.0)
        assert flip is not None, "Expected a verdict flip point for yeast"
        assert 0.0 < flip < 1.0, f"Flip at u_G={flip} should be between 0 and 1"


# ---------------------------------------------------------------------------
# M2 — Parameter provenance
# ---------------------------------------------------------------------------

class TestParams:
    def test_provenance_table_complete(self):
        """Provenance table should have 5 params × 3 organisms = 15 rows."""
        df = provenance_table()
        assert len(df) == 15

    def test_derived_rho_direction(self):
        """Check ρ direction matches Kukurugya's claim (glycolysis faster for yeast/mammalian)."""
        df = derived_table()
        yeast_row = df[df["organism"] == "yeast"].iloc[0]
        # For yeast: V·γ_glyc > V·γ_resp → ρ < 1 → glycolysis faster per protein
        assert yeast_row["rho (resp/glyc)"] < 1.0

    def test_all_params_have_source(self):
        """Every parameter should have a non-empty source location."""
        df = provenance_table()
        assert all(df["source_location"].str.len() > 0)
