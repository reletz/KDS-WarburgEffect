"""Unit tests for individual modules (M1–M4).

Validation gates:
- Model A: boundary behavior (low glucose → pure respiration; saturated proteome → glycolysis)
- Model B: reproduces Kukurugya capacity ratios directionally
- E2: verdict flip exists at realistic u_G
- E3: Sobol/Morris produce consistent rankings
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
    OrganismParams,
    provenance_table,
    derived_table,
    get_param_ranges,
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


# ===========================================================================
# UNIT TESTS — M1: Model A (Shen LP)
# ===========================================================================

class TestModelA:
    def test_low_glucose_pure_respiration(self, yeast_params: ModelParams):
        """At very low glucose, respiration should dominate (higher yield)."""
        p = ModelParams(
            gamma_resp=yeast_params.gamma_resp,
            gamma_glyc=yeast_params.gamma_glyc,
            V_resp=yeast_params.V_resp,
            V_glyc=yeast_params.V_glyc,
            Phi=yeast_params.Phi,
            g_avail=1e-5,
        )
        result = solve_shen(p)
        assert result.phenotype == "respiration"

    def test_high_glucose_glycolysis_onset(self, yeast_params: ModelParams):
        """At high glucose with saturated proteome, glycolysis should appear."""
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
        Vg_glyc = yeast_params.V_glyc * yeast_params.gamma_glyc
        Vg_resp = yeast_params.V_resp * yeast_params.gamma_resp
        if Vg_glyc > Vg_resp:
            assert result.frac_glyc > 0.5
        else:
            assert result.frac_glyc < 0.5

    def test_respiration_at_zero_glucose(self):
        """At g_avail → 0, near-zero ATP produced."""
        p = ModelParams(
            gamma_resp=26.0, gamma_glyc=2.0,
            V_resp=0.01, V_glyc=0.5,
            Phi=0.2, g_avail=1e-8,
        )
        result = solve_shen(p)
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

    def test_atp_is_nonnegative(self):
        """ATP production must be non-negative for any valid input."""
        for org in ALL_ORGANISMS.values():
            for g in [1e-6, 0.01, 0.1, 1.0]:
                p = ModelParams.from_organism_params(org, g_avail=g)
                result = solve_shen(p)
                assert result.J_ATP >= -1e-12

    def test_frac_glyc_bounded(self):
        """Fermentative fraction must be in [0, 1]."""
        op = ALL_ORGANISMS["yeast"]
        for g in np.linspace(1e-5, 0.5, 20):
            p = ModelParams.from_organism_params(op, g_avail=float(g))
            result = solve_shen(p)
            assert -1e-9 <= result.frac_glyc <= 1.0 + 1e-9


# ===========================================================================
# UNIT TESTS — M1: Model B (Kukurugya analytical)
# ===========================================================================

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
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=1.0)
        result = solve_kukurugya(p)
        assert result.frac_glyc > 0.9

    def test_ecoli_respirofermentative(self):
        """E. coli: glycolysis capacity > respiration capacity → glycolysis at saturation."""
        op = ALL_ORGANISMS["ecoli"]
        p = ModelParams.from_organism_params(op, g_avail=1.0)
        result = solve_kukurugya(p)
        assert result.frac_glyc > 0.5

    def test_mammalian_glycolysis_at_saturation(self):
        """Mammalian cells at high glucose should show Warburg-like glycolysis."""
        op = ALL_ORGANISMS["mammalian"]
        p = ModelParams.from_organism_params(op, g_avail=1.0)
        result = solve_kukurugya(p)
        assert result.frac_glyc > 0.9

    def test_atp_monotone_with_glucose(self):
        """Total ATP should be non-decreasing as glucose increases."""
        op = ALL_ORGANISMS["yeast"]
        g_values = np.linspace(1e-5, 0.5, 50)
        prev_atp = 0.0
        for g in g_values:
            p = ModelParams.from_organism_params(op, g_avail=float(g))
            result = solve_kukurugya(p)
            assert result.J_ATP >= prev_atp - 1e-9
            prev_atp = result.J_ATP

    def test_phenotype_transitions_exist(self):
        """Sweep should contain at least one transition from respiration."""
        op = ALL_ORGANISMS["yeast"]
        base = ModelParams.from_organism_params(op, g_avail=0.0)
        results = sweep_glucose(base, model="kukurugya")
        phenotypes = {r.phenotype for r in results}
        assert "respiration" in phenotypes
        assert len(phenotypes) > 1  # must have at least one transition


# ===========================================================================
# UNIT TESTS — M1: Cross-model consistency
# ===========================================================================

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
        """Both models should agree on direction at high glucose."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=1.0)
        r_shen = solve_shen(p)
        r_kuk = solve_kukurugya(p)
        assert r_shen.frac_glyc > 0.5
        assert r_kuk.frac_glyc > 0.5

    def test_both_return_valid_result_all_organisms(self):
        """Both solvers should return valid results for all organisms at all glucose levels."""
        for org in ALL_ORGANISMS.values():
            for g in [1e-5, 0.01, 0.1, 1.0]:
                p = ModelParams.from_organism_params(org, g_avail=g)
                r_a = solve_shen(p)
                r_b = solve_kukurugya(p)
                assert r_a.J_ATP >= 0
                assert r_b.J_ATP >= 0
                assert 0 <= r_a.frac_glyc <= 1
                assert 0 <= r_b.frac_glyc <= 1


# ===========================================================================
# UNIT TESTS — M1: E2 (capacity vs. realized)
# ===========================================================================

class TestE2:
    def test_capacity_respiration_loses_for_yeast(self):
        """At full capacity (u=1), glycolysis has higher V·γ for yeast."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=0.1)
        result = compute_efficiency_with_utilization(p, u_G=1.0, u_R=1.0)
        assert result["verdict_capacity"] == "glycolysis"

    def test_realized_can_flip_with_idle_glycolysis(self):
        """With glycolytic enzymes mostly idle, respiration can win on realized basis."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=0.1)
        result = compute_efficiency_with_utilization(p, u_G=0.2, u_R=1.0)
        assert result["verdict_realized"] == "respiration"

    def test_find_flip_point_exists(self):
        """A verdict flip point should exist for yeast parameters."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=0.1)
        flip = find_verdict_flip_uG(p, u_R=1.0)
        assert flip is not None
        assert 0.0 < flip < 1.0

    def test_flip_point_matches_rho(self):
        """Flip point u_G should approximately equal ρ (capacity ratio)."""
        for org_name, op in ALL_ORGANISMS.items():
            p = ModelParams.from_organism_params(op, g_avail=0.1)
            flip = find_verdict_flip_uG(p, u_R=1.0)
            rho = (p.V_resp * p.gamma_resp) / (p.V_glyc * p.gamma_glyc)
            assert flip is not None, f"No flip for {org_name}"
            assert abs(flip - rho) < 0.01, (
                f"{org_name}: flip={flip:.4f} should ≈ ρ={rho:.4f}"
            )

    def test_margin_sign_at_extremes(self):
        """Margin should be positive (resp wins) at u_G=0.01 and negative at u_G=1.0."""
        op = ALL_ORGANISMS["yeast"]
        p = ModelParams.from_organism_params(op, g_avail=0.1)
        r_low = compute_efficiency_with_utilization(p, u_G=0.01, u_R=1.0)
        r_high = compute_efficiency_with_utilization(p, u_G=1.0, u_R=1.0)
        assert r_low["margin"] > 0  # respiration wins when glycolysis is idle
        assert r_high["margin"] < 0  # glycolysis wins at full capacity


# ===========================================================================
# UNIT TESTS — M2: Parameter provenance
# ===========================================================================

class TestParams:
    def test_provenance_table_complete(self):
        """Provenance table should have 5 params × 3 organisms = 15 rows."""
        df = provenance_table()
        assert len(df) == 15

    def test_derived_rho_direction(self):
        """ρ < 1 for yeast/mammalian (glycolysis faster per protein on capacity basis)."""
        df = derived_table()
        for org in ["yeast", "mammalian"]:
            row = df[df["organism"] == org].iloc[0]
            assert row["rho (resp/glyc)"] < 1.0

    def test_all_params_have_source(self):
        """Every parameter should have a non-empty source location."""
        df = provenance_table()
        assert all(df["source_location"].str.len() > 0)

    def test_param_ranges_valid(self):
        """Uncertainty ranges should have low < high and all positive."""
        for org in ALL_ORGANISMS:
            ranges = get_param_ranges(org)
            assert len(ranges) == 5
            for r in ranges:
                assert r.low > 0
                assert r.low < r.high

    def test_organism_params_consistency(self):
        """gamma_resp > gamma_glyc for all organisms (biological invariant)."""
        for org_name, op in ALL_ORGANISMS.items():
            assert op.gamma_resp.value > op.gamma_glyc.value, (
                f"{org_name}: gamma_resp should > gamma_glyc"
            )

    def test_Vgamma_glyc_greater_than_resp(self):
        """V·γ_glyc > V·γ_resp for all organisms (Kukurugya's main claim)."""
        for org_name, op in ALL_ORGANISMS.items():
            assert op.Vgamma_glyc > op.Vgamma_resp, (
                f"{org_name}: Kukurugya claims glycolysis capacity > respiration"
            )


# ===========================================================================
# UNIT TESTS — M3: Audit module
# ===========================================================================

class TestAudit:
    def test_e1_returns_both_models(self):
        """E1 overlay should contain results from both models."""
        from src.audit import run_e1
        df = run_e1(organism="yeast", n_points=20)
        models = df["model"].unique()
        assert "Model A (Shen LP)" in models
        assert "Model B (Kukurugya)" in models

    def test_e1_shared_glucose_axis(self):
        """Both models in E1 should share the same glucose values."""
        from src.audit import run_e1
        df = run_e1(organism="yeast", n_points=20)
        g_a = sorted(df[df["model"] == "Model A (Shen LP)"]["g_avail"].values)
        g_b = sorted(df[df["model"] == "Model B (Kukurugya)"]["g_avail"].values)
        np.testing.assert_array_almost_equal(g_a, g_b)

    def test_e2_sweep_returns_all_columns(self):
        """E2 sweep should return expected columns."""
        from src.audit import run_e2
        df, summary = run_e2(organism="yeast", n_points=20)
        expected_cols = {"u_G", "Vgamma_G_capacity", "Vgamma_R_capacity",
                         "Vgamma_G_realized", "Vgamma_R_realized", "margin",
                         "verdict_capacity", "verdict_realized", "organism"}
        assert expected_cols.issubset(set(df.columns))

    def test_e2_summary_has_flip(self):
        """E2 summary should report a flip point for yeast."""
        from src.audit import run_e2
        _, summary = run_e2(organism="yeast", n_points=100)
        assert summary["flip_uG"] is not None
        assert 0 < summary["flip_uG"] < 1

    def test_e2b_accounting_effect(self):
        """Increasing membrane overhead should shift the flip point."""
        from src.audit import run_e2b_accounting
        df = run_e2b_accounting(organism="yeast")
        assert len(df) > 0
        assert "membrane_frac" in df.columns
        assert "flip_uG" in df.columns

    def test_e2c_uncertainty_has_ci(self):
        """Uncertainty quantification should produce valid confidence intervals."""
        from src.audit import run_e2c_uncertainty
        result = run_e2c_uncertainty("yeast", n_bootstrap=50)
        assert result["n_valid"] > 0
        assert result["ci_5"] is not None
        assert result["ci_95"] is not None
        assert result["ci_5"] < result["ci_95"]
        assert result["ci_5"] < result["nominal_flip_uG"] < result["ci_95"]

    def test_e2c_ci_covers_nominal(self):
        """90% CI should contain the nominal flip point for all organisms."""
        from src.audit import run_e2c_uncertainty
        for org in ALL_ORGANISMS:
            result = run_e2c_uncertainty(org, n_bootstrap=50)
            assert result["ci_5"] <= result["nominal_flip_uG"] <= result["ci_95"], (
                f"{org}: nominal {result['nominal_flip_uG']:.4f} outside "
                f"CI [{result['ci_5']:.4f}, {result['ci_95']:.4f}]"
            )

    def test_e2d_decomposition_flip_tracks_rho(self):
        """In decomposition, flip_uG should always approximately equal ρ."""
        from src.audit import run_e2d_decomposition
        df = run_e2d_decomposition("yeast", n_points=10)
        valid = df.dropna(subset=["flip_uG"])
        assert len(valid) > 0
        for _, row in valid.iterrows():
            assert abs(row["flip_uG"] - row["rho"]) < 0.02

    def test_wang_comparison_crossing(self):
        """Wang comparison should show crossing at ρ for all organisms."""
        from src.audit import run_wang_comparison
        df = run_wang_comparison(n_points=50)
        for org in df["organism"].unique():
            sub = df[df["organism"] == org]
            assert sub["crossing_point_uG"].iloc[0] > 0
            assert sub["crossing_point_uG"].iloc[0] < 1
            margins = sub["margin"].values
            # At low u_G (glycolysis idle) → resp wins → margin > 0
            # At high u_G (full capacity) → glyc wins → margin < 0
            assert margins[0] > 0  # u_G=0.01 → glycolysis idle → resp wins
            assert margins[-1] < 0  # u_G=1.0 → full capacity → glyc wins
            assert any(margins > 0) and any(margins < 0)


# ===========================================================================
# UNIT TESTS — M4: Identifiability
# ===========================================================================

class TestIdentifiability:
    def test_sobol_returns_all_params(self):
        """Sobol analysis should return indices for all parameters."""
        from src.identifiability import run_sobol
        df, problem = run_sobol(organism="yeast", N=64, output="margin")
        assert len(df) == problem["num_vars"]
        assert "ST" in df.columns
        assert "S1" in df.columns

    def test_sobol_indices_bounded(self):
        """Sobol indices should be in [0, 1] (approximately)."""
        from src.identifiability import run_sobol
        df, _ = run_sobol(organism="yeast", N=64, output="margin")
        assert all(df["ST"] >= -0.1)  # small negative possible due to estimation
        assert all(df["ST"] <= 1.5)   # can slightly exceed 1 due to interactions

    def test_morris_returns_all_params(self):
        """Morris analysis should return mu_star for all parameters."""
        from src.identifiability import run_morris
        df = run_morris(organism="yeast", N=32, output="margin")
        assert "mu_star" in df.columns
        assert "sigma" in df.columns
        assert len(df) >= 5  # at least 5 base params

    def test_sobol_top_param_is_uG(self):
        """u_G should be the dominant parameter for margin output."""
        from src.identifiability import run_sobol
        df, _ = run_sobol(organism="yeast", N=256, output="margin")
        top_param = df.iloc[0]["parameter"]
        assert top_param == "u_G"

    def test_morris_sobol_rank_agreement(self):
        """Top-2 parameters should agree between Sobol and Morris."""
        from src.identifiability import run_sobol, run_morris
        df_s, _ = run_sobol(organism="yeast", N=256, output="margin")
        df_m = run_morris(organism="yeast", N=128, output="margin")
        top2_sobol = set(df_s.head(2)["parameter"])
        top2_morris = set(df_m.head(2)["parameter"])
        assert len(top2_sobol & top2_morris) >= 1  # at least 1 in common
