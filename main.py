"""Warburg Audit — One-shot reproducible pipeline.

Usage: python main.py [--stretch]

Runs M2→M1→M3→M4→M5 (and M6 if --stretch), writes all artifacts.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from src.params import provenance_table, derived_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Warburg Proteome Efficiency Audit")
    parser.add_argument("--stretch", action="store_true", help="Include M6 ecCore validation")
    args = parser.parse_args()

    t0 = time.time()
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)

    # --- M2: Parameter provenance ---
    print("=" * 70)
    print("[M2] Parameter provenance")
    print("=" * 70)
    prov = provenance_table()
    prov.to_csv(results_dir / "provenance.csv", index=False)
    print(prov.to_string(index=False))
    print()
    derived = derived_table()
    print(derived.to_string(index=False))
    print()

    # --- M3: Audit experiments E1 + E2 ---
    print("=" * 70)
    print("[M3] Audit experiments (E1 + E2)")
    print("=" * 70)
    from src.audit import run_full_audit
    audit_results = run_full_audit()
    print()

    # --- M4: Identifiability (E3) ---
    print("=" * 70)
    print("[M4] Identifiability analysis (E3 — Sobol + Morris)")
    print("=" * 70)
    from src.identifiability import run_e3
    e3_results = run_e3()
    print()

    # --- M5: Visualization ---
    print("=" * 70)
    print("[M5] Generating figures")
    print("=" * 70)
    from src.viz import generate_all_figures
    generate_all_figures()
    print()

    # --- Summary ---
    elapsed = time.time() - t0
    print("=" * 70)
    print(f"PIPELINE COMPLETE in {elapsed:.1f}s")
    print("=" * 70)
    print(f"  Results: {results_dir}/")
    print(f"  Figures: figures/")
    print()

    # Key findings
    print("KEY FINDINGS:")
    for s in audit_results["e2_summaries"]:
        flip = s["flip_uG"]
        if flip:
            print(f"  [{s['organism']}] Verdict flips at u_G = {flip:.4f} "
                  f"(capacity ρ = {s['rho_capacity']:.4f})")
        else:
            print(f"  [{s['organism']}] No verdict flip detected "
                  f"(capacity ρ = {s['rho_capacity']:.4f})")

    sobol_df = e3_results["sobol_margin"]
    print("\n  Dispute-settling measurements (top Sobol ST per organism):")
    for org in sobol_df["organism"].unique():
        top = sobol_df[sobol_df["organism"] == org].iloc[0]
        print(f"    {org}: {top['parameter']} (ST = {top['ST']:.4f})")


if __name__ == "__main__":
    main()
