# Computational Biology of Cellular Respiration & ATP Production: Solution Design for Seven Problems

**TL;DR**
- All seven problems are tractable for a 3–4 person undergraduate/Master's team using purely algorithmic/mathematical methods (ODEs, MILP/LP, Gillespie SSA, PDE/Monte Carlo); the highest payoff is in problems 1, 2, 6 (clean reproductions/extensions of published models with public code) and the most ambitious is problem 5 (Monte Carlo on tomography meshes) and problem 7 (genome-scale ETFL-style MILP).
- The single most important reference for each is: (1) Saa & Siqueira 2013 / Bertram et al. 2006; (2) Gillespie 1977 + Cao–Gillespie–Petzold 2006 tau-leaping; (3) Henry et al. 2007 TMFA / Salvy & Hatzimanikatis 2020 ETFL; (4) Selivanov et al. 2009 PLoS CB & Selivanov et al. 2011 PLoS CB; (5) Garcia et al. 2023 J. Gen. Physiol.; (6) Grass et al. 2022 J. Biol. Chem.; (7) Shen et al. 2024 Nat. Chem. Biol. — with the active counter-argument by Kukurugya, Rosset & Titov 2024 PNAS.
- Recommended stack: Python + SciPy/SUNDIALS for ODEs, StochPy/GillesPy2 for SSA, COBRApy + pyTFA + Gurobi/CPLEX academic license for FBA/TFA/ETFL, MCell4 + Blender for spatial — all open-source, all runnable on a single workstation except possibly problem 5 (needs ≥16 GB RAM and overnight runs).

## Key Findings

**Problem-by-problem feasibility summary**

| # | Problem | Core algorithm | Team feasibility | Main complexity bottleneck |
|---|---------|---------------|------------------|----------------------------|
| 1 | BPLS / Magnus–Keizer ODE | RK4(5) / LSODA | High — 4 ODEs | Stiffness in J_ANT, J_uni terms |
| 2 | Gillespie SSA on ATP synthase | Direct SSA + tau-leaping | High | Disparate reaction timescales |
| 3 | TFA/ETFL with TCA | MILP | Medium | MILP scaling with reaction count |
| 4 | Selivanov ROS bistability | Rule-based ODE + continuation | Medium-High | ~400-state ODE generation |
| 5 | Garcia 2023 spatial | MCell4 Monte Carlo on EM mesh | Medium-Low | CPU hours per simulated second |
| 6 | Grass 2022 IRI | Stiff ODE, 67 species | High | Parameter identifiability |
| 7 | Proteome-constrained FBA | ETFL / GECKO MILP | Medium | Genome-scale MILP solve time |

## Details

### 1. ODE-based modeling of mitochondrial ATP homeostasis (BPLS / Magnus–Keizer)

**Biology.** The Magnus–Keizer (MK) model (Magnus & Keizer 1997, 1998, *Am. J. Physiol.*) is a kinetic ODE model for mitochondrial Ca²⁺ handling and oxidative phosphorylation in the β-cell, built on first biophysical principles for the TCA dehydrogenases, the proton pump, F₁F₀-ATPase, ANT, uniporter, and Na⁺/Ca²⁺ antiporter. Bertram, Pedersen, Luciani & Sherman (BPLS, *J. Theor. Biol.* 243:575–586, 2006) reduced MK to **four dynamical variables** — mitochondrial membrane potential ΔΨ, NADHₘ, ADPₘ, and Caₘ — with conservation NADₘ + NADHₘ = NAD_tot and ATPₘ + ADPₘ = A_tot. Saa & Siqueira (*Bull. Math. Biol.* 75:1636, 2013; arXiv:1212.1194) identified that the BPLS approximations for the **ANT flux J_ANT and the uniporter flux J_uni** contained inaccuracies traceable to the fact that BPLS was actually an approximation of Cortassa et al. 2003 rather than the original MK. They proposed corrected sigmoidal/Goldman-Hodgkin-Katz-style flux expressions and showed numerically that the enhanced 4-ODE system possesses a **unique attractor fixed point** over physiological ranges of cytosolic Ca²⁺ (Ca_c) and the glycolytic substrate FBP — consistent with mitochondrial ATP homeostasis.

**Algorithmic solution.** The system is a 4-D, mildly stiff, autonomous nonlinear ODE on the order of milliseconds–seconds. An explicit Runge-Kutta 4(5) (Dormand–Prince) integrator is adequate for the BPLS reduced model; if combining with the full Cortassa 2003 cardiac TCA-coupled extension (~13 ODEs) or with downstream cytosolic Ca²⁺ oscillators, switch to LSODA/CVODE (BDF) because dehydrogenase activations introduce timescale separation between fast (ΔΨ) and slow (Caₘ) variables.

**Computational complexity.** O(N_steps · N_vars) per trajectory with N_vars = 4, completing in <1 s on a laptop; bifurcation/parameter sweeps over (Ca_c, FBP) on a 100×100 grid finish in minutes.

**Extensions feasible without AI.**
1. Bifurcation/continuation analysis using AUTO-07p, PyDSTool's Continuer, or MatCont — verify the unique-attractor claim and search for Hopf bifurcations under elevated FBP.
2. Reintroduce the J_ANT/J_uni "inaccuracies" and reproduce Saa & Siqueira's correction, with side-by-side comparison.
3. Sensitivity analysis (Morris elementary effects + Sobol indices via SALib) over ~20 kinetic parameters.
4. Couple to a glycolytic oscillator (Bertram et al. 2007) to test the original BPLS hypothesis that glycolytic oscillations drive slow Ca²⁺/ATP oscillations.
5. CellML implementation already exists in the Physiome repository (`bertram_2006.cellml`) — students can validate against published figures 3–9.

### 2. Stochastic simulation of ATP synthase kinetics (Gillespie SSA)

**Biology.** Inside a single cristae sub-compartment, copy numbers of ATP synthase (estimated at ~3,070 µm⁻² of high-curvature CM in Garcia et al. 2023), ANT, ADP and ATP are small enough that continuum mass-action breaks down — diffusion-limited gradients of ADP and ATP develop between matrix and cristae (Mannella et al. 2001; Garcia et al. 2019, *Sci. Rep.* 9:18306). The Chemical Master Equation (CME) is then the correct description.

**Algorithmic solution.**
- **Exact:** Gillespie's Direct Method (Gillespie 1977, *J. Phys. Chem.* 81:2340) — sample next-reaction time τ from exponential(a₀) and reaction index from categorical(aⱼ/a₀). O(M) per step where M = # reactions; O(log M) with the Next Reaction Method (Gibson & Bruck 2000).
- **Approximate:** explicit tau-leaping (Gillespie 2001) for systems with many fast firings; OTL (optimized tau-leap, Cao–Gillespie–Petzold 2006) and binomial tau-leap (Chatterjee 2005) handle the "negative population" pathology. Switch automatically when the leap condition is satisfied.
- **CME directly:** Finite State Projection (Munsky & Khammash 2006) is tractable only when the ATP synthase + ADP + ATP joint state space is truncated to a few hundred states.

**Computational complexity.** Direct SSA: O(M · ⟨a₀⟩ · T) total events; for a single cristae with ~100 ATP synthase molecules and turnover ~100 s⁻¹, ~10⁴ events per second simulated — overnight for 1 s of biological time per replicate; embarrassingly parallel across replicates.

**Extensions feasible without AI.**
1. Implement the 6-state Pietrobon–Caplan / Garcia et al. 2023 ATP synthase kinetic diagram in StochPy or GillesPy2 and reproduce the deterministic ODE prediction (well-mixed limit) before adding noise.
2. Compare ⟨ATP(t)⟩ and Var[ATP(t)] from SSA against the chemical Langevin equation (CLE) and the linear noise approximation (van Kampen Ω-expansion) — quantify the noise floor.
3. Test sensitivity of stochastic ATP output to copy number reduction (mimicking Leigh-syndrome–type loss of ATP-synthase dimerization; Siegmund et al. 2018).
4. Use slow-scale SSA (Cao, Gillespie & Petzold 2005) to handle the stiff substrate-enzyme binding step.

### 3. Thermodynamics-based Flux Analysis (TFA/TMFA) and ETFL

**Biology / problem.** Standard FBA can return thermodynamically infeasible loops (futile cycles) — e.g., in the TCA cycle, a reversible succinate-fumarate cycle that hydrolyzes ATP with no net flux. Henry, Broadbelt & Hatzimanikatis (2007, *Biophys. J.* 92:1792, TMFA) introduced linear Gibbs-free-energy constraints with binary variables y_i ∈ {0,1} per reaction such that ΔᵣG'ᵢ + K·yᵢ ≤ K − ε and vᵢ ≤ V_max·(1−yᵢ), forcing flux only in the thermodynamically favorable direction. This is a Mixed-Integer Linear Program. ETFL (Salvy & Hatzimanikatis, *Nat. Commun.* 11, Article 30, 13 January 2020; doi:10.1038/s41467-019-13818-7) extends TFA with **expression coupling** (mRNA→protein→catalytic capacity) and **resource allocation** while keeping the problem MILP rather than MINLP, by discretizing growth rate µ into a binary expansion.

**Algorithmic solution.** Formulate as MILP using COBRApy + pyTFA (Salvy et al. 2019, *Bioinformatics* 35:167) or matTFA. Required inputs: (i) compartment-specific pH, ionic strength, membrane potential; (ii) standard ΔfG° per metabolite from the component contribution method of Noor, Haraldsdóttir, Milo & Fleming (*PLoS Comput. Biol.* 9(7):e1003098, 2013; doi:10.1371/journal.pcbi.1003098) or eQuilibrator; (iii) elemental/charge balance for every reaction. Solve with Gurobi or CPLEX (academic license free).

**Computational complexity.** TFA MILP on the E. coli core model (~95 reactions) solves in milliseconds; iJO1366 (~2,500 reactions) in seconds to minutes; ETFL on E. coli (~5,000 expression variables) in minutes per solve on a modern CPU. The hard step is **variability analysis (TVA)**: 2N LPs/MILPs needed, scaling linearly with reaction count.

**Extensions feasible without AI.**
1. Build a TCA-cycle submodel (~20 reactions) and demonstrate elimination of the succinate↔fumarate futile loop after applying ΔᵣG' constraints.
2. Integrate metabolomics: bound metabolite concentrations between measured ranges and re-solve to tighten the feasible Gibbs energy ranges; replicate the multiTFA reduction reported by Vishal Mahamkali et al. in the multiTFA paper (*Bioinformatics* 2021, btab151; PMC8479682), which states verbatim: "a median reduction of 6.8 kJ/mol in reaction Gibbs free energy ranges, while three out of 12 reactions in glycolysis changed from reversible to irreversible" (specifically ENO, GAPD, and PGM).
3. Couple TFA to a minimal proteome-allocation constraint (foreshadowing problem 7).

### 4. Bistability and ROS in the ETC (Selivanov model)

**Biology.** Selivanov et al. (2009, *PLoS Comput. Biol.* 5:e1000619) modeled the Q-cycle of Complex III using a **rule-based automated DE construction** that enumerates the 400 possible redox states of complex III (combinations of oxidation states of the two cytochrome b hemes, Rieske FeS, cytochrome c₁, Qₒ-site occupant, Qᵢ-site occupant). They demonstrated **bistability**: two stable steady states of Qₒ-bound semiquinone exist for the same parameters, with the choice depending on initial reduction state. A transition to the high-semiquinone state — induced by transient anoxia followed by reoxygenation — explains the paradoxical ROS burst on reperfusion. The follow-up paper (Selivanov et al. 2011, *PLoS Comput. Biol.* 7:e1001115) extended this to the **whole respiratory chain** including Complex I N2-FMN-ubiquinone interactions, proton translocation, ΔΨ, the TCA cycle, and crucially **reverse electron transport (RET)** at Complex I driven by succinate.

**Algorithmic solution.** Rule-based DE generation analogous to BioNetGen (Faeder et al. 2009) — students write ~10 reaction rules and a code-generator emits the hundreds of ODEs. Solve as stiff ODE (BDF / Radau IIA implicit RK). For bifurcation analysis use AUTO-07p or MatCont with continuation in succinate concentration, ΔΨ, or O₂.

**Computational complexity.** N ≈ 400 ODEs for full Complex III; ~600–1000 for the full chain. Stiff Newton iterations dominate — a single trajectory costs seconds, full bistability/hysteresis sweep takes ~1 hour on one core.

**Extensions feasible without AI.**
1. Reproduce Figures 8–9 of Selivanov 2009 (hysteresis loop of Qₒ semiquinone vs succinate) using their published rate constants (Tables 1–3).
2. Add Complex II with the hysteresis/bistability mechanism from the 2019 follow-up (Selivanov et al. *Biochem. Moscow Suppl.* 13:341) — succinate-induced positive feedback at FAD.
3. Quantify ROS production via the Qₒ-semiquinone × O₂ second-order rate and compute the bistability region's sensitivity to pH (Selivanov et al. 2008).
4. Bridge to problem 6: drive the Selivanov ETC model with the ischemic substrate accumulation (succinate ↑, NADH ↑) from the Grass IRI model.

### 5. Spatial constraints, cristae morphology and ATP diffusion (Garcia et al. 2023)

**Biology.** Garcia, Gupta, Bartol, Sejnowski & Rangamani (*J. Gen. Physiol.* 155:e202213263, 2023) built a thermodynamically consistent reaction-diffusion model with six modules (ETC lumped, ATP synthase, ANT based on Metelkin 2006, Pi carrier, leak, buffering) and ran it both as an ODE and as **agent-based Monte Carlo simulations in MCell** on **nine 3-D mitochondrial reconstructions from electron tomography** (Mendelsohn et al. 2021). ATP synthases were localized to high-curvature CM regions (first principal curvature > 70 µm⁻¹) at 3,070 µm⁻². Key result: **ATP production, not ATP export, is the rate-limiting step** under physiological conditions, and **globular mitochondria generate more cytosolic ATP than elongated morphologies** because they have more crista junctions per OM area and more ATP-synthase-bearing high-curvature area. The earlier Garcia 2019 *Sci. Rep.* paper had already demonstrated **sub-organelle ATP gradients** between IMS, matrix, and cytosol — gradients absent in ODE-only models.

**Algorithmic solution.** Two routes:
- **Continuum PDE:** ∂c/∂t = D∇²c + R(c) on the tomographic mesh; finite-element or finite-volume discretization (FEniCS, COMSOL). Reasonable when copy numbers are large.
- **Particle-based Monte Carlo (MCell4 + BioNetGen, Husar et al. 2024 *PLoS Comput. Biol.*):** discrete molecules diffuse via Brownian steps, react on collision with rate-constant-derived probabilities. Mesh built in Blender from EM tomograms (CellBlender pipeline).

**Computational complexity.** MCell scales O(N_particles · N_steps); a single mitochondrion (10⁵–10⁶ particles, Δt ≈ 1 µs) runs ~12–24 CPU-hours per second of biological time. Embarrassingly parallel over replicates (typically 10 per morphology in Garcia 2023).

**Extensions feasible without AI.**
1. Take the 9 published meshes and ATP-synthase placement code from the Rangamani lab GitHub (`RangamaniLabUCSD/spatial_mito_model`) and re-run with altered cristae density or junction radius to quantify how morphology controls cytosolic ATP rate.
2. Compute the curvature-weighted CM area as a single morphological predictor of ATP rate; fit a regression and compare to the published nine-point relation.
3. Compare MCell results against an ODE limit (well-mixed compartments) on the same parameter set — quantify the gradient-induced rate discrepancy reported by Garcia 2019.
4. Implement the Mannella–Pi PDE model (2001, 2013) on the same mesh and compare to particle-based output.

### 6. Ischemia-Reperfusion Injury in cardiomyocytes (Grass et al. 2022)

**Biology.** Grass et al. (*J. Biol. Chem.* 298:101693, 2022) — building on McDougal & Dewey 2017 — built an ODE system spanning **five compartments and 67 molecular species** covering glycolysis, TCA, ETC (Complexes I–IV explicitly), ATP synthase, ANT, mitochondrial buffering, and cytosolic adenylate kinase. The simulation captures: (i) ATP depletion and lactate accumulation during 15 min of ischemia; (ii) a **spike in ΔΨ_m and reverse electron transport at Complex I** at the moment of reperfusion when succinate has accumulated and O₂ returns; (iii) the resulting ROS burst. Their key clinical insight: a **two-step reperfusion protocol with initial O₂ at 5 % of physiological** levels avoids both the ΔΨ spike and RET-driven ROS generation. The succinate–RET mechanism is corroborated by the Krieg/Murphy line of experimental work (Prag, Murphy & Krieg 2023).

**Algorithmic solution.** Stiff nonlinear ODE; use BDF (CVODE/LSODA). Ischemic phase (low O₂) is moderately stiff; reperfusion transient (microseconds for ΔΨ vs minutes for ATP) is severely stiff — fixed-step explicit RK4 will fail; an adaptive implicit method (Radau, BDF) is required.

**Computational complexity.** ~67 ODEs over 30 min biological time integrates in seconds on a laptop. The expensive task is **parameter identification**: ≥200 parameters, ~30 measurable observables (O₂, ATP, NADH, ΔΨ) → use Latin Hypercube + gradient-free optimizer (CMA-ES, scipy.optimize.differential_evolution).

**Extensions feasible without AI.**
1. Reproduce Grass et al. Fig. 4–6 (ATP, ΔΨ_m, ROS trajectories) from their open-source SBML/CellML model on BioModels.
2. Scan O₂ reintroduction profiles parametrically (linear ramp, step, two-step, sinusoidal) and identify the protocol that minimizes the time-integrated RET flux. The Grass paper's headline finding is that initial reperfusion at 5 % O₂ minimizes damage — extend by adding a malonate-mediated SDH inhibition term as predicted by Prag et al. 2023 to test pharmacological cardioprotection in silico.
3. Add cytosolic Ca²⁺ overload via SERCA pump failure and the mPTP (mitochondrial permeability transition pore) opening Hill function (Bertram & Pedersen variants).
4. Couple to the Selivanov Complex I RET kinetics (problem 4) for a mechanistically detailed RET module.

### 7. Proteome-constrained FBA and the Warburg Effect (Shen et al. 2024 Nat. Chem. Biol.)

**Biology.** Shen et al. (*Nat. Chem. Biol.* 20:1123–1132, 2024) challenged the standard explanation of the Warburg effect. The dogma was that aerobic glycolysis is selected because it produces ATP **faster per unit enzyme mass** than respiration (Basan et al. 2015 *Nature* 528:99, Pfeiffer–Schuster–Bonhoeffer 2001 *Science* 292:504). Using ¹³C metabolic flux analysis at genome scale plus quantitative TMTpro proteomics on *S. cerevisiae* (CEN.PK), *Issatchenkia orientalis* SD108, primary CD8⁺ T cells, NCI60 cancer cell lines, and mouse leukemic spleen, they reported (verbatim, bioRxiv 2022.08.10.503479 v2): "respiration is actually several-fold more proteome efficient than glycolysis. Similar results are attained in mammalian tissues and cancer cells." The Nature *News & Views* companion (Nat Chem Biol 20:1108, 2024) summarised it as: "Per unit of enzyme mass, mitochondrial respiration generates energy faster than glycolysis and is thus more proteome efficient." The hypothesis they support instead is **proteome hedging**: cells maintain high glycolytic enzyme expression as spare capacity against future hypoxia, since aerobic-glycolytic yeasts outgrow respiring cells under O₂ limitation but not in O₂-rich aerobic growth.

**Methods used in Shen et al. (NOT ETFL or pcFBA themselves).** Genome-scale ¹³C-MFA via the Maranas group's framework (Gopalakrishnan & Maranas 2015), with custom GEMs (iSace1144 for *S. cerevisiae*, modified Suthers 2020 model for *I. orientalis*); TMTpro 16-plex proteomics (Wühr lab pipeline); a coarse-grained proteome-allocation model with two sectors (glycolytic mass fraction f_G, respiratory mass fraction f_R). ETFL, GECKO and Elsemman whole-cell models are **cited but not run** by the authors. Proteome efficiency is defined operationally as **ε = J_ATP / m_pathway-enzymes** with units mmol ATP · (g protein)⁻¹ · h⁻¹.

**The competing 2024 view.** Kukurugya, Rosset & Titov (*PNAS* 121(46):e2409509121, published November 12, 2024; doi:10.1073/pnas.2409509121) reach the opposite conclusion: their experimental + modeling work reports that **glycolysis produces ATP at 0.54-, 2.1-, and 3.1-fold faster rates per mg pathway protein than respiration in *E. coli*, *S. cerevisiae*, and mammalian cells, respectively** (PNAS Fig. 2 / SI Appendix Fig. S5). Students should treat the Warburg question as actively contested in 2024–2026.

**Algorithmic solution for the student team.** The cleanest student-feasible reproduction is via **ETFL (Salvy & Hatzimanikatis 2020) or GECKO 3.0** on a published GEM. Yeast8 — used as the *S. cerevisiae* base in many ETFL/GECKO studies — contains 3,991 reactions, 1,149 genes, 2,691 compartmentalized metabolites, and 14 compartments (Lu et al. *Nat. Commun.* 2019; yETFL specs, PMC8352978). For mammalian cancer cells use iML1515-style Recon3D or Human-GEM. Add an enzyme-mass constraint Σ m_eᵢ ≤ Φ_max and a per-enzyme catalytic constraint vⱼ ≤ k_catⱼ · [E]ⱼ. Sweep Φ_max and the glucose uptake bound — observe the FBA-predicted switch from pure respiration to aerobic glycolysis. Critically, compare the **proteome efficiency ratio ε_resp/ε_glyc** predicted by ETFL with the Shen et al. measurements and with Kukurugya et al.'s opposing values.

**Computational complexity.** A single ETFL/GECKO MILP on Yeast8 (3,991 reactions, ~1,000 enzyme-constrained reactions after coupling) solves in 30 s–5 min with Gurobi; full nutrient-space sweep (~50 points) takes ~hours. Robust on a single workstation.

**Extensions feasible without AI.**
1. Reproduce Shen et al.'s qualitative finding ε_resp > ε_glyc using GECKO 3.0 (Chen et al. 2023, *Mol. Syst. Biol.*) on Yeast8 — the simplest entry point.
2. Add TFA constraints (problem 3) on top of the enzyme constraints to get the full ETFL formulation; check whether thermodynamic infeasibility removes any of the predicted Warburg flux distributions.
3. Implement the Shlomi et al. 2011 / Vazquez 2010 solvent-capacity (macromolecular crowding) constraint as an alternative explanation and compare.
4. Build the Shen et al. 2-sector toy model (f_G + f_R + f_other = 1, maximize µ subject to ATP demand) analytically and reproduce their Fig. 6 phase diagram — this is the most elegant didactic exercise.
5. Re-run the same analysis under the Kukurugya & Titov 2024 parameterisation and quantify when the ETFL prediction crosses the ε_resp = ε_glyc boundary.

## Recommendations

**Stage 1 (weeks 1–4) — feasibility & infrastructure:** Pick 2 of the 7 problems per pair of students. Recommended pairings — (1,6) for ODE-focused students; (3,7) for optimization/MILP students; (2,5) for stochastic/spatial students; (4) as an ambitious individual project that bridges. Set up: Python ≥3.10, SciPy, CVODE via SUNDIALS or scipy.integrate.LSODA, COBRApy, pyTFA, GillesPy2, MCell4. Get a free academic Gurobi license.

**Stage 2 (weeks 5–10) — baseline reproduction:** For each problem, **reproduce one published figure** before any extension. Specifically: (1) BPLS Fig. 3 steady-state vs Ca_c; (2) Gillespie SSA on a 6-state ATP-synthase reproducing Garcia 2023 ATP-flux; (3) TMFA on E. coli core (Henry 2007 Fig. 2); (4) Selivanov 2009 Fig. 8 bistability hysteresis; (5) Garcia 2023 Fig. 8C globular > elongated ATP rate; (6) Grass 2022 Fig. 4 ATP/ΔΨ during ischemia; (7) Shen 2024 Fig. 6 phase diagram f_G/f_R.

**Stage 3 (weeks 11–14) — one extension each:** From the extension lists above, each pair commits to one quantitative extension. The benchmark for "done" is a falsifiable comparison against an experimental or published numerical result. **Threshold to escalate scope:** if baseline reproduction is not within 20 % of published values by week 8, abandon extensions and write up the reproduction itself.

**Software/data sources (all open):** Bertram CellML at Physiome (`bertram_2006.cellml`), Selivanov models on BioModels, Garcia code at `RangamaniLabUCSD/spatial_mito_model`, Salvy ETFL at `EPFL-LCSB/etfl`, multiTFA at the relevant PMC repo, Grass open-source SBML in the JBC 2022 supplement, Shen 13C-MFA code at `maranasgroup/iSace_GSM`, `maranasgroup/SteadyState-MFA`, `maranasgroup/yeastsMFA`, and `yihuishen/T_cell_MFA`.

## Caveats

- The "BPLS inaccuracies in J_ANT and J_uni" identified by Saa & Siqueira 2013 are not universally accepted; the original BPLS authors did not publish a rebuttal, but the fact that BPLS was an approximation of the Cortassa 2003 cardiac model rather than the original Magnus–Keizer model is a documented historical wrinkle. Students should present both formulations and let the reader judge.
- For problem 7, the Shen et al. 2024 paper does NOT itself use ETFL or pcFBA — it uses empirical proteome-efficiency calculation (flux ÷ enzyme mass) plus a coarse-grained 2-sector model. ETFL/GECKO are appropriate *extension* tools but should not be presented as the original methodology.
- The competing Kukurugya, Rosset & Titov 2024 *PNAS* paper (glycolysis 0.54×–3.1× faster per mg protein than respiration depending on organism) reaches the **opposite** conclusion to Shen et al. using a different experimental + theoretical model; the field is currently unsettled and students should flag this controversy explicitly rather than treat Shen et al. as settled.
- Garcia 2023's claim that "globular mitochondria produce more ATP than elongated ones" is from simulations over nine reconstructions — a small sample, with substantial geometric variability between mitochondria; the morphology–ATP relationship is correlative within the simulation, not causally proven in vivo.
- Selivanov bistability has been reproduced computationally but direct experimental demonstration of two stable Qₒ-semiquinone steady states in intact mitochondria remains indirect.
- Grass et al.'s "5 % O₂ initial reperfusion" recommendation is a model prediction that has not yet been validated in animal models of myocardial infarction; cite it as in-silico-only.
- MCell4 simulations of full mitochondria typically diverge from continuum predictions at low ANT copy numbers (≤10 per cristae junction); the regime where students would learn the most is exactly the regime most expensive to compute.
- The "several-fold" proteome-efficiency advantage for respiration in Shen et al. 2024 is the headline phrasing in both the bioRxiv preprint (2022.08.10.503479) and the published Nat. Chem. Biol. paper; exact per-system fold-values are reported only in the figure bars and supplementary tables and are not summarized as a single global multiplier in the paper text.