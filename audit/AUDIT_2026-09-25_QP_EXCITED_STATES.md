# QDEX audit, second pass: quasiparticle models, excited-state frameworks, documentation

Date: 25 September 2026. Baseline: branch `claude/gallant-wright-6u105o` at
`9ab7b9b` (post-audit fixes `8188931`, CdSe test `c100b68`/`9ab7b9b`).
Scope: scientific soundness of the QP and BSE layers against the literature,
consistency of documentation and code, and a numerical study of the supplied
CdSe 2.0 nm cluster. It builds on `SCIENTIFIC_AUDIT.md` and
`POST_FIX_VERIFICATION.md` and does not repeat their derivations.

## 0. Summary

**What is sound.** The TDA BSE core is correctly assembled:
`A = ΔE_QP + 2K^x(bare) − K^d(W)` for singlets, no exchange for triplets, and a
factor-1 exchange in the spinor basis. The complex-conjugation pattern keeps
the matrix Hermitian, and the length-gauge oscillator strengths use the right
spin prefactors. The Löwdin/Mulliken transition-charge machinery, dense and
matrix-free paths and the post-audit P0 fixes (QP eigenvalues reach the BSE;
diagonal-mode exchange; fraction normalisation) are verified on the CdSe run.
All 100 unit tests pass.

**What limits predictive power.** The quality of an optical gap from QDEX is
controlled by two model choices whose spread is far larger than any
numerical error:

1. **Charged vs. neutral environment inconsistency.** The QP models add a
   dielectric-mismatch (surface polarization) term, but the Resta, DIM and
   `xs` BSE direct kernels have no `eps_out` dependence. The CdSe run shows the
   consequence directly: the exciton binding stays at 0.227 eV while the QP gap
   moves by 0.66 eV between vacuum and `eps_out = 2.4`, so the optical gap
   inherits the full solvent shift. Physically the neutral excitation is
   nearly insensitive to the environment (Brus 1984; Delerue, Lannoo, Allan
   2003; Hybertsen–Louie-type image arguments). The docs now state this
   (`validation/environment_cancellation`), but the code does not yet fix it.
2. **Screening model for W.** On the same QP gap (3.861 eV) the CdSe binding
   energy is 0.23 eV (Resta), 0.32 eV (DIM), 1.38 eV (sBSE monopole RPA) and
   1.68 eV (bare MNOK). There is no internal criterion to choose among them.
   The Resta profile uses bulk ε∞ everywhere (no surface), sBSE is almost
   unscreened (ε_eff ≈ 1.06), and the "Penn" size correction halves ε for this
   dot. That is roughly twice the reduction found in atomistic calculations.

**Bugs fixed in this pass.** `--qp_gap brus`, the CLI default, crashed on every
run (argument order). The `dynamic_z` docstring and CLI help, and the
`evgw`/`qsgw`/`sgw` help strings, overclaimed. The sBSE citation page number
was wrong in code (3054 → 3438). The Davidson solver stalled just above its
default tolerance on every large CdSe problem and aborted (fixed; §4 item 6). Several documentation
contradictions with the code are corrected (Penn formula, MNOK mixing rule).
See §5 for the full list.

## 1. Status of the first audit's findings

| First-audit item | Status now | Evidence |
|---|---|---|
| P0 qsgw eigenvalues discarded | **Fixed** | `cli.py` passes `eps_qp_active`; qsgw-dim CdSe run: nonrigid QP levels enter BSE |
| P0 qsgw zero-ΔW loses bulk shift | **Fixed** | `f_h = f_l = 0.5` fallback (`hardness.py` qsgw loops) |
| P0 signed solvent destroys split | **Fixed** (fractions clipped to [0,1]) | code inspection |
| P0 diagonal BSE drops exchange | **Fixed** | diag_bse CdSe run: Kx present (S₁ 3.653 vs full 3.634 eV) |
| P1 diagonal ignores Löwdin | **Fixed** | code inspection |
| P1 complex population phase | **Fixed** (`abs()**2`) | code; lowdin diagonal exchange still takes `.real` of complex transition charges (only matters for complex spatial MOs) |
| P1 qsgw needs square C | **Fixed** (explicit `ValueError`) | code |
| P1 `xs` "exact" claim | **Docs fixed; code comments still say "exact"** (`build_xs_kernel` docstring, `hardness.py:2139,2188`) | open, cosmetic |
| P1 sbse QP site-diagonal, own bulk reference | **Open** | unchanged |
| P1 `dynamic_z` mislabelled | **Docs fixed earlier; code docstring and CLI help fixed now** | this pass |
| P1 Approach B uses monomer data for absolute edges | **Labelled** (`anchor-reconstructed`) but still reported as IP/EA | open |
| P1 Wannier/bulk convergence untested | **Open**; CdSe active-space study below | §3.3 |
| P2 scale claims | Docs softened | — |
| P2 silent non-convergence | **Partly**: qsgw returns `converged` flag only; Davidson raises and discards all results | open (§4.4) |
| Elementwise `max(0, W_QD − W_bulk)` clipping | **Open** in sgw-resta, qsgw-dim, qsgw-resta | code |
| `eps_out` silently clamped to ≥ 1 | **Open** (`max(1.0, eps_out)` in qsgw) | code |
| Penn formula mismatch | **Docs corrected now** to state the implemented expression; the model issue remains (§2.3) | this pass |
| MNOK mixing rule | **Docs corrected now** (page stated both rules) | this pass |

## 2. Scientific assessment against the literature

### 2.1 What the QP layer is, and how it should be described

None of the `qp_gap` options computes a self-energy. They split into three kinds:

* **Scalar gap interpolation:** `pbe`, `brus`, `gw`/`sgw-anchor`.
* **Classical charging functionals** `½ qᵀΔW q` added to a bulk ΔGW opening:
  `sgw-dim`, `sgw-resta(-pure)`, `sgw` (sBSE).
* **Model fixed points:** `evgw-*` iterates a scalar gap against Penn-scaled
  screening. `qsgw-*` is a static ΔCOHSEX, Hadamard (atom-block ZDO)
  orbital-relaxation model in the Löwdin AO basis.

That is a legitimate and useful design for 10³–10⁴-atom nanocrystals,
similar in spirit to Delerue–Lannoo–Allan's self-energy/polarization
decomposition. It should be sold as **"bulk-calibrated GW with a
finite-size polarization correction"**, not as G₀W₀/evGW/qsGW. The docs now
mostly say this; the CLI option names (`evgw`, `qsgw`) still suggest more.
I recommend keeping the names for compatibility but always printing the
honest label in the provenance, which the code largely does.

Physics checks:

* **Brus model** (`brus`) uses the experimental bulk gap plus the kinetic
  confinement term only. It omits the −1.786 e²/(εR) Coulomb term and the
  polarization terms, so it is neither a QP gap nor an optical gap. It also
  switches to a hyperbolic band for any bulk gap < 2 eV, which includes CdSe.
  Label it "EMA kinetic confinement estimate".
* **Anchor model.** κ = 11.52 eV Å (1 − 1/ε) is 0.8 × e²(1 − 1/ε). The factor
  0.8 is an undocumented averaging choice; Brus' self-polarization for a
  1S envelope gives ≈ 0.79–0.9 depending on ε. Document the origin.
  The anchor residual `A (R₀/R)²` is fixed from the **vacuum** anchor and is
  kept unchanged in solvent. That is an assumption: the residual may itself
  contain polarization.
* **ΔCOHSEX in qsgw-*.** For a Löwdin atom-block ΔW, `Σ_SEX = −½ P∘ΔW` and
  `Σ_COH = ½ diag ΔW` are a consistent static COHSEX *difference* within the
  ZDO assumption. The damping constant Z is fixed at 0.80 when `dynamic_z` is
  off, while the CLI forces `dynamic_z` on for `evgw`/`qsgw`. The CdSe run gives
  Z ≈ 0.94. The effect of Z on the gap is therefore a hidden ~5–15 % knob.

### 2.2 The BSE layer

The TDA is appropriate for the lowest excitons of II–VI dots. TDA errors are
typically ≤ 0.05–0.1 eV for bright band-edge states in closed-shell
semiconductor clusters. The singlet/triplet/spinor factors are correct.
One statement in `excitons/foundations` was wrong (now corrected): triplets lack
K^x because the spin-flip combination makes the e–h *annihilation* matrix
element vanish, not because of "spatial exchange cancellation".

The diagonal-BSE option recovers 98 % of the coupled S₁ shift for this cluster
(S₁ 3.653 vs 3.634 eV). It redistributes oscillator strength, though:
f(S₁) = 0.71 vs 0.39, because the three near-degenerate bright band-edge
states are not mixed. It is fine for peak positions near the edge, not for
line shapes or fine structure.

### 2.3 Screening: where most of the physics error sits

* **Resta profile.** `S(r) = 1/ε + (1 − 1/ε) e^{−k_s r}` with
  `k_s = √(ε − 1)/d_NN` gives λ_s = 1.14 Å for CdSe. It is a Yukawa-type
  interpolation, not Resta's (1977) Thomas–Fermi construction, where
  ε(r) is flat beyond a screening radius R_s fixed by
  sinh(q R_s)/(q R_s) = ε. For the BSE it uses bulk ε∞ at every distance and has
  no surface term. For a 2 nm dot this double-counts interior screening and
  misses surface polarization. Both errors push the binding down.
* **"Penn" scaling.** The implemented
  `ε = 1 + (ε∞ − 1)/[1 + (ΔE/E_g^{PBE})²]` is not Penn's model. Penn gives
  (ε − 1) ∝ (ħω_p/E_P)², with E_P the average gap (≈ 4–5 eV in CdSe). Using the
  0.64 eV PBE band gap as reference gives ε = 2.98 for this dot, against ≈ 4.7
  from the Penn form. Atomistic calculations (Wang–Zunger 1994;
  Delerue–Lannoo–Allan 2003) find reductions of tens of percent at 2 nm. I
  recommend `ε_in − 1 = (ε∞ − 1)[E_P/(E_P + ΔE_conf)]²` with a tabulated E_P,
  or directly the Wang–Zunger form `ε(R) = 1 + (ε∞ − 1)/(1 + (α/R)^l)` with
  material fits.
* **sBSE (Cho–Bintrim–Berkelbach 2022).** The atom-condensed transition
  charges keep only monopoles. Intra-atomic (s–p, p–d) polarization, which
  dominates the electronic polarizability of CdSe, is lost. With 664,000
  transitions the kernel still reports ε_eff ≈ 1.06 for the band-edge exciton
  and a median inter-site ε of 0.99 (< 1). This is a basis-condensation
  artefact, not a converged RPA. The AO-resolved mode or an added
  induced-dipole (charge–dipole) term would be needed. Until then, sBSE
  bindings of > 1 eV for CdSe dots should be treated as unscreened upper bounds.
* **DIM/Thole.** The response is computed properly as an induced-dipole
  problem, but mapped back to W through a normalized scalar per atom (first
  audit, §4). The binding (0.32 eV) lies between Resta and sBSE.

**Recommended direction (single consistent W).** Build one atomistic
charge–dipole response, with monopole charge transfer plus Thole dipoles,
inside a dielectric continuum ε_out. The reaction field can use the image-charge
Green's function of a sphere or a PCM-like surface discretization. Then use the same W in
the QP charging term and in K^d. This gives the Brus/Delerue cancellation by
construction and removes the need for separate "Penn" and "solvent" patches.

### 2.4 What the CdSe run says physically

With `eps_out = 2.4` (a toluene/ligand-shell-like environment), the `gw` model
with the Resta BSE gives S₁ = 2.98 eV. The DIM and sgw-resta QP models give
similar values. This agrees with the experimental first exciton of a
d ≈ 1.6–2.0 nm CdSe dot, about 2.7–3.0 eV. That reference comes from the
Yu–Qu–Guo–Peng (2003) sizing curve; the more recent Aubert–Hens sizing
function (Nano Lett. 22, 1778 (2022)) should replace it once its CdSe numbers
are added to the benchmark page.
The agreement comes partly from cancelling errors, though. The binding (0.23 eV) is
about half the Brus estimate 1.786 e²/(ε∞R) ≈ 0.45 eV for this size. The QP gap
compensates through a solvent term that the BSE does not see. In vacuum the
same model predicts 3.63 eV, which is not a physical prediction. Conclusion:
**QDEX currently reproduces optical gaps only at the ε_out it was tuned
with**, and QP gaps and binding energies individually are not validated.

## 3. CdSe 2.0 nm numerical study

System: Cd₆₈Se₅₅Cl₂₆ (149 atoms, 2753 AOs, DZVP-MOLOPT-PBE-GTH), full MO set,
DFT gap 1.4568 eV, R_eff = 9.221 Å (hull + 1.25 Å). Spin-free unless noted;
25 occupied × 25 virtual active space; dense diagonalization. Timings on 4 CPU
cores. Inputs in `tests/CdSe/`; runs reproducible with the commands below.

### 3.1 Quasiparticle models (Resta BSE kernel, ε_out = 1 unless noted)

| qp_gap | QP gap (eV) | S₁ (eV) | f(S₁) | E_b = QP − S₁ | wall |
|---|---|---|---|---|---|
| pbe | 1.457 | 1.230 | 0.13 | 0.227 | 7 s |
| brus | *crashed before fix* | — | — | — | — |
| gw (anchor) | 3.861 | 3.634 | 0.39 | 0.227 | 8 s |
| gw, ε_out = 2.4 | 3.203 | 2.977 | 0.32 | 0.226 | 7 s |
| gw + `--dynamic_z` | 3.861 | 3.634 | 0.39 | 0.227 | 7 s (flag silently ignored) |
| sgw-resta | 3.762 | 3.535 | 0.38 | 0.227 | 9 s |
| sgw-resta, ε_out = 2.4 | 3.226 | 2.999 | 0.32 | 0.227 | 10 s |
| sgw-resta-pure | 3.498 | 3.271 | 0.35 | 0.227 | 10 s |
| sgw-dim | 3.593 | 3.366 | 0.36 | 0.227 | 9 s |
| evgw-resta | 3.851 | 3.625 | 0.39 | 0.226 | 34 s |
| qsgw-dim | 3.955 | 3.728 | 0.44 | 0.227 | 16 s (3 iterations) |
| qsgw-resta | 4.234 | 4.005 | 0.47 | 0.229 | 49 s (20 iterations) |
| sgw + sbse kernel | 4.598 | 3.219 | 0.41 | 1.379 | 17 s |

Spread of the vacuum QP gap across the "GW-like" models: 3.50–4.60 eV.
With ε_out = 2.4 the models converge to 3.20–3.23 eV.

### 3.2 BSE kernel and framework (gw QP gap 3.861 eV, ε_out = 1)

| variant | S₁ (eV) | f(S₁) | E_b (eV) |
|---|---|---|---|
| Resta (default) | 3.634 | 0.39 | 0.227 |
| DIM | 3.540 | 0.41 | 0.321 |
| sBSE (atom) | 2.482 | 0.32 | 1.379 |
| bare MNOK (`--kernel bse`) | 2.184 | 0.24 | 1.677 |
| `xs-resta` (AO density-pair integrals + Resta) | 3.594 | 0.35 | 0.267 |
| Resta, Löwdin charges | 3.629 | 0.41 | 0.232 |
| triplet (Resta) | T₁ = 3.578 | 0 | 0.283 (S₁ − T₁ = 56 meV, spin-free) |
| diagonal_bse | 3.653 | 0.71 | 0.208 |
| independent_qp | 3.861 | 0.75 | 0 |

### 3.3 Active-space convergence (Resta, gw)

| n_occ × n_virt | transitions | S₁ (eV) | E_b (eV) | wall |
|---|---|---|---|---|
| 25 × 25 | 625 | 3.634 | 0.227 | 8 s |
| 50 × 50 | 2 500 | 3.626 | 0.235 | 9 s |
| 100 × 100 | 10 000 | 3.615 | 0.246 | 212 s (dense) |
| 100 × 100 (Davidson, fixed) | 10 000 | 3.615 | 0.246 | 17 s, 26 iterations (identical to dense) |
| 200 × 200 (Davidson, fixed) | 40 000 | 3.606 | 0.255 | 98 s, 31 iterations; before the fix: aborted after 500 iterations |

The binding grows by ~8–10 meV per doubling of the window (0.227 → 0.255 eV
from 25 × 25 to 200 × 200) and is not converged at the default 25 × 25.
For CdSe dots of this size, use ≥ 100 × 100 with Davidson for binding energies.

### 3.4 Reproduce

```bash
gunzip -k tests/CdSe/MOs_cleaned_20ang.txt.gz
QDEX_RUN_CDSE=1 python -m pytest tests/test_cdse_integration.py   # reference values
cd tests/CdSe && qdex --config config.yaml [--qp_gap ... --eps-out ... --kernel ...]
```

## 4. Implementation findings (new in this pass)

| # | Severity | Finding | Action |
|---|---|---|---|
| 1 | **High** | `--qp_gap brus` (CLI default) always crashed: `estimate_brus_qp_gap(coords, syms, material)` vs signature `(material, coords, symbols)` | **Fixed** (keyword call) + opt-in regression test `tests/test_cdse_integration.py` |
| 2 | High (physics) | BSE direct kernels ignore `eps_out` while QP models include it; optical gap shifts 0.66 eV between vacuum and ε_out = 2.4 in CdSe | Documented with figure; needs a common W (§2.3) |
| 3 | High (physics) | "Penn" scaling uses the PBE band gap as the Penn gap; ε_in(2 nm CdSe) = 2.98 | Docs now state the actual formula and the issue; model change recommended |
| 4 | High (physics) | sBSE monopole RPA gives ε_eff ≈ 1 (median inter-site 0.99 < 1) for CdSe | Report; use AO mode or add dipoles; warn when ε_eff < 1.5 for a semiconductor |
| 5 | Medium | `--dynamic_z` has no effect for `gw`/`sgw-anchor` and is not reported as ignored | Warn or reject |
| 6 | **High** | Davidson could not converge any large CdSe problem. It dropped correction vectors with an *absolute* norm threshold of 10⁻⁵ after orthogonalization. Once residuals approached tol = 10⁻⁵ the subspace stopped growing and the solver stalled at residual ≈ 10⁻⁵ for hundreds of iterations, then raised and discarded all roots (100×100 and 200×200 CdSe runs both aborted). It also started from random complex vectors. | **Fixed**: normalize before orthogonalizing, two-pass Gram–Schmidt, relative threshold, lowest-diagonal initial guess; `tests/test_davidson.py` fails on the old code and passes now. Still recommended: thick restart, real arithmetic for real problems, return unconverged roots with a flag |
| 7 | Medium | YAML key `physics.exchange: true` means *direct e–h* (`include_direct_eh`), not exchange; the supplied CdSe config uses it | Rename in examples; deprecation already printed |
| 8 | Low | "Confinement Energy" printed as QP gap − experimental bulk gap; this also contains polarization and the bulk QP shift | Rename "QP gap − E_g(exp, bulk)" |
| 9 | Low | Exciton type classifier labels every state "Wannier" for a 2 nm dot (d_eh ≈ 10 Å ≈ R) | Use R-normalised criteria or drop label |
| 10 | Low | `dynamic_z` code docstring and CLI help claimed plasmon-pole f-sum rule | **Fixed** |
| 11 | Low | CLI help called `evgw`/`qsgw`/`sgw` GW variants | **Fixed** (honest labels) |
| 12 | Low | sBSE citation JCTC 18, 3054 in code and one doc page (correct: 3438) | **Fixed** |
| 13 | Low | `qsgw-*` clamps `eps_out` to ≥ 1 silently; elementwise `max(0,·)` clipping of ΔW persists | Open |
| 14 | Low | Diagonal Löwdin exchange uses `.real` of complex transition charges | Open (complex spatial MOs only) |

## 5. Documentation

Build: `sphinx-build -W` passes. All 4,027 rendered equations parse in MathJax
3 (checked by script). No raw TeX leaks into topic pages; the only
unrendered formulas are in API docstrings (autodoc, plain text).

Problems found and fixed in this pass:

* **34 bullet/numbered lists rendered as run-on paragraphs** (no blank line before the
  list). These are the "formatting issues" seen in the Auger, SOC, fuzzy-band,
  configuration and NTO pages. Fixed.
* **246 migration "From docs/partN/index.rst:a-b" headings and 73 stray `---`
  rules** were visible to readers. Removed; provenance stays in
  `audit/documentation_inventory.json`.
* **All section tables of contents were alphabetical** ("Foundations" was
  10th of 13 in Quasiparticles; "Overview" 12th in Recombination). They are now in
  teaching order: foundations → models → implementation → selection →
  configuration.
* Wrong callable `qdex.solver.solve` → `qdex.solver.ExcitonSolver.solve`. All
  32 documented callables now resolve by import.
* Resta "Penn" formula and MNOK mixing rule corrected to match the code, with notes.
* Overclaims corrected: the triplet-exchange explanation, "TDA … negligible loss", "Rashba … without
  empirical parameters", and the on-site/bulk Wannier–Mott limits in `interactions/architecture`.
* **14 new schematic/data figures** (`docs/_static/figures/*.svg`, generated
  by `docs/figures/make_figures.py`):
  - the workflow and the QP model map;
  - the anchor model curve, dielectric sphere and screening profiles;
  - the BSE kernel diagrams, excitation frameworks and qsGW loop;
  - SOC spinor construction, the NAMD workflow, Auger channels and TA;
  - the CdSe energy ladder and the environment-cancellation plot.

Still to do (editorial, not mechanical):

* Pages are split well by topic, but many still read as pasted chapter
  fragments. Examples: numbered headings "1. The Two-Particle Excitation Problem"
  inside a page called Foundations, and marketing phrases ("exceptionally
  robust", "without empirical parameters", "recovering the bulk Wannier–Mott
  binding"). Each page needs a short editorial pass with the template of
  the first audit (§III.2): assumptions → equations → implemented form →
  code contract → limits.
* The CdSe benchmark page quoted a local macOS path. It now points to `tests/CdSe` and carries the result tables.

## 6. Recommended priorities

1. One environment-aware W used in both QP and BSE (charge–dipole + continuum).
   This is the single change with the largest effect on predictive power.
2. Replace the "Penn" interpolation by a Penn-gap or Wang–Zunger form. Fit
   ε(R) against a small set of atomistic RPA calculations.
3. Make sBSE polarizability include dipoles (or default to AO mode) and
   report ε_eff sanity checks.
4. Convergence hygiene: thick-restart Davidson with a real path, active-space
   extrapolation for E_b, raise on `converged=False` unless `--allow-unconverged`.
5. A validation table of 3–4 CdSe sizes against experiment (Aubert–Hens 2022
   sizing) and against one full GW/BSE reference (e.g. MOLGW/BerkeleyGW on
   Cd₃₃Se₃₃).

## 7. Follow-up: environment-consistent W implemented (`qp_gap: env`)

Implemented recommendation §6.1 with the agreed choices:
* bulk ε∞ interior screening (Resta);
* the anchor residual off by default (opt-in `--env-anchor-residual`);
* a spherical cavity;
* `eps_out` = optical n².

`qdex/environment.py` builds the dielectric-sphere reaction field
W^refl_AB (multipole series; checked against the Born term, the Kelvin image
and the matched-medium limit in `tests/test_environment.py`). The same matrix
gives state-resolved QP shifts ±½ q_pᵀW^refl q_p, added to the bulk GW–PBE
opening, and is added to the Resta or xs-resta direct kernel.

Found and fixed along the way: with `soc_flag: true`, every model that passes
QP energies instead of a rigid scissor (`qsgw-*`, and now `env`) built the SOC
spinors from **DFT** energies. The SOC BSE therefore had no QP correction at all
(qsgw-dim: SOC S₁ would sit near the DFT gap). The spinor Hamiltonian now uses
the QP energies, and the SOC gap printout no longer double counts the scissor.

CdSe 2 nm, full comparison suite (`benchmarks/compare_models.py --profile full`),
32 runs, 0 FAIL:

| model | ε_out | QP gap | S₁ | binding |
|---|---|---|---|---|
| gw | 1.0 → 2.24 | 3.861 → 3.237 | 3.634 → 3.010 | 0.227 (fixed) |
| env | 1.0 → 2.24 → 6.2 | 3.852 → 3.115 → 2.727 | 2.501 → 2.503 → 2.500 | 1.351 → 0.612 → 0.227 |
| env + residual | 1.0 → 2.24 | 4.041 → 3.304 | 2.690 → 2.692 | 1.351 → 0.612 |

The optical gap is now solvent independent (dS₁/dQP = 0.00). Its value is set
by the interior physics: without the residual it is 0.2 eV below the
experimental window (2.70–2.95 eV, Yu et al. 2003 sizing curve); with the
monomer-calibrated residual it sits at the lower edge. The vacuum binding
(1.35 eV) is plausible in magnitude for an unscreened ~2 nm cluster, but it has not
been checked against a full GW–BSE reference.
**Correction to recommendation 2.** A probe with reduced interior
permittivity inside the `env` model shows that a size-dependent interior ε
cancels in S₁ in the same way as the surface term. Its contributions are:
* QP gap: ε_in = 4.7 → 3.047 eV and ε_in = 2.98 → 2.880 eV, plus an interior
  ΔW term of +0.1 and +0.26–0.33 eV respectively;
* binding: 0.64 and 0.70 eV, against 0.61 eV at ε_in = 6.2.

The net S₁ change is only a few tens of meV. Any static, classical change of W used
consistently in QP and BSE enters S₁ as ½(q_e − q_h)ᵀδW(q_e − q_h). The optical
gap is therefore set by DFT gap + bulk ΔGW − bulk-screened binding. Only
non-classical terms can move it: the size dependence of the short-range self-energy
(which the anchor residual approximates), the DFT gap itself, and the uncertainty
of the experimental size assignment. The decisive next test is a G₀W₀@PBE
calculation of this same cluster. The `sgw-*` models reach ≈ 3.0 eV in toluene only
because their interior term enters the QP gap without the matching binding.
