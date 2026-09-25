# QDEX scientific and architectural audit

> **Later verification:** The user subsequently supplied a CdSe 2.0 nm
> example folder, and the implementation and documentation were updated.
> See [post-fix verification](POST_FIX_VERIFICATION.md) and the
> [CdSe benchmark page](../docs/validation/cdse_benchmark.rst). The findings
> below describe the original audited snapshot and are retained as history.

Audit date: 25 September 2026. Baseline commit: `658c45ff461b48993a149c8d1e71d7e221e79894`, **including the pre-existing uncommitted working-tree changes**. File fingerprints are in `source_fingerprints.json`; line references below refer to that snapshot. Production code and existing documentation were not edited. This is an audit and reorganization proposal, not a certification of predictive accuracy.

## Assessment

QDEX contains useful reduced-interaction, static-screening, and exciton machinery, but several published-theory claims exceed what the implementation establishes. In its present form it should be described as a family of **bulk-calibrated quasiparticle correction models plus approximate static BSE/TDA**, with experimental orbital-relaxation extensions. The audited `sgw`, `evgw`, and `qsgw` paths are not interchangeable with the identically named conventional MBPT methods. `xs` is an analytical AO-density-pair integral approximation, not full four-index GW/BSE.

The most urgent issues are lost QP eigenvalues in the CLI; undefined/unstable frontier partitioning; mismatched full/diagonal BSE behavior; incomplete complex-orbital handling; and inconsistent claims about on-site screening, solvent cancellation, bulk convergence, and model provenance. The dimensional conversions inspected are largely consistent; the dominant errors are operator definitions, normalization, data flow, and unsupported physical extrapolation.

### Evidence and scope

* Direct inspection: `hardness.py`, QP dispatch/alignment in `cli.py`, Hamiltonian construction/actions, solver, Löwdin helpers, integral wrapper and C++ density-pair extraction. Documentation inventory covers all 13 RST files, 5,688 lines, in 298 disjoint migration units.
* Existing checks: **45 tests passed**: 16 in `test_physics_models.py`, 6 in `test_transition_modes.py`, 23 in `test_xs_kernel.py`. Logs are alongside this report. These are selected CPU suites, not the complete suite or a GPU validation. The installed compiled extension was used; it was not rebuilt from the dirty C++ source.
* Independent synthetic probes: `reproduce_findings.py`, `probe_results.json`, and `probes.log`. They demonstrate implementation behavior, not CdSe material accuracy.
* Literature: primary papers were checked through publisher records and accessible preprints. Some historical publisher pages expose only abstracts; no claim of equation-by-equation reproduction of paywalled papers is made.
* No 2 nm CdSe geometry, converged DFT orbitals, surface/passivation model, or reference GW/BSE data were supplied. Thus a physical CdSe spectrum, 55:45 split, or binding energy cannot be validated from this request alone. NAMD/Auger/analysis documentation is inventoried and mapped; their entire scientific implementations were not independently audited.

## I. Scientific soundness and theoretical foundations

### 1. What is established, and what is a QDEX model?

| Component | Established foundation | Actual QDEX status |
|---|---|---|
| Anchor interpolation | Bulk calibration and dielectric polarization motivate a constant plus finite-size corrections | Specific `11.52`, length regularization, residual exponent, anchor transfer and edge split are model choices |
| DIM | Self-consistent induced dipoles with short-range damping | Dipole response is subsequently normalized and mapped to a fitted pair-screening profile; this is not a computed microscopic inverse dielectric operator |
| Resta/Penn | Semiconductor dielectric models and gap-dependent electronic polarizability | QDEX's nearest-neighbor screening length and gap/radius interpolation need their own derivation/calibration |
| `sbse` | Simplified product-density Coulomb algebra and static response inversion | Related static screening approximation; additional atomic aggregation, solvent model, tables, and QP projection differ from the paper's complete method |
| `sgw` | Environmental self-energy differences can be modeled through changes in screening | Code projects only a clipped diagonal screening contrast and adds a bulk gap offset; it does not evaluate the paper's frequency-dependent self-energy |
| `dynamic_z` | Quasiparticle weight is a self-energy derivative | Scalar damping ansatz without an evaluated spectral response or enforced sum rule |
| `evgw-*` | Iterated eigenvalue corrections | Scalar gap/screening fixed point, not updates of a full set of eigenvalues in G and W |
| `qsgw-*` | Static Hermitian effective Hamiltonians can relax orbitals | Self-consistent atom-block ZDO-like ΔCOHSEX model plus bulk projector; not standard QSGW |

The Cho–Bintrim–Berkelbach method explicitly approximates orbital-pair products, then evaluates a dynamical self-energy with pole parameters; its name does not justify calling an arbitrary static ΔW projection sGW. Its paper is **JCTC 18, 3438–3446 (2022)**, DOI `10.1021/acs.jctc.2c00087`, not page 3054 as repeatedly given in QDEX. [Primary paper](https://arxiv.org/pdf/2109.04421), [published record](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.2c00087).

Thole supports damped interacting polarizabilities; it does not establish QDEX's `p_eff/max(p_eff)` mapping to W or semiconductor bulk calibration. Resta establishes semiconductor Thomas–Fermi screening with appropriate boundary conditions; merely using an exponential profile does not reproduce that construction. Bechstedt-related model dielectric functions similarly motivate model screening, but do not validate the specific QDEX formulas by association. [Thole](https://www.sciencedirect.com/science/article/pii/0301010481851762), [Resta](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.16.2717), [Cappellini et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.47.9892).

Delerue, Lannoo and Allan found substantial self-energy/electron-hole cancellation in the Si nanocrystals they studied. That is support for treating charged and neutral environmental corrections consistently, not proof of universal vxc cancellation or a CdSe parameterization. [Original paper abstract](https://pubmed.ncbi.nlm.nih.gov/11018909/).

### 2. Anchor interpolation

`estimate_gw_qp_gap` (`hardness.py:277`) implements

\[
\Delta(R)=\Delta_b+\frac{\kappa_{out}}{R+\ell}+A(R_0/R)^p,
\quad A=\Delta_0-\Delta_b-\frac{\kappa_{vac}}{R_0+\ell}.
\]

It reproduces the supplied **vacuum** anchor at R=R0 and tends to Δb for p>1. In a different environment, the anchor changes by `(κout−κvac)/(R0+ℓ)`, as intended. There is no uniqueness theorem: infinitely many interpolation functions satisfy two endpoints. Neither p=2 nor ℓ=1 Å is determined by those endpoints. A can have either sign; monotonicity is not guaranteed because

\[
\frac{d\Delta}{dR}=-\frac{\kappa_{out}}{(R+\ell)^2}
-pA\frac{R_0^p}{R^{p+1}}.
\]

For κout<0 or A<0, a turning point may be physical or may signal a poor fit. Require intermediate-size reference data and report the interpolation uncertainty. The code already clamps sub-anchor radii, or raises with `strict=True`; that avoids the residual divergence but does not make the clamped result physically predictive.

11.52 eV Å is about 0.8 times the Coulomb constant. It is not a unit conversion. Its physical averaging/calibration and whether any renormalization is already absorbed must be documented before adding another Z factor. The geometry routine adds 1.25 Å to the nuclear convex-hull equivalent-volume radius. Anchor radii and target radii must use the same definition; the database comment about relaxed convex hulls does not establish this.

### 3. ΔW, vxc cancellation, and the distinction between charging energy and COHSEX

With a fixed reference orbital and energy convention,

\[
\delta E_p^{QP}\simeq Z_p\langle p|\delta\Sigma-\delta v_{xc}|p\rangle
\]

is already an approximation. Subtracting a bulk reference additionally changes orbitals, occupied density matrix, screening, energy arguments and possibly Z. Local densities similar to bulk can make δvxc small in a sufficiently large, well-passivated interior. It does **not** cancel identically for polar surfaces, strain, reconstruction, ligand-induced dipoles, composition gradients, traps, different pseudopotentials/functionals, or sub-nanometer clusters. A bulk **gap** correction does not give separate absolute valence and conduction corrections.

At static COHSEX level, for one spin channel,

\[
\Sigma_{SEX}(r,r')=-\rho_{occ}(r,r')W(r,r';0),\qquad
\Sigma_{COH}(r,r')=\tfrac12\delta(r-r')[W(r,r;0)-v(r,r')|_{r'\to r}].
\]

At fixed occupied density, the difference is `−ρocc ΔW + ½δ ΔW(r,r)`. If the occupied density changes, even the SEX difference includes `−δρ Wbulk` in addition to `−ρQD ΔW`. A state projection of SEX involves **transition products with all occupied states**, not simply a state population multiplied by the total occupied populations. The population-only SEX equation in Part 3:308 is therefore also an extra approximation.

By contrast,

\[
U_p=\tfrac12\sum_{AB}q_A^p\Delta W_{AB}q_B^p
\]

is a classical polarization/charging functional of a normalized charge distribution. It can be useful, but is not the general COHSEX expectation value. In DIM/Resta the code treats its signed value as a removal/addition **magnitude**: positive UH lowers the occupied edge, positive UL raises the virtual edge. That convention must be separated from the signed eigenvalue shifts of the orbital-relaxation path. The `sgw` path instead uses `½Σ_A q_A ΔW_AA`, which is another operator altogether.

A useful exact limiting check within the discrete approximation is ΔW=c for every pair. Then `ΣSEX=−c Πocc`, `ΣCOH=c I/2`; occupied/virtual shifts are −c/2 and +c/2. The charging model reproduces those magnitudes for normalized q. Agreement in that limit does not establish agreement for spatially varying ΔW.

Elementwise clipping `maximum(0,WQD−Wbulk)` removes legitimate signed changes and is not projection onto a positive-semidefinite response. In DIM/Resta the contrast has zero diagonal but can have positive off-diagonal entries. A nonzero symmetric matrix with zero trace cannot be positive semidefinite; for two sites `[[0,d],[d,0]]` has eigenvalues ±d. Nonnegative orbital populations still give nonnegative charging energies for such a matrix, but arbitrary transition densities need not. Do not interpret that clipping as a passivity/stability guarantee.

### 4. DIM is a response model followed by a heuristic W mapping

`build_dim_screening_factors` (`hardness.py:645`) uses atomic polarizabilities in Bohr³, builds `M=α⁻¹+T`, inverts M, and applies three **uniform** fields. With code's tensor convention `T=(fT I−3 fD r̂r̂)/r³`, the sign is compatible with this M. The tensor written in Part 3:332 has the opposite undamped sign while retaining the minus-field convention; documentation must choose one convention consistently.

The code then takes the trace-averaged site response, sets `ηA=clip(pA/max(p),0.05,1)`, and defines

\[
\epsilon_{AB}=1+(\epsilon_b-1)\sqrt{\eta_A\eta_B},\quad
S_{AB}=\epsilon_{AB}^{-1}+(1-\epsilon_{AB}^{-1})e^{-k_{AB}R_{AB}}.
\]

Uniform fields determine a small set of response moments, not the site-to-site screened potential. For a linear polarizable model a more direct construction uses charge-generated local fields F and an induced-energy correction `W=v−F† M⁻¹ F`, with compatible damping and self terms. Three uniform-field solutions are not a substitute for those source fields.

Normalizing by the most polarizable site makes every site's screening depend on an extreme site, possibly a ligand. It discards absolute response scale and does not force interior ηA→1. Heteronuclear sites can retain different η values indefinitely. The asymptotic coefficient in the implemented pair formula is `1/εAB`, not necessarily `1/εbulk`. This prevents a general proof of the claimed DIM bulk limit. The separated-pair probe illustrates unequal η (1 and 0.645833); it is **not** a thermodynamic bulk simulation.

The raw table has Cd=48 and Se=31 Bohr³. Thus the request's narrative that Se is the more polarizable site is not even reflected in this particular table. Free-atom values do not establish ionic-in-crystal response or the sign of frontier asymmetry. Validate against cluster tensor polarizabilities, boundary fields, and a consistent bulk limit. Use solves with three RHS for the present uniform-field algorithm instead of an explicit inverse; diagnose stability rather than silently using a pseudoinverse after singularity.

### 5. Resta/Penn model and solvent consistency

The code's one-shot Resta dielectric is

\[
\epsilon_{eff}=1+\frac{\epsilon_b-1}{1+[(E_g^{DFT,QD}-E_g^{DFT,b})/\max(0.1,E_g^{DFT,b})]^2},
\]

when the QD gap exceeds the bulk PBE gap (`hardness.py:1484`). Iterated variants use differences from the bulk GW gap and normalize by that GW gap. This is **not** the documented `(Rp/R)²` expression. If ordinary confinement scales as R⁻², the squared gap-difference formula asymptotically changes ε by R⁻⁴, whereas the written radius formula changes it by R⁻². This is a functional discrepancy, not merely different notation.

`ks=√(εeff−1)/dNN` has correct inverse-length units, but is a nearest-neighbor-based interpolation, not a density-derived Thomas–Fermi wavevector. Inverted/small bulk PBE gaps make the artificial 0.1 eV floor particularly consequential. Define the intended model and calibrate it; do not swap formulas solely to match prose.

The solvent term `C(1/εout−1/εb)/√(rAB²+R²)` has the correct energy dimension and a plausible central-charge 1/R scale. It is a softened pair-distance ansatz. The exact dielectric-sphere reaction Green function depends on both radial positions relative to the boundary and their angle, not solely on their separation. It is insufficient for surface-localized traps, anisotropic particles, shells, or ligand-specific polarization.

Negative reaction contrast when εout>εb is physically possible: stronger external screening can reduce the charging gap. It must not be erased by a sign floor or converted to a huge edge partition. Also specify whether εout is electronic/optical or static; vertical electronic excitation and solvent-relaxed charging need different response assumptions.

Crucially, `build_resta_mnok` explicitly excludes εout from the **BSE direct** interaction; DIM does likewise. Their QP models include a reaction correction. Thus changing solvent shifts the QP gap without the corresponding direct electron-hole correction. For a simple polarization model,

\[
\delta E_{opt}=\tfrac12 q_h^T\delta Wq_h+\tfrac12 q_e^T\delta Wq_e-q_h^T\delta Wq_e
=\tfrac12(q_h-q_e)^T\delta W(q_h-q_e),
\]

before exchange, relaxation and dynamical complications. The code's split treatment does not reproduce this cancellation. Part 4:297 says dielectric mismatch enhances Kd, but its Resta/DIM implementation does not include that environmental term. If intentionally omitting environmental BSE polarization, state that scope explicitly and do not claim a unified QP/BSE environment model. The `sbse` branch has a different solvent subtraction using **2R**, so it does not cure cross-model inconsistency.

### 6. Approach B and absolute edge alignment

If UH,UL are positive magnitudes and their sum is nonzero, `fH=UH/(UH+UL)` is a defined phenomenological allocation. It can encode localization in a chosen population representation. It is not a unique MBPT replacement for separately calculated valence and conduction self-energies. `fL=1−fH` should hold, but current code separately divides both by `max(1e−12,UH+UL)`.

Applying this ratio to Δbulk is not automatically double counting: Δbulk can supply the bulk reference and ΔW the finite-size difference. It does assume, without evidence, that a *finite-size* removal/addition ratio also partitions the bulk correction. With fixed projectors, changing fH at fixed Δbulk adds a common energy offset: `Hbulk=Δbulk Πvirt−fH Δbulk I`. It does not independently break orbital symmetry or alter the gap when fractions sum to one. It does change absolute alignment, and becomes uncontrolled as both finite-size terms vanish. Symmetry problems arise when choosing one arbitrary member of a degenerate frontier manifold changes the scalar correction or when the approximate operators are not covariant.

Recommended design: keep separately calibrated signed bulk edge shifts, and add signed finite corrections independently:

\[
\delta\epsilon_H=\delta\epsilon_{v,b}-Z_H U_H,\qquad
\delta\epsilon_L=\delta\epsilon_{c,b}+Z_L U_L,
\quad\delta\epsilon_{c,b}-\delta\epsilon_{v,b}=\Delta_b.
\]

If only the bulk gap is known, declare the midpoint/occupied-fixed convention as a gauge; do not call it an absolute IP/EA prediction. The CLI's Approach B still reconstructs a PBE vacuum reference from monomer data whenever a 14-field material entry exists (`cli.py:1294`). “No monomer anchors” is therefore true at most for the gap calculation, not its reported absolute levels.

For degeneracy, use a subspace operator and diagonalize within the manifold, or a projector trace for a deliberately averaged population model. Averaging q before a nonlinear quadratic functional is not the same as averaging individual charging energies. State which invariant quantity is modeled; track subspace projectors/principal angles, not a single sorted column.

### 7. Dynamic Z and sum rules

The definition is `Zp=[1−∂ω ReΣp(Ep)]⁻¹`. Static COHSEX itself has no frequency derivative. A scalar static shift cannot determine it. For example, a causal one-pole model `Σ(ω)=g²/(ω−b)` gives `−Σ'(E)=g²/(E−b)²`: the shift and derivative depend differently on the pole separation. The proposed ratio can arise under restrictive one-pole assumptions with a chosen effective separation, but those assumptions are not derived from the available scalar ΔW energy.

Hybertsen–Louie constructs dielectric dynamics constrained by moments/sum rules. Godby–Needs fits response information at zero and an imaginary frequency. Neither is the scalar QDEX expression. [Hybertsen–Louie](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.34.5390), [Godby–Needs](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.1169), [primary study describing the two-frequency fit](https://infoscience.epfl.ch/server/api/core/bitstreams/5f57b8de-177d-446e-a5ad-6c32691f86de/content).

In QDEX (`hardness.py:1138`), `material_name` is unused, the default pole energy is **15 eV**, no `MATERIAL_PLASMON_ENERGIES` table exists, and the denominator is `max(1.01,εeff)−1`, not `max(1,εeff−1)`. The returned Z is clipped to **[0.50,0.99]**, not the documentation's [0.50,0.98]. A zero correction returns 0.99. Large-dot ΔW→0 would make the unclipped expression tend to **1**, not 0.80–0.85; a ΔW renormalizer must not be conflated with the full bulk quasiparticle residue.

For a pole energy called an energy, the plasma formula requires ℏ: `Ep=ℏ√(ne²/(ε0 me))` in SI, or `ℏ√(4πne²/me)` in Gaussian units. Atomic units suppress ℏ; density and conversion must then be explicit. A sum rule constrains an integral of the response, including electron count and spectrum. None is evaluated here. Rebrand this option as empirical/model dynamic damping until a causal frequency-dependent W, poles/residues, Σ and its derivative are implemented and checked.

### 8. On-site and bulk limits

“W→v=1/r” needs correction. The bare Coulomb singularity dominates sufficiently short real-space separation, often giving **W/v→1**; the regular correlation/reaction contribution need not vanish. An on-site **integral** samples finite distances and can be screened. Exact on-site v and W need not be equal. Unscreened exchange in BSE is a separate choice of operator; it does not require setting all on-site direct correlation to zero.

MNOK is finite at zero separation: `γAA=ηA` in eV for α=1. Resta/DIM impose SAA=1, so WAA=γAA rather than an actual 1/r divergence. The `sbse` matrix inverse can screen diagonal elements, and its `sgw` bulk reference uses `JAA/εb`; on-site conventions are not uniform across the hierarchy.

The long-distance limit for a finite object and the bulk limit do not commute. For finite R, the far field ultimately depends on the external medium. To recover a bulk exciton, take R large with the relative electron-hole coordinate finite on the exciton scale and far from the surface, then examine the macroscopic screened tail.

A static BSE can reduce to the Wannier equation under additional hypotheses: direct parabolic band extrema with a known reduced mass, smooth Bloch factors, converged transition space, long-range electronic screening, and controlled exchange/local-field corrections. These are not guaranteed by a `1/(εr)` tail alone. The basis being atom-centered does not force Frenkel character; correlation and orbital scales determine the exciton type.

**Diagonal BSE cannot generally recover a finite bulk Wannier binding energy.** In a volume-normalized extended pair state its diagonal Coulomb shift falls with system size. A localized relative exciton requires coherent mixing of many transitions. Only the coupled BSE path is a plausible route, with active-space and size convergence. The fixed active windows and the screening cutoff of 1,000 virtual states need independent convergence studies.

For Resta one-shot gap corrections, the internal contrast tends to zero only if the DFT gap and structural inputs converge to the same bulk reference. DIM's normalized heterogeneous response does not guarantee the required kernel limit. `sgw` uses a different on-site bulk reference and likewise supplies no proof. `qsgw` actually fails the zero-contrast bulk-offset check below. Smooth Wannier–Mott recovery and universal “90% binding from diagonal terms” should be presented as unvalidated claims, not consequences of the implementation.

### 9. AO orbital relaxation, Hermiticity and invariance

For a positive-definite overlap and complete orthonormal MO set,

\[
C^L=S^{1/2}C,\quad H^L=S^{-1/2}H^{AO}S^{-1/2}
=C^L\epsilon C^{L\dagger},\quad C^{QP}=S^{-1/2}U.
\]

These transformations in the real-valued code are structurally correct. `P=2 Cocc Cocc†` and `Πocc=P/2` make `−P⊙ΔW/2` consistent with a closed-shell one-spin SEX operator **within the assumed density-density integral approximation**. COH is `½ diag(diag ΔW)`; writing just `diag(ΔW)` ambiguously denotes a vector.

For real symmetric ΔW and Hermitian P, the Hadamard SEX is Hermitian. It is not the general four-index exchange contraction in a Löwdin basis. Orthogonalization does not erase differential overlap. Exact two-electron transformations involve four factors of S⁻¹/². Defining a ZDO model directly on labeled orthogonalized sites is possible, but is an additional model, not a rigorous transformation of atomic-site W.

Translation/rotation of a distance-based atom-block model is distinct from arbitrary basis covariance. The atom-block Hadamard model is covariant under unitary rotations **within an atom block** if labels and overlap transform consistently; it is not generally invariant under rotations mixing blocks. AO-density-only `xs` loses cross-orbital products required for general angular/basis covariance. “Hermitian” therefore does not establish “exact” or “basis invariant.”

Standard QSGW obtains a static Hermitian potential from a frequency-dependent Σ, conventionally symmetrizing matrix elements evaluated at the two state energies; it replaces the reference exchange-correlation contribution consistently. QDEX neither evaluates that Σ nor performs that construction. Its frozen HDFT plus bulk offset plus ZΔCOHSEX has no general starting-point-independence guarantee. [Kotani, van Schilfgaarde and Faleev](https://arxiv.org/pdf/cond-mat/0611002).

Additional numerical conditions: reject incomplete C for a claimed full-AO reconstruction, or solve in the supplied MO subspace; use overlap rank revelation instead of flooring negative/near-null eigenvalues; use conjugate transposes; validate closed-shell occupancy; and require density/Hamiltonian residual convergence in addition to gap change. Sorted-column sign flips do not track crossings or degenerate manifolds.

## II. Theory, documentation and code consistency

### 1. Units, basis and data flow

| Quantity/path | Checked result | Required qualification |
|---|---|---|
| `coords_ang / ANG_PER_BOHR`, `pdist` | Distances become Bohr | `ANG_PER_BOHR=0.52917721 Å/bohr`; input geometry must be Å |
| `eta_eV / HA_TO_EV`, MNOK inverse length, then `×HA_TO_EV` | Internal Hartree/Bohr dimensional path is consistent | Eta is an empirical energy, not a geometrical inverse length until divided by the Hartree scale |
| `ks_au=√(ε−1)/dNN_au`, `ks_ang=ks_au/ANG_PER_BOHR` | Exponent `ks_ang R_ang` is dimensionless | Choice of ks is a fitted profile, not a unit error |
| `14.3996 eV Å/√(R²+RQD²)` | Energy dimension is correct | This is the Coulomb constant, about 27.2114×0.529177; does not validate the Green function |
| `11.52 eV Å/(R+ℓ)` | Energy dimension is correct | `11.52/14.3996≈0.8000` is an effective prefactor/assumption |
| BSE matrices and `scissor_ev` | Eigenvalues and kernels are passed in eV | The chosen QP correction is a rigid scissor in current excitation construction |
| `xs` C++ output | Hartree density-pair integrals; Python converts to eV | The exactness claim must be limited to computed density-pair integrals |

**MNOK formula mismatch.** Both the request and Part 4:172 state `aAB=2/(ηA+ηB)`. Code (`hardness.py:525`, repeated in screening builders) uses `aAB=[1/ηA+1/ηB]/2` after expressing η in Hartree. They agree only if ηA=ηB. For Cd/Se values 3.4994/5.4795 eV, the two length parameters differ appreciably. The one-center limit `γAA=ηA` holds in both conventions, so on-site tests miss this. Choose a published convention, state its units and heteronuclear mixing rule, and test a heteronuclear pair.

**`xs` representation.** `compute_two_electron_ao` (`integrals.py:25`) calls C++ `compute_aabb_coulomb`, which stores only `(μμ|νν)` (`libint/integrals_core.cpp:404–437`). The subsequent Löwdin coefficient product creates approximate MO four-index integrals; it does not recover omitted `(μν|λσ)` elements. Since the raw Gaussian integrals are computed for original AOs and their density-pair matrix is later combined with **Löwdin** MO coefficient products, even that approximation is not the exact same object as the paper's `(μ'μ'|ν'ν')` in the orthogonalized basis. The user-visible description “exact analytical four-center integrals retaining full anisotropy” needs a precise qualifier: analytical **AO density-pair** integrals, with selected angular components; missing off-diagonal AO products can matter for p/d/f orbitals and basis rotations. If full ERIs are desired, expose the necessary four-index contractions or a proper density fitting/auxiliary basis pathway. If reduced cost is the goal, keep the present operator but label it accurately and benchmark against a full-integral reference.

**Several claimed APIs do not exist.** The docs name `estimate_anchor_scaled_qp_gap`, `_apply_qp_correction`, `build_damped_mnok_matrix`, `solve_diagonal_bse`, and `qdex.solver.davidson`; actual entry points are `estimate_gw_qp_gap`, inline CLI dispatch, `build_gamma`, `ExcitonHamiltonian.independent_transition_energies` through `ExcitonSolver.solve`, and `qdex.davidson.davidson`. `--ci_threshold` is documented but not a parser flag; the parser exposes `--e_thresh`/`--f_thresh`. Documentation also refers to a material plasmon-energy table that does not exist. The Part 3 `qsGW` page says 30 iterations; code defaults to 25.

### 2. Reproducible implementation findings and corrections

Priorities indicate risk of wrong scientific results, not implementation effort. Source references are current working-tree lines; the independent probes are in `probe_results.json`.

| Priority | Finding and evidence | Exact correction or gate |
|---|---|---|
| **P0** | **`qsgw-*` eigenvalues are discarded.** The estimators return `(scissor, provenance, C_qp, eps_qp)`, while `cli.py:1091`/`:1162` assigns only `C=C_qp`. Later the solver receives original `eps` plus a rigid `scissor_ev` (`cli.py:1762`, `:1805`); Fuzzy Bands likewise use a rigid correction. This combines relaxed wavefunctions with stale eigenvalues. | Carry a `QPState` containing `C_qp`, `eps_qp`, `eps_dft`, energy reference and flags to BSE, PDOS, Auger and SOC. For `qsgw`, pass `eps_qp` and `scissor_ev=0`, except in an explicitly documented extra correction stage. Preserve DFT energies separately for `independent_dft`. Test nonrigid QP level spacing end to end. |
| **P0** | **`qsgw` loses the bulk correction at zero ΔW.** With `penn_scaling=False`, `εout=εbulk`, 2-site orthonormal C, `qsgw-resta` returns 0.0 eV scissor instead of the CdSe database's 1.27 eV bulk gap offset. `tot_sig=max(1e−12,0)` makes `fH=fL=0` (`hardness.py:2004–2006`), so Hbulk is zero. | If finite contrast is below a scale-aware tolerance, use explicit calibrated bulk edge shifts; if unavailable, a declared 50:50 gap-only convention. Enforce `fH+fL=1` for finite meaningful allocations. Add ΔW→0 and εout=εbulk tests for both orbital models. |
| **P0** | **Signed solvent response destroys the edge split.** Synthetic `sgw-resta-pure` at εout=10>εbulk=6.2 returns scissor 0.9931 eV and `fH=fL≈−1.3844×10¹¹` because `tot_conf=max(1e−12, negative)` (`hardness.py:1587–1590`). The gap can be sensible while reported IP/EA and bulk projector become unbounded. | Never normalize signed sums with a positive floor. Separate a positive *reference allocation* from signed finite corrections; represent independent signed δEH/δEL. Validate dielectric positivity (`ε>0`) and reject unsupported εout values rather than silently clamp to 1. If retaining fractions, make them optional and undefined when denominator magnitude is too small or signs conflict. |
| **P0** | **Diagonal BSE omits exchange when direct attraction is disabled.** `ExcitonHamiltonian` computes diagonal Kx only inside `if self.include_exchange` (`:398–439`), although full `kernel_actions` computes Kx regardless of that flag. For a 3 AO, 2 atom synthetic basis with `charge_type=lowdin`, direct disabled, full-matrix diagonal `[5.9515,5.4656]` vs diagonal-mode energies `[3.87,5.27]`. | Compute Kx whenever `excitation_mode=diagonal_bse` and spin selection permits it; gate only Kd by `include_direct_eh`. Compare each diagonal-mode energy with `diag(Afull)` across singlet, triplet, UKS, spinor, both charge types and both integral representations. |
| **P1** | **Diagonal mode ignores requested Löwdin charge type.** The atom-centered diagonal path builds Mulliken populations/transition charges regardless of `charge_type` (`exciton_hamiltonian.py:134–165`, `:425–436`). With direct enabled and lowdin selected, full vs diagonal differs by 0.1205 eV in the same synthetic example. | Route diagonal densities and exchange through the same selected representation as full mode, or reject unsupported combinations. Ensure off-diagonal and diagonal code share one charge builder/convention. |
| **P1** | **Complex frontier population is phase dependent.** DIM/Resta and `sgw` compute `C**2`, not `abs(C)**2` (`hardness.py:1244`, `:1471`, `:1087`). A harmless global phase of π/4 changes a synthetic one-shot Resta scissor from 3.0403 to 2.1552 eV. The array assignment casts away the imaginary component. | For state charge `qA=Σμ∈A|Cμp^L|²`; use conjugation in pair/transition products. Propagate complex dtype through spinor paths. Add global-phase and degenerate-unitary invariance tests. |
| **P1** | **`qsgw` requires complete square MOs without checking.** `HDFT=C^L diag(eps) C^{L T}` yields a 3×3 matrix for a 3×2 supplied C; later multiplying with `P`/`Q` fails because `C_curr` is 3×2. The probe raises a broadcasting `ValueError`. | Require `C.shape==(nAO,nAO)` and `C†SC≈I` with all eigenvalues, or explicitly formulate an active-subspace Hamiltonian and frozen complement. Check overlap rank, MO energy count, occupancy and spin before allocating dense matrices. |
| **P1** | **`xs` claim overstates exactness and scope.** The C++ backend evaluates and stores `(μμ|νν)` only; `tests/test_xs_kernel.py` tests that restricted matrix. | Rename/docstring the operator `density_pair_coulomb`; add a small full-ERI comparison for physically relevant singlet/triplet matrix elements and orbital rotations before presenting it as a benchmark-quality exact-integral path. |
| **P1** | **`sbse` QP path is site-diagonal-only and inconsistent with its own screening reference.** `estimate_sgw_qp_gap` compares `diag(WQD)` with `diag(Jbare)/εbulk`, clamps differences and weighs populations linearly (`hardness.py:1070–1097`). DIM/Resta keep on-site interactions unscreened and use a distance-dependent Wbulk. This can create a positive on-site “finite-size correction” even when the actual bulk-reference diagonal convention should be the same. | Define one screened bulk reference built with the same operator and basis as WQD; compare full matrices/operators without elementwise sign clipping; derive the self-energy projection. Keep Cho's actual dynamic sGW as a separately named method only if implemented. |
| **P1** | **`dynamic_z` is mislabeled as sum-rule PPM.** The code has no frequency response or oscillator-strength moment; `material_name` is ignored, `ωp=15` is hardcoded default, and its denominator/clamps differ from docs. | Rename the current option `model_z`/“static-shift damping,” log Ω and clamps, and remove “f-sum-rule-derived” claims. A real PPM needs fitted response poles/residues and `∂Σ/∂ω` from the same model. |
| **P1** | **Reference alignment uses monomer data in Approach B.** `cli.py:1294–1309` reconstructs vacuum PBE edges from monomer HOMO/LUMO when available, then applies microscopic fractions. This conflicts with the documented anchor-free absolute IP/EA claim, and DFT orbital energies have an arbitrary zero unless a vacuum alignment exists. | Separate `gap_correction`, `bulk_edge_alignment`, and `vacuum_reference` outputs. Report absolute IP/EA only when an independent vacuum level or calibrated band-edge anchors are present. Label monomer-dependent output accordingly. |
| **P1** | **No general Wannier–Mott convergence test.** Existing tests check small matrices and Resta/DIM pair behavior, not size scaling of coupled BSE binding and exciton radius. | Benchmark a size series with mass, dielectric and transition-space convergence. Inspect `Ebind(R)`, relative wavefunction, optical gap, and environmental cancellation. Do not assert a guaranteed bulk limit for diagonal BSE. |
| **P2** | **Predictive scale is overstated.** For 10,000 atoms one dense 10,000² float64 pair matrix is 0.8 GB, not “<10 MB.” DIM stores/inverts a 30,000² dense matrix, at least 7.2 GB for one real matrix and far more in practice; dense inversion is cubic. QSGW's nAO² matrices become much larger. | Publish per-mode time/memory complexity and measured maximum system sizes. Replace explicit inverse with matrix solves or iterative response; support sparse/matrix-free methods. Do not advertise every model for 10,000 atoms. |
| **P2** | **Convergence and exceptional paths can silently mislead.** evGW/qsgw can return `converged=False`; qsGW checks only gap change, not density or eigenvectors. `max_iter=0` leaves outputs uninitialized. Singular DIM M falls back to a pseudoinverse. | Reject `max_iter<1`, nonfinite/unstable response and invalid dielectric data; expose convergence failure as an error unless user explicitly permits partial results. Log residuals and condition estimates. |
| **P2** | **Configuration docs and parser drift.** Missing function/flag names and differing defaults listed above; YAML `_apply_config` silently accepts any parser attribute under many section names but errors on unknown keys, so “complete configuration” needs executable schema checks. | Generate reference tables from parser/schema and API introspection, validate examples in CI, add links to authoritative topic pages. |

The priority list separates *fixes to the implemented model* from research choices. It does not propose changing numerical constants or physical model equations without reference benchmarks. A superficial sign flip in qsGW, clipping negative solvent terms, or inserting 55:45 fractions would hide rather than solve the underlying problem.

### 3. Boundary cases and how they should behave

* **Vacuum, εout=1:** anchor coefficient uses vacuum mismatch. DIM/Resta reaction term is positive for εbulk>1. `sbse` applies no solvent correction. These limits are code paths, not proof they describe the same physical boundary problem.
* **Matched media, εout=εbulk:** reaction term should vanish. Resta pure one-shot returns precisely the bulk database correction in the synthetic test, whereas qsgw-resta returns zero because of its fraction normalization. For finite-size internal dielectric contrast, matching only the exterior still permits a confinement term.
* **High dielectric, εout>εbulk:** negative reaction polarization is physically possible. Current fraction calculation becomes pathological. Model optical vs relaxed solvent response separately.
* **R→0/single site:** radius routine assigns 1 Å for a single site and 1.25 Å surface allowance for finite hulls; these are computational regularizations, not a controlled zero-radius limit. `sbse` uses a softened 2R solvent denominator. Degenerate geometry falls back to a max-distance radius. State an applicability domain in number of core shells and R/R0, and reject extrapolation by default for scientific runs.
* **Open shell/SOC:** the qsgw operator assumes P=2Πocc and real matrices. It does not accept unrestricted α/β occupations or complex spinor coefficients. `estimate_sgw_*` frontier magnitudes use one spin-free HOMO/LUMO. Do not silently apply the same correction to radicals or SOC-split degeneracies; build spin-resolved/spinor density matrices or restrict the model.
* **Degenerate frontiers:** individual HOMO/LUMO columns rotate arbitrarily. Scalar q, Z and f then change although the subspace is identical. Use subspace-invariant traces or diagonalize the correction operator in the degenerate manifold; compare subspace overlaps in iterative tracking.

### 4. CdSe 2.0 nm dimensional and numerical benchmark

The supplied core geometry is internally consistent: R=10 Å gives `4πR³/3=4188.79 Å³`. The remaining figures mix incompatible references:

| Item | Supplied request | Current `MATERIAL_DB['CDSE']` / implication |
|---|---|---|
| Bulk PBE gap | 1.14 or 1.25 eV | **0.64 eV** |
| Bulk GW gap | 2.76 eV | **1.91 eV** |
| Bulk opening | “1.51 eV” | 2.76−1.25=1.51; 2.76−1.14=**1.62**; repository 1.91−0.64=**1.27** eV |
| Dielectric | 6.20 | **6.20** |
| Exterior | 2.40 | a user parameter; consistent with the requested probe |
| Plasmon energy | 16 eV | no CdSe-specific value used by `compute_dynamic_z`; default **15 eV** |
| DFT gap + opening | 2.5–2.8 plus 1.5–1.9 eV | arithmetic gives **4.0–4.7 eV**, not the requested 3.3–3.7 eV QP range |
| Bulk exciton | 15 meV, 5.6 nm | at ε=6.2, 15 meV implies μ/m0≈0.0424 and a*=**7.74 nm**; a*=5.6 nm implies Ebind≈**20.7 meV** in the same hydrogenic model |

For a **model radius exactly 10 Å** and the current CdSe database, the anchor formula gives Δbulk=1.27 eV, `κout=2.94194 eV Å`, `A=0.580917 eV`, and total scissor ≈**1.69804 eV**. Combining that with an assumed DFT gap of 2.5–2.8 eV yields 4.20–4.50 eV. This is an arithmetic check of one QDEX formula, not a computed CdSe nanocrystal result; the effective radius from an actual atomic geometry need not equal 10 Å because of the 1.25 Å allowance. The request's 55:45 asymmetry and Z=0.82–0.88 are unvalidated target ranges. Their determination requires the actual frontier orbitals, reference gaps, choice of surface/passivation and a dynamic self-energy calculation. An exciton binding of 250–350 meV in toluene also requires a finite-size BSE/environment comparison; the present source tree alone cannot establish it.

A reproducible CdSe benchmark should specify phase, exact stoichiometric atom list and diameter convention, surface ligands/pseudopotentials, DFT functional and SOC status, MO/basis/energy alignment, dielectric frequency convention, GW reference method, active BSE space and convergence tolerances. Report both the repository database values and any alternative literature parameter set as separate named benchmarks rather than silently mixing them.

## III. Complete topic-based documentation reorganization proposal

The current suite has 13 RST files and 5,688 lines, including the 8 historical parts. The [inventory](documentation_inventory.md) routes **298 disjoint line/byte units** to proposed pages; [machine-readable hashes](documentation_inventory.json) give each original file and unit a SHA-256. This establishes a preservation plan, not a completed move. Existing pages remain intact until the equations and false API claims above are adjudicated. When migrating, preserve each old URL with a redirect or stable stub and check every anchor/cross-reference.

### 1. Proposed top-level navigation

```text
QDEX
├── Start
│   ├── Installation and capabilities
│   ├── Five-minute quickstart: spectra and dynamics
│   └── Configuration map and CLI/YAML precedence
├── Electronic structure
│   ├── Orbitals, AO overlap, orthonormality, population conventions
│   ├── PDOS, surface/core partitions, IPR, COOP
│   ├── Fuzzy bands, k paths, QP energy axes
│   └── Gaussian orbital cubes
├── Relativity and spin
│   ├── SOC origin, material effects, pseudopotentials
│   ├── GTH projectors and spinor Hamiltonian
│   └── Sparse assembly, UKS and spin-preserving excitations
├── Two-body interactions
│   ├── Operator/representation map: v, W, exchange, direct attraction
│   ├── MNOK and analytic AO density-pair integrals
│   ├── Resta/Penn, DIM/Thole, ZDO-RPA and simplified BSE kernels
│   ├── Environment, on-site behavior, length scales, unit conventions
│   └── Model applicability and cost
├── Quasiparticles
│   ├── GW fundamentals, DFT reference and bulk data
│   ├── Anchor-scaled gap model
│   ├── Static ΔW/charging models: DIM, Resta, sBSE path
│   ├── Absolute band edges and frontier allocation
│   ├── Model Z and energy-iteration fixed points
│   ├── AO orbital-relaxation model
│   └── Model comparison and limitations
├── Excitons and spectroscopy
│   ├── BSE/TDA, singlet/triplet/UKS/spinor sectors
│   ├── Four excitation frameworks
│   ├── Direct/exchange contractions and solvers
│   ├── Transition dipoles, oscillator strengths, spectra
│   └── Finite-size to bulk tests
├── Exciton analysis
│   ├── Transition densities and NTOs
│   ├── Plasser–Dreuw descriptors
│   └── Cubes and dashboards
├── Dynamics
│   ├── Trajectory data, cross-frame overlaps, phase/crossing tracking
│   ├── NACs, timestep separation, electronic propagation
│   ├── PME, CPA-FSSH-EDC, DISH
│   ├── Decoherence and spectral densities
│   └── Cooling kinetics, distributions and dashboards
├── Recombination and time-resolved spectroscopy
│   ├── Einstein radiative rates, gap law, FCWD, traps, PLQY
│   ├── Transient absorption and bleach kinetics
│   ├── Auger theory, screened matrix elements, spin/SOC
│   └── Static/dynamic Auger rates and ECSH
├── Workflows and validation
│   ├── Ground-state and SOC tutorials
│   ├── QP edge, four-framework and exciton-analysis tutorials
│   ├── Carrier-cooling/Auger tutorials
│   ├── CdSe 2 nm and perovskite benchmark cards
│   └── Convergence, invariance and known limitations
└── Reference
    ├── CLI/YAML schema and defaults
    ├── Material database provenance and units
    ├── API signatures and data structures
    ├── Citations and method lineage
    └── Glossary of symbols and energy conventions
```

The ordering is topic-based, not an eight-step series. One reader can follow only the experimental observables and tutorials; a computational chemist can follow operators and solver approximations; a developer can move from an equation to the precise implementation and tests. The interaction pages are shared: GW and BSE pages link to **one canonical definition** of each v/W representation, preventing contradictory on-site/solvent equations.

### 2. Template every topic page must follow

Every topic page should put a compact “theory ↔ code” panel next to the derivation rather than postponing implementation until the chapter end:

1. **Observable, assumptions, and validity domain:** e.g. vertical QP gap, optical root, NAMD population; state spin, reference energy, static/dynamic screening, size regime, environment and units.
2. **Derivation and matrix representation:** retain full mathematical steps and define every index, complex conjugate, occupancy factor, spin factor, and unit. Show exactly where a continuous operator becomes an atom or AO matrix.
3. **Implemented approximation:** present the **actual** equation from the current code, then a tagged distinction from literature methods or a proposed future equation. Display clipping, floors, active-space truncation and default parameters.
4. **Implementation contract:** module, exact callable signature, input shape/units, returned fields, CLI flag, YAML block, default, and path through dispatch. Generate signatures/defaults where practical from parser and code, using `api_signatures.json` as a starting inventory.
5. **Worked miniature example and limiting checks:** include one hand-reproducible matrix or geometry, symmetry/unit limits, configuration and expected output. Link to tests/benchmark data with method/version and reproducibility metadata.
6. **Known limitations and provenance:** distinguish a peer-reviewed result, a calibrated adaptation, and a QDEX heuristic. State references and dataset uncertainties near the equation they support.

For example, the ΔW page should put `U=½qᵀΔWq`, the actual `estimate_sgw_resta_qp_gap(..., penn_scaling=True, dynamic_z=False, ...)` signature, the `physics.qp_gap: sgw-resta` YAML block, and `--qp_gap sgw-resta` CLI invocation in the same subsection. It should state how the returned scissor is applied by `cli.py` and separately describe the experimental `qsgw-resta` orbitals/eigenvalues path. The BSE page should similarly pair `A=D+Kx−Kd` with `ExcitonHamiltonian.kernel_actions`, `ExcitonSolver.solve`, `--include-direct-eh`, `--excitation-mode`, and the corresponding YAML keys.

### 3. Source-to-topic preservation and implementation pairings

| Existing material to preserve | Proposed canonical topic and concrete implementation link |
|---|---|
| `getting_started/*`, `docs/index.rst` | Start + Reference/config; parser in `qdex.cli.main`, YAML mapping `_apply_config`; install/quickstart examples retained verbatim as historical examples until validated |
| Part 1 §§1–2 (KS basis, orthonormality, PDOS, population, surface/core) | Electronic structure/orbitals and populations; `qdex.io_utils`, `qdex.lowdin`, `qdex.pdos_coop.compute_pdos_and_coop`; `--charge_type`, `fuzzy.pdos_atoms` |
| Part 1 §§3–6 (IPR, COOP, Fuzzy, cubes) | Electronic structure/localization, bonding, bands and cubes; `qdex.orbital_analysis`, `qdex.pdos_coop`, `qdex.fuzzy_bands.run_fuzzy_bands_and_pdos`, `qdex.exciton_cube`; preserve equations, high-symmetry paths and output illustrations |
| Part 2 (Dirac reduction through UKS) | Relativity/foundations, pseudopotentials, spinor/UKS and assembly; `qdex.soc_utils.compute_spinor_subspace`, `compute_spinor_subspace_uks`; `--soc_flag`, `--gth_file`, `physics.soc_window_ev` |
| Part 3 §§1–5 (GW, shared engine, DFT gap, scaling) | Quasiparticles/foundations plus Two-body interactions/architecture; `qdex.hardness` and `qdex.solver.ExcitonSolver` routing; include method map and measured cost table |
| Part 3 §§6–7 (anchor, ΔW, DIM, Resta, sBSE, asymptotics, Approach B) | Quasiparticles/anchor, ΔW, alignment; Two-body/screening and environment; exact functions `estimate_gw_qp_gap`, `estimate_sgw_dim_qp_gap`, `estimate_sgw_resta_qp_gap`, `estimate_sgw_qp_gap`, `build_dim_screening_factors`, `build_sbse_kernel`; `physics.qp_gap`, `dynamic_z`, `eps_out` |
| Part 3 §8 (model Z, eigenvalue and orbital iteration) | Quasiparticles/model Z, gap iteration, orbital iteration; `compute_dynamic_z`, `estimate_evgw_dim_qp_gap`, `estimate_evgw_resta_qp_gap`, `estimate_qsgw_dim_qp_gap`, `estimate_qsgw_resta_qp_gap`; `--dynamic_z`, `--update_orbitals` |
| Part 3 §§9–15 (environment, edge alignment, database, matrices, presets, analysis, flags) | Two-body/environment; Quasiparticles/alignment/model selection; Reference/materials; Workflows/QP presets; Electronic structure/QP axes; `MATERIAL_DB`, `get_cluster_size_metrics`, `cli` provenance writer, `qdex.fuzzy_bands.build_qp_energies`; preserve **all** tables and per-material rows with a version/provenance column |
| Part 4 §§1–4 (TDA, spin, representation, kernels, bulk limits) | Excitons/foundations; Two-body/representations and screening; Validation/bulk limit; `qdex.integrals.compute_two_electron_ao`, `qdex.hardness.build_gamma`, `build_resta_mnok`, `build_xs_kernel`, `build_sbse_kernel`; `--2e-integrals`, `--kernel` |
| Part 4 §§5–10 (four frameworks, dipoles, Davidson, matrices, presets, flags) | Excitons/frameworks and solvers; Spectroscopy/dipoles; Workflows/exciton presets; `ExcitonHamiltonian.independent_transition_energies`, `kernel_actions`, `ExcitonSolver.solve`, `qdex.davidson.davidson`, `qdex.oscillator`; `--excitation-mode`, `--nroots`, `--tol` |
| Part 5 (transition density, seven Plasser–Dreuw metrics, NTOs, cubes, dashboard) | Exciton analysis/transition density, descriptors, NTOs, cubes, dashboard; `qdex.exciton_analysis.ExcitonAnalyzer`, `qdex.nto.compute_nto_pairs`, `qdex.exciton_cube`; `analysis.nto`, `analysis.plot`, cube YAML |
| Part 6 §§1–6 (pipeline, timestep, PME/FSSH/DISH, NAC, crossing, decoherence) | Dynamics/pipeline, propagation, PME/FSSH/DISH, NAC tracking, decoherence; `qdex.namd.precompute.align_phases_and_crossings`, `qdex.namd.integrator`, `master_equation.propagate_pme_tensor`, `surface_hopping.run_namd_dynamics`; `namd.engine`, `namd.dt_fs` |
| Part 6 §§7–11 (spectral density, recombination, PLQY, TA, analysis, config) | Spectroscopy/spectral density and transient absorption; Recombination/radiative/gap law/FCWD/traps/PLQY; Dynamics/analysis; `qdex.namd.analysis`, `qdex.namd.transient_absorption.compute_transient_absorption`, `qdex.hardness.compute_radiative_rates`, `compute_fcwd_rate`; retain all plots/diagrams and extraction recipes |
| Part 7 six tutorials | Workflows by scientific task, cross-linked from each topic; preserve complete input YAML, commands, outputs and observations; add an executable example validation job |
| Part 8 §§1–9 (Auger mechanisms, trions, rates, kernels, line shapes, SOC, ECSH, config) | Recombination/Auger foundations, matrix elements, screening, SOC, static/dynamic rates, ECSH and configuration; `qdex.auger.calculate_auger_rates`, `compute_auger_matrix_element`, `qdex.namd.surface_hopping`; `auger.*` YAML and CLI flags |
| `docs/api/index.rst` | Reference/API, supplemented by links from every topic and exact signatures for `qdex` modules, not a replacement for physical explanations |

The [298-row line inventory](documentation_inventory.md) gives the more granular destination for every heading, code block, table and inter-heading paragraph. Its routing includes repeated material from different parts; retain both copies during the initial move, record their provenance, then consolidate only after the scientific equations are reconciled and links point to one canonical derivation. Preserve the original page content as archival source through migration so no mathematical, tabular, or benchmark material disappears accidentally.

### 4. Editorial and technical migration gates

1. Freeze a source snapshot and compare the inventory hashes. Create topic pages by **moving complete disjoint units** in the manifest; ensure total migrated source bytes/units are accounted for before editorial rewrites. Keep original URLs redirected or stubbed. Verify Sphinx cross-references, formulas, tables, images and external reference links.
2. Put a “literature / QDEX adaptation / proposed” badge at each equation and a versioned citation at each material table. Reconcile the contradictory CdSe parameters, physical limits, model names and code signatures identified above before making new user-facing “recommended” pages.
3. Generate or check CLI/YAML defaults against parser and `_apply_config`. Execute every quickstart/preset example on a small supplied system. Keep a unit/shape glossary: AO versus Löwdin AO, atomic sites, Hartree/eV, bohr/Å, static/electronic/external permittivity, vacuum versus periodic energy zero.
4. Add a **validation atlas** linking tests and material benchmarks: heteronuclear MNOK, spin/SOC Hermiticity and phase invariance, full-versus-diagonal BSE, zero-contrast bulk offset, εout sweep, active-space/size convergence, solvent cancellation, model-Z behavior, memory scaling. State unavailable benchmarks rather than assigning expected ranges as successful results.
5. Set up three navigation entrances: “I need an absorption spectrum” for experimentalists; “Which approximation did this matrix element use?” for computational chemists; and “Which function and data structure implements this equation?” for developers. Each path lands on the same canonical physics and implementation page, with recipes and API links adapted to the reader.

### Recommended order of work after the audit

First correct the P0 data-flow and diagonal/operator bugs and add the limiting tests. Then make the model names and equations honest, establish a single environment/screening convention, and build a real CdSe reference data set. Only after these are stable should the lossless documentation move become the new main navigation. This order avoids giving a polished documentation structure to formulas that are still inconsistent with the code.
