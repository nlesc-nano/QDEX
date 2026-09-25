From GW–BSE to quantum dots
===========================

Part of :doc:`/quasiparticles/index`.

This page derives the quasiparticle (QP) and exciton layer of QDEX from many-body perturbation
theory. It shows which approximations are made, why each is needed for nanocrystals of 10³–10⁴
atoms, and what each one keeps or loses. The formulas of the individual models are on
:doc:`/quasiparticles/models`.

1. Quasiparticles in GW
-----------------------

A quasiparticle energy is the energy to add an electron to, or remove one from, the system. It
solves the Dyson equation, written here to first order around Kohn–Sham (KS) orbitals
:math:`\psi_n` with energies :math:`\varepsilon_n`:

.. math::

   \varepsilon_n^{\mathrm{QP}} = \varepsilon_n
   + Z_n\,\langle n|\,\Sigma(\varepsilon_n) - v_{xc}\,|n\rangle,
   \qquad
   Z_n = \Big[1 - \partial_\omega \langle n|\Sigma(\omega)|n\rangle\big|_{\varepsilon_n}\Big]^{-1}.

In the GW approximation the self-energy is the product of the one-particle Green's function and
the screened Coulomb interaction,

.. math::

   \Sigma(\mathbf r,\mathbf r';\omega) = \frac{i}{2\pi}\int d\omega'\,
   G(\mathbf r,\mathbf r';\omega+\omega')\,W(\mathbf r,\mathbf r';\omega')\,e^{i\eta\omega'},
   \qquad
   W = \epsilon^{-1} v,\quad \epsilon = 1 - v\chi_0 .

:math:`\chi_0` is the independent-particle response (RPA). Everything that distinguishes a
quasiparticle from a KS level is in W: how strongly the rest of the system screens an added
charge.

**Static COHSEX limit.** When the relevant frequencies of W (the plasmons, 10–20 eV) are large
compared with the level spacings of interest, W can be taken at ω = 0. Σ then splits into two
terms (Hedin 1965):

.. math::

   \Sigma^{\mathrm{SEX}}(\mathbf r,\mathbf r') = -\sum_{m\in\mathrm{occ}}
   \psi_m(\mathbf r)\psi_m^*(\mathbf r')\,W(\mathbf r,\mathbf r';0),
   \qquad
   \Sigma^{\mathrm{COH}}(\mathbf r,\mathbf r') = \tfrac12\,\delta(\mathbf r-\mathbf r')\,
   \big[W(\mathbf r,\mathbf r;0) - v(\mathbf r,\mathbf r)\big].

With diagonal elements

.. math::

   \mathrm{SEX}_n = -\sum_{m\in\mathrm{occ}} (nm|W|mn),
   \qquad
   \mathrm{COH}_n = \tfrac12\int d\mathbf r\,|\psi_n(\mathbf r)|^2\,[W - v](\mathbf r,\mathbf r).

**Screened exchange (SEX)** is the exchange with the occupied states, screened by W. It is
non-local and depends on how orbital n overlaps the occupied manifold. **Coulomb hole (COH)** is
the interaction of a charge with its own polarization cloud at the same point, which is local.
Static COHSEX overestimates gaps somewhat in bulk. Its value here is that it is exact for the
classical (image-charge) part of the finite-size correction; see section 4.

2. Excitons: the Bethe–Salpeter equation
----------------------------------------

The neutral excitations follow from the BSE. In the Tamm–Dancoff approximation (TDA), for
transitions i → a between occupied i, j and empty a, b states,

.. math::

   A_{ia,jb} = (\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}})\,\delta_{ij}\delta_{ab}
   + 2\,(ia|v|jb) - (ij|W|ab),

for singlets; triplets have no exchange term. The last term, the direct term, is the screened
electron–hole attraction. Written for a single pair, it is the mutual interaction of the electron
and hole densities, :math:`\int\!\!\int |\psi_a(\mathbf r)|^2\,W(\mathbf r,\mathbf r')\,|\psi_i(\mathbf r')|^2`.

The **same W** appears in Σ and in the BSE kernel. That is not a convention: the QP gap and the
exciton binding are two faces of one screening process. Section 4 shows what goes wrong when
they are computed with different W.

3. Why full GW–BSE is not used for quantum dots
-----------------------------------------------

A 3 nm CdSe dot with ligands has 10³ atoms and about 10⁴ basis functions. G₀W₀ scales as N⁴ in
common implementations, needs the dielectric matrix of the whole cluster, and many empty states.
A single calculation is possible on large machines, but not the routine calculation of size series,
ligand shells and solvents, and not molecular dynamics.

QDEX therefore separates what is **bulk-like** and transferable from what is **finite-size** and
must be computed for each dot. Following Delerue, Lannoo and Allan (PRL 84, 2457 (2000);
PRL 90, 076803 (2003)), write the screened interaction of the dot as

.. math::

   W_{\mathrm{QD}} = W_{\mathrm{bulk}} + \Delta W .

The self-energy splits the same way:

.. math::

   \Sigma_{\mathrm{QD}} - v_{xc} = \big[\Sigma_{\mathrm{bulk}} - v_{xc}\big] + \Delta\Sigma[\Delta W] .

* The **bulk part** is the GW correction of the bulk crystal. It is taken from a bulk GW
  calculation as the gap opening :math:`\Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW}} - E_g^{\mathrm{PBE}}`
  (1.27 eV for CdSe).
* The **finite-size part** :math:`\Delta\Sigma` is computed for each dot, in the static COHSEX form:

.. math::

   \Delta\Sigma_n = \underbrace{\tfrac12 \langle n|\Delta W(\mathbf r,\mathbf r)|n\rangle}_{\Delta\mathrm{COH}_n}
   \;\underbrace{-\sum_{m\in\mathrm{occ}}(nm|\Delta W|mn)}_{\Delta\mathrm{SEX}_n}.

The static form is well justified for ΔW. ΔW is dominated by the long-range surface polarization,
whose characteristic frequency is the plasmon, much higher than the band-edge level spacings. The
remaining dynamical effect is kept to first order through the quasiparticle weight Z
(:doc:`/quasiparticles/dynamic_z`).

**What ΔW contains:**

1. **Surface polarization (dielectric mismatch).** The dot has ε∞ ≈ 6 inside and a lower ε_out
   outside (1 in vacuum, 2.24 in toluene). An added charge polarizes the interface, and its image
   acts back on it.
2. **Reduced interior screening.** A small dot has fewer states to screen with: its
   polarizability per volume is lower than in bulk.
3. **The solvent,** which enters through ε_out.

4. The classical limit and why W must be shared
-----------------------------------------------

Take ΔW smooth on the scale of an atom, and let the occupied states act as a complete set for it,
:math:`\sum_{m\in\mathrm{occ}}\psi_m(\mathbf r)\psi_m^*(\mathbf r') \approx \delta(\mathbf r - \mathbf r')`
within the valence manifold. Then:

.. math::

   \Delta\mathrm{SEX}_n \to -\langle n|\Delta W(\mathbf r,\mathbf r)|n\rangle\ \ (n\ \mathrm{occupied}),
   \qquad \Delta\mathrm{SEX}_n \to 0\ \ (n\ \mathrm{empty}),

and

.. math::

   \Delta\Sigma_{\mathrm{empty}} = +\tfrac12\langle \Delta W(\mathbf r,\mathbf r)\rangle_e,\quad
   \Delta\Sigma_{\mathrm{occ}} = -\tfrac12\langle \Delta W(\mathbf r,\mathbf r)\rangle_h,\quad
   \Delta E_g^{\mathrm{QP}} = \Sigma^{\mathrm{pol}}_e + \Sigma^{\mathrm{pol}}_h .

This is the classical result. The gap opens by the **self-image** energies of the electron and the
hole (Brus, J. Chem. Phys. 80, 4403 (1984)).

In the BSE, the same ΔW changes the direct term by the **mutual image**, the interaction of the
electron with the image of the hole, :math:`\langle\Delta W(\mathbf r_e,\mathbf r_h)\rangle`.
For a dielectric sphere the image potential is a multipole series
(:doc:`/quasiparticles/anchor`):

* The **l = 0 (Born) term** is constant inside the dot. It is the same in the self- and the mutual
  image, so it cancels exactly between the QP gap and the binding.
* The **l ≥ 1 terms** survive in the self-image but average out in the mutual image.

The optical gap therefore keeps only a small classical size effect, about 0.1 e²/R. For CdSe:

.. list-table::
   :header-rows: 1

   * -
     - 1.2 nm
     - 2.0 nm
   * - self-image added to the QP gap
     - 2.54 eV
     - 1.46 eV
   * - mutual image added to the binding
     - 2.28 eV
     - 1.29 eV
   * - net effect on S₁
     - +0.26 eV
     - +0.17 eV

**Consequences.**

* **The QP gap is strongly size and solvent dependent; S₁ much less.** Tight-binding GW-BSE for Si
  nanocrystals shows the same near-cancellation (Delerue et al. 2000).
* **QP and BSE must use one W.** If the QP gap contains ΔW but the BSE kernel is the bulk W, the
  whole self-image ends up in S₁. For CdSe that is 1.3 eV at 2 nm in vacuum, and S₁ then shifts with
  the solvent almost as much as the QP gap. QDEX enforces the shared W; see
  :doc:`/quasiparticles/models`.
* **A better classical W alone cannot change S₁ much.** A size-dependent ε(R), atomistic dipoles or a
  smoother dielectric boundary change the self- and the mutual image almost equally.

5. What survives the cancellation
---------------------------------

The size dependence of S₁ beyond "KS gap + bulk correction" comes from the parts of ΔΣ that are
not classical:

* **Non-classical screened exchange.** The completeness assumption above fails in a small dot. The
  real :math:`\Delta\mathrm{SEX}_n = -\sum_{m\in\mathrm{occ}}(nm|\Delta W|mn)` depends on how each
  orbital overlaps the occupied states, and it has no counterpart in the direct BSE term. QDEX
  evaluates it for every orbital (one-shot ΔCOHSEX). At 1.2 nm it raises the ``sgw-resta`` QP gap
  from 6.24 to 6.57 eV (evGW: 6.03 eV) and S₁ by 0.43 eV, before the anchor calibration.
* **Energy dependence of the bulk correction (band stretching).** Confined states are built from
  bulk states away from the band edges, where the bulk GW correction differs from its band-edge value.
  To first order this adds :math:`s\,E_{\mathrm{conf}}^{\mathrm{KS}}(R)`, with E_conf the KS confinement
  energy. For CdSe the evGW anchor gives a negative net value (the gap opens less than
  classically), consistent with GW making the bands heavier near the edge.
* **Electron and hole distributions.** They are included through the actual orbitals.
* **Dynamical effects.** To first order they reduce the QP shift by Z and the binding by the same
  factor (Bechstedt et al., PRL 78, 1528 (1997)). QDEX applies one Z to both.

The non-classical remainder that the model does not capture is calibrated on an evGW\@PBE0
calculation of the smallest cluster (the anchor) and scaled to other sizes with
:math:`E_{\mathrm{conf}}^{\mathrm{KS}}(R)/E_{\mathrm{conf}}^{\mathrm{KS}}(R_0)`.
Because the anchor starts from PBE0 while the QDEX orbitals are PBE, the residual also absorbs the
part of the PBE → PBE0 starting-point change that eigenvalue self-consistency does not remove
(evGW keeps the PBE0 orbitals).

6. The approximations and why
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Approximation
     - Why
     - What it keeps / loses
   * - Bulk + ΔW split, bulk GW correction from a table
     - The bulk self-energy is transferable; only ΔW depends on the dot.
     - Keeps the bulk gap exactly; the energy dependence of the bulk correction enters only through
       the anchor residual.
   * - Static ΔCOHSEX
     - ΔW is dominated by long-range polarization with plasmon frequencies well above the band-edge
       level spacings.
     - Exact classical limit and non-classical screened exchange; dynamics to first order via Z.
   * - Atom-condensed interactions (MNOK or xs)
     - The interaction of transition densities is summarized by atom-pair (MNOK) or AO density-pair
       (xs) integrals; standard in sTDA-type BSE.
     - Keeps the long-range physics exactly and the short range approximately; the choice matters
       only for small clusters (1.2 nm: 0.6 eV in S₁; 2 nm: 0.05 eV).
   * - Model dielectric for W (Resta profile, DIM/Thole dipoles) instead of an RPA χ₀
     - An RPA over the cluster costs as much as the GW it replaces. The monopole RPA of sBSE
       underscreens badly (ε_eff ≈ 1 for CdSe).
     - Keeps bulk screening at long range and the size and surface dependence; loses local-field
       detail.
   * - Penn-gap ε_in(R) (Resta)
     - One oscillator whose frequency follows from the valence plasmon and ε∞; confinement opens
       its gap.
     - Moderate interior reduction (ε_in = 3.97 at 1.2 nm, 5.06 at 2 nm).
   * - Plasmon-pole Z from the same ε
     - The dynamics of ΔW to first order, without a frequency grid.
     - Z ≈ 0.94–0.97 (vacuum), 0.98–0.99 (toluene); the same Z scales ΔW in the kernel.
   * - Anchor calibration (evGW of the smallest cluster)
     - Catches what the model misses (band stretching, remaining non-classical terms) at one size.
     - Exact at the anchor; the E_conf scaling to other sizes is an assumption, removable with
       ``qp_anchor_residual: off``.
   * - TDA and static BSE kernel
     - Standard for band-edge excitons of II–VI dots; errors ≲ 0.05–0.1 eV.
     - No coupling to de-excitations; dynamical kernel effects only through Z.
