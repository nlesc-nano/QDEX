Approximations for ΔW
=====================

Part of :doc:`/quasiparticles/index`. Background: :doc:`gw`.

Every QP model in QDEX has the form derived in :doc:`gw`:

.. math::

   \varepsilon_n^{\mathrm{QP}} = \begin{cases}
   \varepsilon_n - f_b\,\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n[\Delta W] + r_H\,s(R), & n \in \mathrm{occ}, \\
   \varepsilon_n + (1 - f_b)\,\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n[\Delta W] + r_L\,s(R), & n \in \mathrm{virt},
   \end{cases}

where:

* :math:`\varepsilon_n` is the KS (DFT) energy;
* :math:`\Delta_{\mathrm{bulk}}` is the tabulated bulk GW gap opening, split asymmetrically by :math:`f_b` (the anchor-derived bulk valence fraction: :math:`f_b = d_{h0} / (d_{h0} + d_{l0}) \approx 41.2\%` for CdSe, matching first-principles bulk GW literature);
* :math:`Z_n\,\Delta\Sigma_n[\Delta W]` is the finite-size self-energy correction of :math:`\Delta W`, weighted by the plasmon-pole renormalization factor :math:`Z_n`;
* :math:`r_{H/L}` (:math:`r_H` for occupied / valence states, :math:`r_L` for empty / conduction states) is the **anchor residual**, scaled by the confinement decay function :math:`s(R)` (:doc:`anchor`).

What is the anchor residual :math:`r_{H/L}`?
--------------------------------------------

At the anchor monomer :math:`R_0`, high-level benchmark calculations (evGW@PBE0 with complete basis sets) provide reference quasiparticle levels for the HOMO (:math:`\varepsilon_H^{\mathrm{ref}}`) and LUMO (:math:`\varepsilon_L^{\mathrm{ref}}`). This defines the benchmark band-edge shifts from PBE to QP:
:math:`d_{h0} = -(\varepsilon_H^{\mathrm{ref}} - \varepsilon_H^{\mathrm{PBE}})` (positive downward shift) and :math:`d_{l0} = \varepsilon_L^{\mathrm{ref}} - \varepsilon_L^{\mathrm{PBE}}` (positive upward shift).

The anchor residual :math:`r_{H/L}` is the signed difference between this benchmark shift and the shift predicted by the microscopic model at the anchor geometry :math:`R_0`:

.. math::

   r_H = d_{h0} - \left( f_b\,\Delta_{\mathrm{bulk}} + Z_H\,\Delta\Sigma_H[\Delta W](R_0) \right), \qquad
   r_L = d_{l0} - \left( (1 - f_b)\,\Delta_{\mathrm{bulk}} + Z_L\,\Delta\Sigma_L[\Delta W](R_0) \right).

The total gap residual is :math:`r_{\mathrm{gap}} = r_H + r_L`. Physically, :math:`r_{H/L}` captures all physical effects absent from a static dielectric model:

1. **Dynamic screening beyond the single-plasmon pole:** Dynamical vertex corrections and multiexcitonic/plasmon satellite screening not captured by static ΔCOHSEX.
2. **Starting-point and hybrid DFT effects:** The benchmark uses an evGW@PBE0 reference, while the nanocrystal KS orbitals are computed at PBE. The starting-point difference (band stretching from exact exchange) is naturally absorbed into :math:`r_{H/L}`.
3. **Microscopic chemical asymmetry & passivation:** Atomistic coordination, surface chlorine/ligand passivations, and localized atomic multipoles at the molecular boundary that continuous dielectric models smoothen out.

As dot radius :math:`R` increases toward the bulk crystal, the decay factor :math:`s(R) \in [0, 1]` smoothly turns off this molecular correction (:math:`s(R_0) = 1` and :math:`s(R \to \infty) = 0`), ensuring seamless convergence to the bulk GW limit.

The models differ in how much of ΔW they build, in order of increasing detail:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - ``qp_gap``
     - ΔW
     - ΔΣ
   * - ``brus``
     - none (ΔW = 0)
     - none; effective-mass confinement on the experimental bulk gap
   * - ``gw`` (two-anchor)
     - reaction field of a dielectric sphere with ε∞ inside
     - classical self-image, residual fixed at the anchor
   * - ``sgw-resta``
     - Resta profile with a size-dependent ε_in(R) + sphere reaction field
     - ΔCOHSEX for every orbital
   * - ``sgw-dim``
     - atom-resolved screening from polarizable dipoles + sphere reaction field
     - ΔCOHSEX for every orbital

``evgw-*`` and ``qsgw-*`` iterate the Resta and DIM models (section 7). The same W, including ΔW, is
the BSE kernel (:doc:`/excitons/screened_kernel`).

1. No anchor: ``brus``
----------------------

**ΔW = 0.** The dot is treated as bulk material with a kinetic confinement energy (Brus 1984,
without the polarization terms):

.. math::

   E_g^{\mathrm{QP}}(R) = E_g^{\mathrm{exp}}(\mathrm{bulk}) + \frac{\hbar^2\pi^2}{2\mu R^2},
   \qquad
   E_g^{\mathrm{QP}}(R) = \sqrt{E_g^2 + 2E_g\,\frac{\hbar^2\pi^2}{2\mu R^2}}\quad (E_g < 2\ \mathrm{eV}).

* :math:`E_g^{\mathrm{exp}}(\mathrm{bulk})`: the experimental bulk gap. It replaces the KS gap plus
  Δ_bulk.
* :math:`\hbar^2\pi^2/2\mu R^2`: the kinetic energy of an electron–hole pair of reduced mass μ in a
  sphere with infinite walls. The hyperbolic form corrects the non-parabolicity of narrow-gap
  materials.

The KS levels are shifted rigidly so that their gap equals :math:`E_g^{\mathrm{QP}}(R)`.

**Why.** It is the limit in which the dot has bulk screening everywhere. It needs no DFT gap, only
two bulk parameters.

**What it misses.**

* The surface polarization, so the QP gap is too small and does not depend on the solvent.
* The effective-mass kinetic term overestimates confinement in small dots.

With ΔW = 0 the consistent BSE kernel is the bulk W (``kernel: resta``). The missing surface
polarization then largely cancels in S₁ (:doc:`/excitons/cancellation`), so S₁ is more reliable than
the QP gap.

2. Two anchors: ``gw``
----------------------

**ΔW = the reaction field of a dielectric sphere.** The dot is a sphere of radius R with the bulk
ε∞ inside and ε_out outside. The potential at r of the polarization induced by a unit charge at r′
is

.. math::

   \Delta W(\mathbf r,\mathbf r') = G(\mathbf r,\mathbf r') = \frac{e^2}{R}\sum_{l\ge0}
   \frac{(\epsilon_\infty-\epsilon_{\mathrm{out}})(l+1)}{\epsilon_\infty\,[l\epsilon_\infty+(l+1)\epsilon_{\mathrm{out}}]}
   \Big(\frac{rr'}{R^2}\Big)^l P_l(\cos\theta).

* l = 0 is the Born term, constant inside the sphere.
* l ≥ 1 are the higher image multipoles. They grow toward the surface.

**ΔΣ in the classical limit** (:doc:`gw`, section 6). Each carrier sees its own image. Averaged over
the 1S envelope, the gap opens by

.. math::

   P(R) = F(\epsilon_\infty,\epsilon_{\mathrm{out}})\,\frac{e^2}{R},\qquad
   F = \Big\langle\sum_{l\ge0}\frac{(\epsilon_\infty-\epsilon_{\mathrm{out}})(l+1)}
   {\epsilon_\infty[l\epsilon_\infty+(l+1)\epsilon_{\mathrm{out}}]}\Big(\frac rR\Big)^{2l}\Big\rangle_{1S},

with F = 0.937 for CdSe in vacuum.

**The two anchors.** The QP correction interpolates between the bulk (Δ_bulk) and an evGW
calculation of the smallest cluster (radius R₀):

.. math::

   \Delta_{\mathrm{GW}}(R) = \Delta_{\mathrm{bulk}} + P(R;\epsilon_\infty,\epsilon_{\mathrm{out}})
   + A\,s(R),\qquad
   A = \Delta_{\mathrm{evGW}}(R_0) - \Delta_{\mathrm{bulk}} - P(R_0;\epsilon_\infty,1).

* A collects everything that is not the classical polarization of a sharp sphere with bulk ε∞: the
  reduced interior screening, the non-classical screened exchange, and the energy dependence of the
  bulk correction. For CdSe, A = −0.42 eV.
* :math:`s(R) = E_{\mathrm{conf}}(R)/E_{\mathrm{conf}}(R_0)` is the KS confinement energy relative to
  the anchor (:doc:`anchor`).
* Each band edge has its own curve, so the HOMO and LUMO reproduce the evGW frontier shifts at R₀.
* Every other orbital adds its own image term relative to the frontier orbital.

**Why.** Tight-binding GW for Si nanocrystals finds that the finite-size self-energy is dominated
by this surface polarization, with the *bulk* ε∞ inside (Delerue, Lannoo and Allan 2000, 2003).
The model is exact at R₀ by construction, exact in the bulk, and exact classical electrostatics at
large R. It costs nothing.

**What it misses.**

* **Everything non-classical is in one number, A.** Its size dependence, s(R), is assumed, not
  computed.
* **Only the radius enters.** Shape, facets, ligands and the actual orbitals enter only through R
  and the KS gap.
* **No screened exchange** beyond the classical limit.

3. Why an atomistic ΔW: Resta and DIM
-------------------------------------

The two-anchor model treats the dot as a uniform dielectric with bulk ε∞ and a sharp surface. Four
things are missing:

1. **Reduced interior screening.** A small dot has a larger gap, so it screens less:
   :math:`\epsilon_{\mathrm{in}}(R) < \epsilon_\infty` (Wang and Zunger, PRL 73, 1039 (1994)). Part of
   ΔW is therefore not at the surface but everywhere inside the dot.
2. **No screening at short range.** Electrons cannot screen a charge on the scale of a bond.
   :math:`W(r)\to v(r)` for r below the nearest-neighbour distance, and only at large r does
   :math:`W\to v/\epsilon`. A uniform 1/ε is wrong exactly where the orbitals overlap.
3. **The full ΔW(r, r′).** ΔSEX (:doc:`gw`, section 5) samples ΔW between points where orbital n and
   the occupied orbitals overlap, not only on the diagonal. The non-classical screened exchange
   needs ΔW as a matrix, with the right short-range behaviour.
4. **Geometry.** Facets, vertices, ligands, non-spherical shapes and core/shell structures screen
   differently from a sphere.

Resta covers 1–3 with one number per dot, ε_in(R). DIM covers 1–4 with one polarizability per atom.
Both have :math:`W_{\mathrm{QD}}\to W_{\mathrm{bulk}}` for large dots, so ΔW → 0 and the bulk limit is
exact. Both use the same sphere reaction field as ``gw`` for the environment (section 6). Both then
evaluate ΔΣ with the actual orbitals (section 4).

What follows is written on atom pairs A, B with the bare interaction :math:`\gamma_{AB}` (:doc:`representation`).

Resta: ``sgw-resta``
~~~~~~~~~~~~~~~~~~~~

.. figure:: /_static/figures/screening_profile.svg
   :width: 100%
   :alt: screening profile

   Resta screening for CdSe (ε∞ = 6.2, d_NN = 2.60 Å): unscreened at short range, 1/ε∞ at long range.

**Screening profile** (Resta, PRB 16, 2717 (1977)). Screening switches on over the Thomas–Fermi
length:

.. math::

   S_\epsilon(r) = \frac1\epsilon + \Big(1 - \frac1\epsilon\Big)e^{-k_s r},\qquad
   k_s = \frac{\sqrt{\epsilon-1}}{d_{NN}},\qquad
   W_{AB} = S_\epsilon(r_{AB})\,\gamma_{AB}.

* :math:`S_\epsilon \to 1` for r → 0: no screening on site and within a bond.
* :math:`S_\epsilon \to 1/\epsilon` for r → ∞: macroscopic screening.
* :math:`d_{NN}` is the nearest-neighbour distance, which sets the length over which the valence
  electrons screen.

**Interior dielectric constant** (Penn model, one oscillator). ε − 1 = (ħω_p/E_P)². Confinement
opens the average gap E_P by the KS confinement energy ΔE:

.. math::

   \epsilon_{\mathrm{in}}(R) = 1 + (\epsilon_\infty - 1)\Big[\frac{E_P}{E_P + \Delta E}\Big]^2,\qquad
   E_P = \frac{\hbar\omega_p}{\sqrt{\epsilon_\infty - 1}},\qquad
   \Delta E = E_g^{\mathrm{PBE}}(\mathrm{QD}) - E_g^{\mathrm{PBE}}(\mathrm{bulk}).

* ħω_p is the free-electron plasmon of the bulk valence (s, p) density: 14.1 eV for CdSe, so
  E_P = 6.2 eV.
* ε_in = 3.97 at 1.2 nm and 5.06 at 2 nm for CdSe.

**ΔW:**

.. math::

   W^{\mathrm{QD}}_{AB} = S_{\epsilon_{\mathrm{in}}}(r_{AB})\,\gamma_{AB},\qquad
   W^{\mathrm{bulk}}_{AB} = S_{\epsilon_\infty}(r_{AB})\,\gamma_{AB}.

**Why.** It is the simplest W that is right at both ends, unscreened at short range and bulk at long
range. One parameter, fixed by the dot's own KS gap, carries the size dependence.

**What it misses.** The dot is uniform inside, so there is no surface or ligand dependence of the
screening. ``sgw-resta-pure`` keeps ε_in = ε∞, leaving only the surface term.

DIM: ``sgw-dim``
~~~~~~~~~~~~~~~~

**Induced dipoles** (Applequist 1972; Thole 1981). Each atom, ligands included, carries a tabulated
polarizability α_A. In a uniform field E the induced dipoles solve

.. math::

   (\boldsymbol\alpha^{-1} + \mathbf T)\,\mathbf p = \mathbf E,

with Thole-damped dipole tensors T. An atom surrounded by other polarizable atoms responds more,
and an under-coordinated surface atom responds less.

**Local screening.** The response of each atom relative to the most polarizable one,
:math:`\eta_A = p_A/\max_B p_B \in [0.05, 1]`, sets a pair dielectric constant:

.. math::

   \epsilon_{AB} = 1 + (\epsilon_\infty - 1)\sqrt{\eta_A\eta_B},\qquad
   W^{\mathrm{QD}}_{AB} = S_{\epsilon_{AB}}(r_{AB})\,\gamma_{AB} .

The Resta profile is kept, so short range stays unscreened.

**Why.** The reduction of screening is placed where it occurs: at surface atoms, vertices and
ligands. It follows the actual geometry, so it is the model for non-spherical dots, core/shell
structures and ligand-rich surfaces.

**What it misses.**

* The polarizabilities are static and tabulated per element; they do not change with confinement
  (``evgw-dim`` rescales them with the gap, section 7).
* The mixing of η into ε is a model choice.

For spherical CdSe dots, Resta and DIM agree within 0.1 eV in the QP gap and S₁ at 2 nm after
calibration.

4. ΔΣ for all orbitals: one-shot ΔCOHSEX
----------------------------------------

For Resta and DIM, :math:`\Delta\Sigma_n` of :doc:`gw` (section 5) is evaluated for every orbital. In
the Löwdin basis :math:`c = S^{1/2}C`, with ΔW expanded to AO blocks
:math:`\Delta W_{\mu\nu} = \Delta W_{A(\mu)B(\nu)}` and :math:`P = 2\,c_{\mathrm{occ}}c_{\mathrm{occ}}^{\mathsf T}`:

.. math::

   \Delta\mathrm{COH}_n = \tfrac12\sum_\mu c_{\mu n}^2\,\Delta W_{\mu\mu},
   \qquad
   \Delta\mathrm{SEX}_n = -\tfrac12\sum_{\mu\nu} c_{\mu n}c_{\nu n}\,P_{\mu\nu}\,\Delta W_{\mu\nu}.

These are the continuous formulas with the integrals represented on the basis (:doc:`representation`).

* **Classical limit.** For a constant ΔW = c: ΔCOH = c/2 for every orbital, ΔSEX = −c (occupied)
  and 0 (empty), so the gap opens by c.
* **Beyond it.** The actual ΔW varies over the dot. At 1.2 nm (``sgw-resta``):

  - ΔSEX = −3.01 eV for the HOMO and −0.22 eV for the LUMO;
  - ΔCOH = +1.55 and +1.67 eV.

  The screened exchange raises the QP gap from 6.24 eV (classical) to 6.57 eV (evGW: 6.03 eV).
* **Classical option.** ``qp_selfenergy: classical`` replaces ΔΣ_n by ±½ q_nᵀΔW q_n, with q_n the
  atomic populations.
* **Cost.** One cached eigendecomposition of S and three matrix products for all orbitals.

5. Quasiparticle weight Z
-------------------------

.. math::

   Z_n = \Big[1 + \frac{|\Delta\Sigma_n|}{\tilde\omega}\Big]^{-1},\qquad
   \tilde\omega = \frac{\omega_p}{\sqrt{1 - 1/\epsilon_{\mathrm{eff}}}},

from one plasmon pole of the same dielectric model. ε_eff is ε_in (Resta) or the median pair
screening (DIM). Z ≈ 0.94–0.97 in vacuum and 0.98–0.99 in toluene (:doc:`dynamic_z`).

6. The environment term
-----------------------

The environment enters every model through the same sphere reaction field as ``gw``, with ε∞
inside:

.. math::

   W^{\mathrm{add}}_{AB} = G(\mathbf r_A,\mathbf r_B),\qquad
   \Delta W_{AB} = \max\big(0,\,W^{\mathrm{QD}}_{AB} - W^{\mathrm{bulk}}_{AB}\big) + W^{\mathrm{add}}_{AB}.

The interior part is clipped at zero: a dot never screens more than the bulk.
``qp_solvent_term: born`` replaces G by the older softened Born form
:math:`(1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty)\,e^2/\sqrt{r_{AB}^2 + R^2}`, which has no
position dependence of the self-image.

7. Iterated variants
--------------------

* **``evgw-resta`` / ``evgw-dim``.** ΔW is iterated against the QP gap, then section 4 is applied
  with the converged ΔW. With the current gap :math:`E_g^{(k)}`:

  - Resta recomputes ε_in with :math:`\Delta E = E_g^{(k)} - E_g^{\mathrm{GW}}(\mathrm{bulk})`;
  - DIM scales the polarizabilities by :math:`E_g^{\mathrm{GW}}(\mathrm{bulk})/E_g^{(k)}`.

  The loop is damped and stops when the gap changes by less than 10⁻⁴ eV. A larger QP gap screens
  less, so ΔW and the gap grow together.
* **``qsgw-resta`` / ``qsgw-dim``.** The same ΔW and COHSEX operators as a matrix in the Löwdin basis,
  diagonalized self-consistently:

  .. math::

     H^{\mathrm{eff}} = H^{\mathrm{KS}} + H^{\mathrm{bulk}} + Z\big(\Sigma^{\mathrm{SEX}} + \Sigma^{\mathrm{COH}}\big),\qquad
     \Sigma^{\mathrm{SEX}} = -\tfrac12 P\circ\Delta W,\quad \Sigma^{\mathrm{COH}} = \tfrac12\,\mathrm{diag}\,\Delta W .

  Orbitals and energies are updated until convergence. This adds orbital relaxation to ΔCOHSEX, at
  the cost of one diagonalization per iteration.

Both are self-consistent in ΔW only; the bulk part stays the tabulated Δ_bulk.

8. Anchor residual
------------------

Each Resta and DIM model is run once on the anchor cluster, Cd₁₆Se₁₃Cl₆ in vacuum, and compared with
evGW\@PBE0. The HOMO and LUMO errors r_H and r_L are stored and added to all occupied and all empty
orbitals, scaled by

.. math::

   s(R) = \frac{E_g^{\mathrm{PBE}}(\mathrm{QD}) - E_g^{\mathrm{PBE}}(\mathrm{bulk})}
               {E_g^{\mathrm{PBE}}(\mathrm{anchor}) - E_g^{\mathrm{PBE}}(\mathrm{bulk})} \in [0,1]

(1 at the anchor, 0.41 at 2 nm, 0.26 at 3.2 nm for CdSe). The gap residuals are −0.54 eV
(``sgw-resta``), −0.44 eV (``sgw-dim``) and −0.83 eV (``evgw-resta``), the same sign and size as the
two-anchor A. They contain what ΔCOHSEX misses and the PBE → PBE0 starting point of the reference
(:doc:`anchor`). ``qp_anchor_residual: off`` removes them.

9. Absolute levels (IP/EA)
--------------------------

The gap correction is split between HOMO and LUMO with the per-edge anchor curves, which are exact at
the anchor (41 %/59 % for CdSe) and tend to the bulk split for large dots. ``qp_edge_split: model``
uses the model's own split.

10. Summary of the approximations
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Approximation
     - Why
     - What it keeps / loses
   * - Bulk + ΔW split, Δ_bulk from a table
     - The bulk self-energy is transferable; only ΔW depends on the dot.
     - Exact bulk gap; the energy dependence of the bulk correction only through the residual.
   * - Static ΔCOHSEX
     - ΔW is long-range polarization with plasmon frequencies far above the level spacings.
     - Exact classical limit plus non-classical screened exchange; dynamics to first order via Z.
   * - Model dielectric (sphere, Resta, DIM) instead of an RPA χ₀
     - An RPA over the cluster costs as much as the GW it replaces.
     - Bulk screening at long range, size and surface dependence; no local-field detail beyond the
       atom.
   * - Penn ε_in(R) (Resta)
     - One oscillator whose gap opens with confinement.
     - Uniform interior reduction; no surface resolution.
   * - Thole dipoles (DIM)
     - Screening follows coordination and geometry.
     - Surface and ligand resolution; tabulated static polarizabilities.
   * - Plasmon-pole Z
     - Dynamics of ΔW to first order, without a frequency grid.
     - Z ≈ 0.94–0.99.
   * - Anchor residual
     - Catches what the model misses at one size.
     - Exact at the anchor; the E_conf scaling to other sizes is an assumption.

11. Which model
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Purpose
     - Recommendation
   * - Optical spectra, size series
     - ``sgw-resta`` (spherical dots) or ``sgw-dim`` (shape, ligands, core/shell), mnok, in the
       solvent. ``gw`` for the cheapest classical estimate. Add SOC for the band edge.
   * - QP gap, IP/EA
     - Any calibrated Resta/DIM model or ``gw``; all reproduce evGW at the anchor.
   * - Orbital relaxation
     - ``qsgw-*``; S₁ within about 0.1 eV of ΔCOHSEX, at much higher cost for large dots.
   * - Reproduce old results
     - ``qp_selfenergy: classical``, ``qp_anchor_residual: off``, ``qp_residual_scaling: power``,
       ``qp_solvent_term: born``, ``qp_polarization: legacy``.

**Other options,** kept for reference and not consistent in the sense of
:doc:`/excitons/cancellation`:

* ``pbe``: no correction.
* ``gw`` with ``qp_polarization: legacy``: the older κ/(R + ℓ) curve with the bulk kernel.
* ``sgw``: site-diagonal ΔW on the sBSE monopole RPA. It underscreens CdSe (ε_eff ≈ 1).
