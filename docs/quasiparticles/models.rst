Approximations for ΔW
=====================

Part of :doc:`/quasiparticles/index`. Background: :doc:`gw`.

Every QP model in QDEX has the form derived in :doc:`gw`:

.. math::

   \varepsilon_n^{\mathrm{QP}} = \begin{cases}
   \varepsilon_n - f_b\,\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n[\Delta W], & n \in \mathrm{occ}, \\
   \varepsilon_n + (1 - f_b)\,\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n[\Delta W], & n \in \mathrm{virt},
   \end{cases}

where

* :math:`\varepsilon_n` is the KS (PBE) energy;
* :math:`\Delta_{\mathrm{bulk}}` is the bulk GW gap opening from ``MATERIAL_DB`` (bulk QSGW), and
  :math:`f_b` its valence share (41.2 % for CdSe);
* :math:`Z_n\,\Delta\Sigma_n[\Delta W]` is the finite-size self-energy of ΔW, weighted by the
  quasiparticle weight :math:`Z_n`.

No term is fitted to a cluster calculation: each approximation is used as it is.

The models differ in how much of ΔW they build, in order of increasing detail:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - ``quasiparticles.model``
     - ΔW
     - ΔΣ
   * - ``none``
     - none (ΔW = 0)
     - none; KS energies + bulk GW correction
   * - ``brus``
     - none (ΔW = 0)
     - none; effective-mass confinement on the experimental bulk gap
   * - ``sgw-resta``
     - Resta profile with a size-dependent ε_in(R) + sphere reaction field
     - ΔCOHSEX for every orbital
   * - ``sgw-dim``
     - atom-resolved screening from polarizable dipoles + sphere reaction field
     - ΔCOHSEX for every orbital

``evgw-*`` and ``qsgw-*`` iterate the Resta and DIM models (section 7). The same W, including ΔW, is
the BSE kernel (:doc:`/excitons/kernel`).

1. No finite-size correction: ``none`` and ``brus``
---------------------------------------------------

**ΔW = 0.** The dot is treated as bulk material. Both models are consistent with the bulk W in the
BSE (``kernel: resta`` or ``dim``). The missing surface polarization then largely cancels in S₁
(:doc:`/excitons/cancellation`), so S₁ is more reliable than the QP gap.

**``none``** keeps the KS energies and adds only the bulk correction:

.. math::

   \varepsilon_n^{\mathrm{QP}} = \varepsilon_n - f_b\,\Delta_{\mathrm{bulk}}\ (\mathrm{occ}),\qquad
   \varepsilon_n^{\mathrm{QP}} = \varepsilon_n + (1-f_b)\,\Delta_{\mathrm{bulk}}\ (\mathrm{virt}).

The confinement is the KS confinement of the actual orbitals. With ``excitations.mode: sbse`` and a
bulk kernel this is the sBSE: by the Delerue–Lannoo–Allan cancellation the surface polarization drops
out of the neutral excitation, and neither the QP gap nor the BSE needs it.

**``brus``** replaces the KS gap by the Brus kinetic confinement on the experimental bulk gap (without
the polarization terms):

.. math::

   E_g^{\mathrm{QP}}(R) = E_g^{\mathrm{exp}}(\mathrm{bulk}) + \frac{\hbar^2\pi^2}{2\mu R^2},
   \qquad
   E_g^{\mathrm{QP}}(R) = \sqrt{E_g^2 + 2E_g\,\frac{\hbar^2\pi^2}{2\mu R^2}}\quad (E_g < 2\ \mathrm{eV}),

with μ the reduced mass; the hyperbolic form corrects the non-parabolicity of narrow-gap materials.
The KS levels are shifted rigidly to this gap.

**What they miss.**

* The surface polarization, so the QP gap is too small and does not depend on the solvent.
* The non-classical screened exchange of a small dot, which does not cancel in S₁.
* ``brus`` only: the effective-mass kinetic term overestimates confinement in small dots.

2. The classical ΔW: reaction field of a dielectric sphere
----------------------------------------------------------

The dot is a sphere of radius R with the bulk ε∞ inside and ε_out outside. The potential at r of the
polarization induced by a unit charge at r′ is

.. math::

   W^{\mathrm{add}}(\mathbf r,\mathbf r') = G(\mathbf r,\mathbf r') = \frac{e^2}{R}\sum_{l\ge0}
   \frac{(\epsilon_\infty-\epsilon_{\mathrm{out}})(l+1)}{\epsilon_\infty\,[l\epsilon_\infty+(l+1)\epsilon_{\mathrm{out}}]}
   \Big(\frac{rr'}{R^2}\Big)^l P_l(\cos\theta).

* l = 0 is the Born term, constant inside the sphere.
* l ≥ 1 are the higher image multipoles. They grow toward the surface.
* Atom positions are capped at r/R = 0.9, where the series diverges at the boundary.

In the classical limit (:doc:`gw`, section 6) each carrier sees its own image. Averaged over the 1S
envelope, the gap opens by :math:`P(R) = F(\epsilon_\infty,\epsilon_{\mathrm{out}})\,e^2/R`, with
F = 0.937 for CdSe in vacuum (Brus 1984).

**Why.** Tight-binding GW for Si nanocrystals finds that the finite-size self-energy is dominated by
this surface polarization, with the *bulk* ε∞ inside (Delerue, Lannoo and Allan 2000, 2003). It is
exact classical electrostatics at large R. The Resta and DIM models use it as the environment part of
ΔW (section 6).

3. Why an atomistic ΔW: Resta and DIM
-------------------------------------

The sphere treats the dot as a uniform dielectric with bulk ε∞ and a sharp surface. Four things are
missing:

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
exact. Both add the sphere reaction field of section 2 for the environment and evaluate ΔΣ with the
actual orbitals (section 4).

What follows is written on atom pairs A, B with the bare interaction :math:`\gamma_{AB}` (:doc:`/integrals/representation`).

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

For spherical CdSe dots, Resta and DIM agree within about 0.1 eV in the QP gap and S₁ at 2 nm.

4. ΔΣ for all orbitals: one-shot ΔCOHSEX
----------------------------------------

For Resta and DIM, :math:`\Delta\Sigma_n` of :doc:`gw` (section 5) is evaluated for every orbital. In
the Löwdin basis :math:`c = S^{1/2}C`, with ΔW expanded to AO blocks
:math:`\Delta W_{\mu\nu} = \Delta W_{A(\mu)B(\nu)}` and :math:`P = 2\,c_{\mathrm{occ}}c_{\mathrm{occ}}^{\mathsf T}`:

.. math::

   \Delta\mathrm{COH}_n = \tfrac12\sum_\mu c_{\mu n}^2\,\Delta W_{\mu\mu},
   \qquad
   \Delta\mathrm{SEX}_n = -\tfrac12\sum_{\mu\nu} c_{\mu n}c_{\nu n}\,P_{\mu\nu}\,\Delta W_{\mu\nu}.

These are the continuous formulas with the integrals represented on the basis (:doc:`/integrals/representation`).

* **Classical limit.** For a constant ΔW = c: ΔCOH = c/2 for every orbital, ΔSEX = −c (occupied)
  and 0 (empty), so the gap opens by c.
* **Beyond it.** The actual ΔW varies over the dot. At 1.2 nm (``sgw-resta``):

  - ΔSEX = −3.01 eV for the HOMO and −0.22 eV for the LUMO;
  - ΔCOH = +1.55 and +1.67 eV.

  The screened exchange raises the QP gap from 6.24 eV (classical) to 6.57 eV (evGW: 6.03 eV).
* **Classical option.** ``quasiparticles.selfenergy: classical`` replaces ΔΣ_n by ±½ q_nᵀΔW q_n, with q_n the
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

The environment enters the Resta and DIM models through the sphere reaction field of section 2, with
ε∞ inside:

.. math::

   W^{\mathrm{add}}_{AB} = G(\mathbf r_A,\mathbf r_B),\qquad
   \Delta W_{AB} = \max\big(0,\,W^{\mathrm{QD}}_{AB} - W^{\mathrm{bulk}}_{AB}\big) + W^{\mathrm{add}}_{AB}.

The interior part is clipped at zero: a dot never screens more than the bulk.
``solvent_term: born`` replaces G by the older softened Born form
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

8. Summary of the approximations
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Approximation
     - Why
     - What it keeps / loses
   * - Bulk + ΔW split, Δ_bulk from a table
     - The bulk self-energy is transferable; only ΔW depends on the dot.
     - Exact bulk gap; no energy dependence of the bulk correction (band stretching).
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

9. Which model
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Purpose
     - Recommendation
   * - Optical spectra, size series
     - ``sgw-resta`` (spherical dots) or ``sgw-dim`` (shape, ligands, core/shell), mnok, in the
       solvent. Add SOC for the band edge.
   * - S₁ without any finite-size QP term
     - ``none`` with ``excitations.mode: sbse`` and a bulk kernel.
   * - QP gap, IP/EA
     - Any Resta/DIM model; ``qsgw-*`` for the most complete static treatment.
   * - Orbital relaxation
     - ``qsgw-*``; S₁ within about 0.1 eV of ΔCOHSEX, at much higher cost for large dots.
   * - Reproduce old results
     - ``selfenergy: classical``, ``solvent_term: born``.

**Other options,** kept for reference:

* ``pbe``: no correction.
* ``gw``: the earlier two-anchor interpolation (classical sphere + one fitted residual). Still used by
  the NAMD precompute; not recommended for new work.
* ``sgw``: site-diagonal ΔW on the sBSE monopole RPA. It underscreens CdSe (ε_eff ≈ 1).
