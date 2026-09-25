QP models with a shared W
=========================

Part of :doc:`/quasiparticles/index`. Background: :doc:`/quasiparticles/theory`.

Every QP model in QDEX that defines a screened interaction W passes that W to the BSE, so the QP
correction and the electron–hole attraction come from one W (theory, section 4). This page gives
the formulas of the models as implemented.

.. list-table::
   :header-rows: 1
   :widths: 22 40 38

   * - ``qp_gap``
     - W of the QP correction
     - BSE kernel
   * - ``gw`` (two-anchor)
     - bulk Resta W + reaction field of a dielectric sphere
     - ``resta-sphere``: the same W
   * - ``sgw-resta``, ``evgw-resta``, ``qsgw-resta``
     - Resta W with the Penn ε_in(R) + sphere reaction field
     - ``qp``: the same W
   * - ``sgw-dim``, ``evgw-dim``, ``qsgw-dim``
     - DIM/Thole W + sphere reaction field
     - ``qp``: the same W

Combining a model with a different kernel is an error; ``--allow-inconsistent-kernel`` exists only
to reproduce old results. ``two_electron_integrals`` (``mnok`` or ``xs``) selects the representation
of that one W for both the QP correction and the kernel.

1. Common ingredients
---------------------

**Bare interaction (MNOK).** For atoms A and B with hardness η:

.. math::

   \gamma_{AB} = \frac{1}{\sqrt{r_{AB}^2 + a_{AB}^2}},\qquad a_{AB} = \tfrac12\big(\eta_A^{-1} + \eta_B^{-1}\big).

**Resta screening profile** with dielectric constant ε, and the bulk reference:

.. math::

   S_\epsilon(r) = \frac1\epsilon + \Big(1 - \frac1\epsilon\Big)e^{-k_s r},\qquad k_s = \frac{\sqrt{\epsilon-1}}{d_{NN}},
   \qquad W^{\mathrm{bulk}}_{AB} = S_{\epsilon_\infty}(r_{AB})\,\gamma_{AB}.

**Solvent / surface term.** The reaction field of a dielectric sphere with ε∞ inside and ε_out
outside, the same Green function as the two-anchor model (section 7):

.. math::

   W^{\mathrm{add}}_{AB} = G(\mathbf r_A,\mathbf r_B) = \frac{e^2}{R}\sum_{l\ge0}
   \frac{(\epsilon_\infty-\epsilon_{\mathrm{out}})(l+1)}{\epsilon_\infty\,[l\epsilon_\infty+(l+1)\epsilon_{\mathrm{out}}]}
   \Big(\frac{r_A r_B}{R^2}\Big)^l P_l(\cos\theta_{AB}).

Its diagonal grows toward the surface through the l ≥ 1 multipoles. The earlier softened Born form,
:math:`(1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty)\,e^2/\sqrt{r_{AB}^2 + R^2}`, has no position
dependence of the self-image; it is available as ``qp_solvent_term: born``.

**The finite-size part of W:**

.. math::

   \Delta W_{AB} = \max\big(0,\; W^{\mathrm{QD}}_{AB} - W^{\mathrm{bulk}}_{AB}\big) + W^{\mathrm{add}}_{AB} .

The models differ only in :math:`W^{\mathrm{QD}}`.

2. ``sgw-resta``: Resta screening with a Penn interior dielectric constant
--------------------------------------------------------------------------

The interior of a dot screens less than the bulk because confinement opens the gaps of the states
that do the screening. In the Penn model (one oscillator), ε − 1 = (ħω_p/E_P)², and confinement
shifts the average gap E_P by ΔE:

.. math::

   \epsilon_{\mathrm{in}}(R) = 1 + (\epsilon_\infty - 1)\Big[\frac{E_P}{E_P + \Delta E}\Big]^2,\qquad
   E_P = \frac{\hbar\omega_p}{\sqrt{\epsilon_\infty - 1}},\qquad
   \Delta E = E_g^{\mathrm{PBE}}(\mathrm{QD}) - E_g^{\mathrm{PBE}}(\mathrm{bulk}),

.. math::

   W^{\mathrm{QD}}_{AB} = S_{\epsilon_{\mathrm{in}}}(r_{AB})\,\gamma_{AB}.

ħω_p is the free-electron plasmon of the bulk valence (s, p) density: 14.1 eV for CdSe, so
E_P = 6.2 eV. For CdSe, ε_in = 3.97 at 1.2 nm and 5.06 at 2 nm.

*Why:* the Resta profile keeps the correct bulk limit at long range and becomes unscreened at short
range, and one parameter (ε_in) carries the size dependence. It is the cheapest W that is right at
both ends. ``sgw-resta-pure`` keeps ε_in = ε∞ (surface term only).

3. ``sgw-dim``: atomistic dipole screening
------------------------------------------

Each atom carries a polarizability α_A. The induced dipoles in a uniform field solve

.. math::

   (\boldsymbol\alpha^{-1} + \mathbf T)\,\mathbf p = \mathbf E,

with Thole-damped dipole tensors T (Applequist 1972; Thole 1981). The response of each atom
relative to the most polarizable one, :math:`\eta_A = p_A/\max_B p_B \in [0.05, 1]`, measures its
local screening. Surface atoms with fewer neighbours respond less. The pair dielectric constant and
W are:

.. math::

   \epsilon_{AB} = 1 + (\epsilon_\infty - 1)\sqrt{\eta_A\eta_B},\qquad
   W^{\mathrm{QD}}_{AB} = S_{\epsilon_{AB}}(r_{AB})\,\gamma_{AB} .

*Why:* the reduction of screening is placed where it physically occurs, at under-coordinated
surface atoms, and follows the actual geometry and ligand shell. For CdSe the interior contrast is
smaller than with the Penn scaling. After calibration, Resta and DIM agree within 0.1 eV in the QP gap
and S₁ at 2 nm.

4. QP energies of all orbitals: one-shot ΔCOHSEX
------------------------------------------------

The finite-size self-energy of every orbital n (theory, section 3) is evaluated in the Löwdin basis
c = S^½ C, with ΔW expanded to AO blocks :math:`\Delta W_{\mu\nu} = \Delta W_{A(\mu)B(\nu)}` and the
density matrix :math:`P = 2\,c_{\mathrm{occ}} c_{\mathrm{occ}}^{\mathsf T}`:

.. math::

   \Delta\mathrm{COH}_n = \tfrac12\sum_\mu c_{\mu n}^2\,\Delta W_{\mu\mu},
   \qquad
   \Delta\mathrm{SEX}_n = -\tfrac12\sum_{\mu\nu} c_{\mu n}c_{\nu n}\,P_{\mu\nu}\,\Delta W_{\mu\nu},
   \qquad
   \Delta\Sigma_n = \Delta\mathrm{COH}_n + \Delta\mathrm{SEX}_n .

Quasiparticle weight from one plasmon pole of the same dielectric model
(:doc:`/quasiparticles/dynamic_z`):

.. math::

   Z_n = \Big[1 + \frac{|\Delta\Sigma_n|}{\tilde\omega}\Big]^{-1},\qquad
   \tilde\omega = \frac{\omega_p}{\sqrt{1 - 1/\epsilon_{\mathrm{eff}}}} .

Here ε_eff is ε_in (Resta) or the median inter-atom screening (DIM).

QP energies, with :math:`\Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW}}(\mathrm{bulk}) - E_g^{\mathrm{PBE}}(\mathrm{bulk})`:

.. math::

   \varepsilon_n^{\mathrm{QP}} = \varepsilon_n - \tfrac12\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n + r_H\,s(R)\quad (n\ \mathrm{occupied}),

.. math::

   \varepsilon_n^{\mathrm{QP}} = \varepsilon_n + \tfrac12\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n + r_L\,s(R)\quad (n\ \mathrm{virtual}).

**Anchor residual.** r_H and r_L are the model's HOMO and LUMO errors against the evGW anchor,
calibrated once per material and model (:doc:`/quasiparticles/anchor`). The size scaling is:

.. math::

   s(R) = \frac{E_g^{\mathrm{PBE}}(\mathrm{QD}) - E_g^{\mathrm{PBE}}(\mathrm{bulk})}
               {E_g^{\mathrm{PBE}}(\mathrm{anchor}) - E_g^{\mathrm{PBE}}(\mathrm{bulk})} \in [0, 1].

s(R) is 1 at the anchor, 0.41 at 2 nm and 0.26 at 3.2 nm for CdSe.

For a constant ΔW = c this gives ΔCOH = c/2 for every orbital, ΔSEX = −c (occupied) and 0 (empty),
so the gap opens by c: the classical limit. The actual ΔW varies over the dot, and the screened
exchange then depends on each orbital's overlap with the occupied states. At 1.2 nm
(``sgw-resta``) ΔSEX is −3.01 eV for the HOMO and −0.22 eV for the LUMO, and ΔCOH is +1.55 and
+1.67 eV. This part does not cancel against the binding.

* **Classical option.** ``qp_selfenergy: classical`` replaces ΔΣ_n by the classical charging term
  ±½ q_nᵀ ΔW q_n, with q_n the atomic populations.
* **Cost.** One eigendecomposition of S (cached) and three matrix products for all orbitals.

5. The BSE kernel
-----------------

.. math::

   W^{\mathrm{BSE}} = W^{\mathrm{bulk}} + \bar Z\,\big(W^{\mathrm{QD}} + W^{\mathrm{add}} - W^{\mathrm{bulk}}\big),
   \qquad \bar Z = \tfrac12 (Z_{\mathrm{HOMO}} + Z_{\mathrm{LUMO}}),

with the same Z as the QP levels. The same pole that reduces the QP shift also makes the
electron–hole interaction dynamical; to first order the two cancel in the neutral excitation
(Bechstedt et al. 1997). The TDA Hamiltonian is

.. math::

   A_{ia,jb} = (\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}})\delta_{ij}\delta_{ab}
   + 2\sum_{AB} q^{ia}_A\,\gamma_{AB}\,q^{jb}_B - \sum_{AB} q^{ij}_A\,W^{\mathrm{BSE}}_{AB}\,q^{ab}_B,

with transition charges q (Mulliken by default).

6. Iterated variants
--------------------

* **``evgw-resta`` / ``evgw-dim``.** W is iterated against the gap, then section 4 is applied with
  the converged W. With the current QP gap :math:`E_g^{(k)}`:
  - Resta recomputes ε_in with :math:`\Delta E = E_g^{(k)} - E_g^{\mathrm{GW}}(\mathrm{bulk})`;
  - DIM scales the atomic polarizabilities by :math:`E_g^{\mathrm{GW}}(\mathrm{bulk})/E_g^{(k)}`.

  The loop is damped and stops when the gap changes by less than 10⁻⁴ eV.
* **``qsgw-resta`` / ``qsgw-dim``.** The same ΔW and the same COHSEX operators, but as a matrix in the
  Löwdin AO basis, diagonalized self-consistently:

  .. math::

     H^{\mathrm{eff}} = H^{\mathrm{DFT}} + H^{\mathrm{bulk}} + Z\big(\Sigma^{\mathrm{SEX}} + \Sigma^{\mathrm{COH}}\big),\qquad
     \Sigma^{\mathrm{SEX}} = -\tfrac12 P\circ\Delta W,\quad \Sigma^{\mathrm{COH}} = \tfrac12\,\mathrm{diag}\,\Delta W .

  Orbitals and energies are updated until convergence. This adds orbital relaxation to the one-shot
  ΔCOHSEX, at the cost of one full diagonalization per iteration.

7. The two-anchor ``gw`` model
------------------------------

``gw`` interpolates the PBE → QP correction between bulk GW and the evGW anchor. The finite-size
term is the self-image energy of a dielectric sphere,
:math:`P(R) = F(\epsilon_\infty,\epsilon_{\mathrm{out}})\,e^2/R`. The non-classical residual A is
fixed at the anchor and scaled by s(R). Each band edge has its own curve
(:doc:`/quasiparticles/anchor`).

The kernel is the bulk Resta W plus the reaction field of the same sphere,

.. math::

   G(\mathbf r,\mathbf r') = \frac{e^2}{R}\sum_{l\ge0}
   \frac{(\epsilon_\infty-\epsilon_{\mathrm{out}})(l+1)}{\epsilon_\infty\,[l\epsilon_\infty+(l+1)\epsilon_{\mathrm{out}}]}
   \Big(\frac{rr'}{R^2}\Big)^l P_l(\cos\theta).

The HOMO and LUMO shifts come from the anchor curves. Every other orbital adds its own image term
relative to the frontier orbital. ``gw`` has no screened-exchange term: it is exact at the anchor by
construction and classical elsewhere.

8. The xs representation
------------------------

With ``two_electron_integrals: xs`` the same W is represented on exact AO density-pair integrals:

.. math::

   W_{\mu\nu} = \frac{W_{A(\mu)B(\nu)}}{\gamma_{A(\mu)B(\nu)}}\,(\mu\mu|\nu\nu) + W^{\mathrm{add}}_{A(\mu)B(\nu)}.

ΔCOHSEX and the kernel then use AO populations. The two representations agree at long range. At
short range xs gives more binding: −0.6 eV in S₁ at 1.2 nm, −0.05 eV at 2 nm.

9. Absolute levels (IP/EA)
--------------------------

The model's gap correction is split between HOMO and LUMO with the per-edge anchor curves, which
are exact at the anchor (41/59 % for CdSe) and tend to the bulk split for large dots. The split is
applied to vacuum-referenced PBE levels. ``qp_edge_split: model`` uses the model's own split.

10. Other options
-----------------

* ``pbe``: no correction.
* ``brus``: effective-mass kinetic confinement on the experimental gap. It defines no W and uses an
  independent kernel.
* ``gw`` with ``qp_polarization: legacy``: the older κ/(R + ℓ) curve with the bulk kernel.
* ``sgw``: site-diagonal ΔW on the sBSE monopole RPA. It underscreens CdSe (ε_eff ≈ 1) and is kept
  for reference only.

These models are not consistent in the sense of section 4.

11. Which model
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Purpose
     - Recommendation
   * - Optical spectra, size series
     - ``sgw-resta`` or ``sgw-dim`` (ΔCOHSEX + anchor residual), mnok, in the solvent; or ``gw``
       with ``resta-sphere`` for the cheapest classical estimate. Add SOC for the band edge.
   * - QP gap, IP/EA
     - Any calibrated Delta-W model or ``gw``; all reproduce evGW at the anchor.
   * - Orbital relaxation
     - ``qsgw-*``; S₁ within about 0.1 eV of the one-shot ΔCOHSEX value, at much higher cost for large dots.
   * - Reproduce old results
     - ``qp_selfenergy: classical``, ``qp_anchor_residual: off``, ``qp_residual_scaling: power``,
       ``qp_polarization: legacy``.
