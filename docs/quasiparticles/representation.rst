Representation of the integrals: MNOK or ZDO xs
===============================================

Part of :doc:`/quasiparticles/index`.

The formulas of :doc:`gw` and :doc:`models` contain integrals of the form

.. math::

   (pq|X|rs) = \iint d\mathbf r\,d\mathbf r'\,\rho_{pq}(\mathbf r)\,X(\mathbf r,\mathbf r')\,\rho_{rs}(\mathbf r'),
   \qquad \rho_{pq} = \psi_p^*\psi_q ,

with X = ΔW for the self-energy. The BSE uses the same integrals with X = W (direct term) and X = v
(exchange term). Computing them exactly needs the full four-index tensor, which is out of reach for
10⁴ basis functions. QDEX offers two approximations, selected by ``two_electron_integrals``. The
choice applies to ΔW in the QP correction and to W and v in the BSE at the same time.

MNOK: atom-condensed densities (``mnok``)
-----------------------------------------

**Densities.** Each orbital-pair density is condensed onto the atoms:

.. math::

   \rho_{pq}(\mathbf r) \to \sum_A q^{pq}_A\,\delta(\mathbf r - \mathbf R_A),\qquad
   q_A^{pq} = \sum_{\mu\in A}\sum_\nu C_{\mu p}S_{\mu\nu}C_{\nu q}\ \ (\text{Mulliken}),
   \qquad q_A^{pq} = \sum_{\mu\in A} c_{\mu p}c_{\mu q}\ \ (\text{Löwdin}, c = S^{1/2}C).

**Bare interaction.** The Mataga–Nishimoto–Ohno–Klopman form, damped so that the on-site value is
the atomic hardness η_A:

.. math::

   \gamma_{AB} = \frac{1}{\sqrt{r_{AB}^2 + a_{AB}^2}},\qquad a_{AB} = \tfrac12\big(\eta_A^{-1} + \eta_B^{-1}\big),
   \qquad \gamma_{AA} = \eta_A .

**Integrals.**

.. math::

   (pq|X|rs) \approx \sum_{AB} q^{pq}_A\,X_{AB}\,q^{rs}_B,\qquad X_{AB} = S_{AB}\,\gamma_{AB}\ \text{or}\ \Delta W_{AB},

with :math:`S_{AB}` the screening factor of the model. For ΔCOHSEX, ΔW_AB is expanded to AO blocks
and contracted with Löwdin coefficients (:doc:`models`, section 4).

**Limits.**

* **Monopoles only.** Each atom carries a point charge. Atomic dipoles and higher moments of the
  density are lost, and with them the shape of p and d orbitals and intra-atomic transitions
  (for example s → p on one atom, whose transition charge on that atom is zero).
* **Parameters.** The short range is set by the hardness η_A, not by the basis. It is the
  semi-empirical part of the interaction.
* **Population dependence.** Mulliken and Löwdin charges differ, most for diffuse basis sets.
  Mulliken populations can be negative.
* **Exact at long range.** For r_AB ≫ a_AB, γ_AB → 1/r_AB. The surface polarization, the solvent term
  and the long-range screening do not depend on the partition.

**Cost.** :math:`N_{\mathrm{at}}^2` storage and contractions; routine for 10⁴ atoms.

ZDO xs: AO density pairs (``xs``)
---------------------------------

**Densities.** Zero differential overlap (ZDO) in the Löwdin-orthogonalized basis: only
products of an AO with itself are kept,

.. math::

   \rho_{pq}(\mathbf r) \to \sum_\mu c_{\mu p}c_{\mu q}\,|\chi_\mu(\mathbf r)|^2,\qquad
   (\mu\nu|\lambda\sigma) \approx \delta_{\mu\nu}\delta_{\lambda\sigma}\,(\mu\mu|\lambda\lambda).

**Bare interaction.** The integrals :math:`(\mu\mu|\nu\nu)` are computed exactly over the contracted
Gaussian basis (Libint2).

**Screened interaction.** The screening of the model is carried over from the atom pairs as a ratio,
and the environment term is added:

.. math::

   W_{\mu\nu} = \frac{W_{A(\mu)B(\nu)}}{\gamma_{A(\mu)B(\nu)}}\,(\mu\mu|\nu\nu) + W^{\mathrm{add}}_{A(\mu)B(\nu)} .

**Integrals.**

.. math::

   (pq|X|rs) \approx \sum_{\mu\nu} c_{\mu p}c_{\mu q}\,X_{\mu\nu}\,c_{\nu r}c_{\nu s}.

QP populations and BSE densities both use Löwdin coefficients in this representation.

**Limits.**

* **ZDO.** Overlap densities :math:`\chi_\mu\chi_\nu` (μ ≠ ν) are dropped: bond charges, and
  intra-atomic transition densities between different AOs. The integrals are evaluated in the
  non-orthogonal basis but contracted with Löwdin coefficients, which is the usual ZDO compromise.
* **Screening borrowed from the atoms.** The ratio W_AB/γ_AB is the same for all AO pairs of two
  atoms. There is no AO-resolved screening.
* **Basis dependence.** The short range follows the basis. Tight functions give larger
  :math:`(\mu\mu|\mu\mu)` than η_A, and diffuse functions smaller.
* **More binding at short range.** Compared with MNOK, S₁ is 0.6 eV lower at 1.2 nm and 0.05 eV lower
  at 2 nm. The two agree at long range.

**Cost.** Several :math:`n_{\mathrm{ao}}^2` matrices (about 1.4 GB each at 13k basis functions) and the
:math:`(\mu\mu|\nu\nu)` integrals.

Which to use
------------

* **Large dots (≥ 2 nm):** ``mnok``. The two agree within 0.05 eV and MNOK is much cheaper.
* **Small clusters and molecules,** where the short range dominates: ``xs``. It removes the hardness
  parameters and keeps the AO shape of the densities.
* **Anchor residuals** are calibrated per representation. Changing it needs no action, but the
  calibration exists only for the combinations in ``qdex/data/dw_anchor_residuals.json``.
