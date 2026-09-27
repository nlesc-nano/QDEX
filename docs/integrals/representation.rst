Representation of the integrals: MNOK or ZDO xs
===============================================

Part of :doc:`/integrals/index`.

The formulas of :doc:`/quasiparticles/gw` and :doc:`/quasiparticles/models` contain integrals of the form

.. math::

   (pq|X|rs) = \iint d\mathbf r\,d\mathbf r'\,\rho_{pq}(\mathbf r)\,X(\mathbf r,\mathbf r')\,\rho_{rs}(\mathbf r'),
   \qquad \rho_{pq} = \psi_p^*\psi_q ,

with X = ΔW for the self-energy. The BSE uses the same integrals with X = W (direct term) and X = v
(exchange term). Computing them exactly needs the full four-index tensor, which is out of reach for
10⁴ basis functions. QDEX offers two approximations, selected by ``integrals.representation``. The
choice applies to ΔW in the QP correction and to W and v in the BSE at the same time.

MNOK: atom-condensed densities (``mnok``)
-----------------------------------------

**Densities.** Each orbital-pair density is condensed onto the atoms:

.. math::

   \rho_{pq}(\mathbf r) \to \sum_A q^{pq}_A\,\delta(\mathbf r - \mathbf R_A),\qquad
   q_A^{pq} = \sum_{\mu\in A}\sum_\nu C_{\mu p}S_{\mu\nu}C_{\nu q}\ \ (\text{Mulliken}),
   \qquad q_A^{pq} = \sum_{\mu\in A} c_{\mu p}c_{\mu q}\ \ (\text{Löwdin}, c = S^{1/2}C).

**Bare interaction.** The Mataga–Nishimoto–Ohno–Klopman (MNOK) family,

.. math::

   \gamma_{AB} = \big(r_{AB}^{\beta} + a_{AB}^{\beta}\big)^{-1/\beta},\qquad
   a_{AB} = \tfrac12\big(\eta_A^{-1} + \eta_B^{-1}\big),\qquad \gamma_{AA} = \eta_A ,

in atomic units. What QDEX uses:

* **Exponent β = 2** (Ohno–Klopman) for every MNOK interaction: the exchange K\ :sup:`x`, the bare γ
  inside the screened W and ΔW, and the direct term. β = 1 is the Mataga–Nishimoto form. There is no
  separate exponent for exchange. (The key ``excitations.kernel_scaling``, historically called
  ``alpha``, is a scale factor of the uniform ``bse`` kernel, not an exponent.)
* **Damping** by the mean of the inverse hardnesses, so that the on-site value is η_A.
* **Hardness** η_A of Ghosh and Islam (``HARDNESS_DICT``), defined as ½ ∂²E/∂N², i.e. half of IP − EA
  (Cd 3.50 eV, Se 5.48 eV). The on-site Coulomb integral of a monopole, (ii|ii) = IP − EA, would be
  2η_A.

Grimme's sTDA uses the same functional form with other choices (:doc:`/excitons/stda`): separate
fitted exponents for exchange (α = 1.42 + 0.48 a_x) and direct term (β = 0.20 + 1.83 a_x), damping by
the mean hardness, and the on-site value 2η_A. These apply only in ``mode: stda``.

**Sensitivity** (Cd₆₈Se₅₅Cl₂₆, vacuum, 25 × 25; changes applied to all MNOK interactions):

.. list-table::
   :header-rows: 1

   * - Choice
     - S₁ sBSE (eV)
     - S₁ ``sgw-resta`` (eV)
     - singlet–triplet, sBSE (meV)
   * - β = 2, γ_AA = η (current)
     - 2.800
     - 3.043
     - 57
   * - β = 1 (Mataga–Nishimoto)
     - 2.854
     - 3.093
     -
   * - β = 2, γ_AA = 2η
     - 2.779
     - 3.038
     - 99

The exponent shifts S₁ by about 0.05 eV; the on-site convention hardly moves S₁ but nearly doubles
the exchange (singlet–triplet) splitting, which is short-range. Neither choice has been validated
against a reference for dots yet.

**Integrals.**

.. math::

   (pq|X|rs) \approx \sum_{AB} q^{pq}_A\,X_{AB}\,q^{rs}_B,\qquad X_{AB} = S_{AB}\,\gamma_{AB}\ \text{or}\ \Delta W_{AB},

with :math:`S_{AB}` the screening factor of the model. For ΔCOHSEX, ΔW_AB is expanded to AO blocks
and contracted with Löwdin coefficients (:doc:`/quasiparticles/models`, section 4).

**Limits.**

* **Monopoles only.** Each atom carries a point charge. Atomic dipoles and higher moments of the
  density are lost, and with them the shape of p and d orbitals and intra-atomic transitions
  (for example s → p on one atom, whose transition charge on that atom is zero).
* **Parameters.** The short range is set by the hardness η_A and the exponent β, not by the basis.
  It is the semi-empirical part of the interaction, and it matters most for exchange (see the
  sensitivity table).
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
