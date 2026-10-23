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

in atomic units, with the on-site value γ_AA = s η_A. QDEX uses (``integrals`` section):

.. list-table::
   :header-rows: 1
   :widths: 30 16 54

   * - Key
     - Default
     - Meaning
   * - ``mnok_exponent``
     - 2
     - β for every MNOK interaction: exchange K\ :sup:`x`, the bare γ in the screened W and ΔW, the
       direct term (2 = Ohno–Klopman, 1 = Mataga–Nishimoto)
   * - ``mnok_exponent_exchange``
     - = ``mnok_exponent``
     - β of the exchange interaction only
   * - ``mnok_onsite``
     - ``ip_ea``
     - on-site value: ``ip_ea`` gives γ_AA = 2η_A = IP − EA, the Coulomb integral of a monopole;
       ``eta`` gives γ_AA = η_A = (IP − EA)/2, the earlier convention

η_A is the hardness of Ghosh and Islam (``HARDNESS_DICT``), defined as ½ ∂²E/∂N² (Cd 3.50 eV,
Se 5.48 eV), and the damping is the mean of the inverse on-site values. (``excitations.kernel_scaling``,
historically ``alpha``, is a scale factor of the uniform ``bse`` kernel, not an exponent.) Grimme's sTDA
uses the same form with its own fitted exponents and 2η (:doc:`/excitons/stda`); those apply only in
``mode: stda``.

**In the output.** Every excited-state run prints the integrals it used, after the QP section:

.. code-block:: text

   --- Two-electron integrals ---
     Representation : MNOK, atom pairs, Mulliken transition charges
     gamma_AB       = (R^beta + a_AB^beta)^(-1/beta),  a_AB = (1/gamma_AA + 1/gamma_BB)/2
     beta (direct)  : 2 (Ohno-Klopman)
     beta (exchange): 2 (Ohno-Klopman)
     On-site        : gamma_AA = IP - EA = 2 eta_A  (integrals.mnok_onsite: ip_ea)
     gamma_AA       : Cd 7.00 eV, Cl 11.73 eV, Se 10.96 eV
     Exchange  K^x  : bare interaction (unscreened)
     Direct    K^d  : W of the QP model (the same W as in the QP correction)  (excitations.kernel: qp)

With ``mode: stda`` the block gives Grimme's γ\ :sup:`J` and γ\ :sup:`K` with the values of a_x, α and
β instead; with ``xs`` it names the exact integrals.

**Benchmark against exact integrals.** The ``xs`` representation evaluates K\ :sup:`x` from exact
(μμ|νν) integrals, with no hardness or exponent, and uses the same distance-dependent screening. It is
the reference for the short range, where MNOK is semi-empirical. The triplet (no K\ :sup:`x`) tests the
direct term, the singlet–triplet splitting ΔST tests exchange. sBSE (bulk Resta W), vacuum,
25 × 25, Löwdin charges as in xs; deviations from xs:

.. list-table::
   :header-rows: 1

   * - MNOK choice
     - 1.2 nm: S₁ / T₁ / ΔST
     - 2.0 nm: S₁ / T₁ / ΔST
   * - xs reference (absolute)
     - 3.246 / 3.043 eV / 203 meV
     - 2.759 / 2.663 eV / 96 meV
   * - β = 2, γ_AA = η (earlier default)
     - +0.44 / +0.49 eV / −50 meV
     - +0.04 / +0.07 eV / −39 meV
   * - β = 1, γ_AA = η
     - +0.56 / +0.61 eV / −55 meV
     - +0.09 / +0.14 eV / −50 meV
   * - β = 2, γ_AA = 2η (default)
     - −0.15 / −0.17 eV / +19 meV
     - +0.01 / −0.01 eV / +13 meV
   * - β = 1, γ_AA = 2η
     - −0.02 / −0.03 eV / +10 meV
     - +0.05 / +0.05 eV / −3 meV
   * - β = 3, γ_AA = 2η
     - −0.20 / −0.22 eV / +23 meV
     - −0.00 / −0.02 eV / +17 meV

* **The on-site value decides exchange.** With η, ΔST is about half the exact value on both clusters;
  with 2η it is within 20 meV. 2η is also the physically correct monopole integral.
* **The exponent tunes the direct term.** With 2η, β = 2 reproduces the 2 nm cluster within 10 meV;
  the 35-atom cluster prefers β = 1. β = 2 is kept as the default for dots of production size.
* **A separate exchange exponent changes ΔST by less than 6 meV** once the on-site value is 2η.
* **Mulliken charges** (the default of ``integrals.charges``) give the same result at 2 nm (ΔST 99 meV)
  but underestimate exchange in the 1.2 nm cluster (ΔST 51 meV); Löwdin charges are closer to xs for
  very small clusters.

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
  It is the semi-empirical part of the interaction; the defaults are chosen against the exact ``xs``
  integrals (see the benchmark above).
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

One-centre extension of xs (``one_centre``)
--------------------------------------------

**Densities.** xs drops every product of two different AOs, also when both sit on the same atom. Those
on-site products (Se 4p·4s, Cd 4d·5p, Pb 6s·6p, …) carry the atomic transition dipoles that dominate the
exchange of the band-edge transitions. ``one_centre: true`` keeps them (NDDO-type densities):

.. math::

   \rho_{pq}(\mathbf r) \to \sum_A\sum_{\mu,\nu\in A} c_{\mu p}c_{\nu q}\,\lambda_\mu(\mathbf r)\lambda_\nu(\mathbf r),
   \qquad \lambda = \chi\,S_{AA}^{-1/2}\ \text{(AOs of each atom orthonormalized among themselves)},

with the Löwdin coefficients :math:`c = S^{1/2}C` of xs. The orthonormalization within each atom keeps the
charges :math:`\sum_\mu c_{\mu p}c_{\mu q}`; products of AOs on different atoms are still dropped.

**Interaction.** The interaction of two one-centre products is exact, :math:`(\lambda_\mu\lambda_\nu|\lambda_\kappa\lambda_\lambda)`
(Libint2, ``compute_onecentre_coulomb``), for the on-site blocks and for atom pairs closer than
``one_centre_cutoff``. Farther apart the xs integrals :math:`(\mu\mu|\nu\nu)` are kept and, with
``one_centre_tail: dipole`` (default), the one-centre densities also interact through their atomic dipoles:

.. math::

   E_{AB} = -q_A\,\mathbf t_{AB}\cdot\mathbf p_B + q_B\,\mathbf t_{AB}\cdot\mathbf p_A
            + \mathbf p_A\cdot\frac{1-3\hat n\hat n}{R_{AB}^3}\cdot\mathbf p_B,\qquad
   \mathbf t_{AB} = \frac{\mathbf R_B-\mathbf R_A}{R_{AB}^3},

   q_A^{pq} = \sum_{\mu\in A} c_{\mu p}c_{\mu q},\qquad
   \mathbf p_A^{pq} = \sum_{\mu,\nu\in A} c_{\mu p}c_{\nu q}\,\langle\lambda_\mu|\mathbf r-\mathbf R_A|\lambda_\nu\rangle .

The charge–charge part is the xs kernel. The dipoles come from the s·p and p·d products on each atom; they
are what the on-site blocks alone leave out at 3–8 Å. The point-multipole form is applied from the nearest neighbours on
(2.5–3 Å, where the densities overlap); against exact integrals it reproduces the full one-centre result
there (table below). ``one_centre_cutoff`` replaces it by exact blocks for the closest pairs if needed. QDEX writes the model as the xs kernel plus a
block-sparse correction

.. math::

   (pq|X|rs) = (pq|X|rs)_{\rm xs} + \sum_{AB,\ R_{AB}\le R_c} D^{pq}_A\cdot\Delta_{AB}\cdot D^{rs}_B,\qquad
   D^{pq}_A[\mu\nu] = \bar c_{\mu p}c_{\nu q},

   \Delta_{AB} = (\lambda_\mu\lambda_\nu|\lambda_\kappa\lambda_\lambda) - \delta_{\mu\nu}\delta_{\kappa\lambda}(\chi_\mu\chi_\mu|\chi_\kappa\chi_\kappa).

Δ (with the dipole tail) is the bare interaction. In the exchange K\ :sup:`x` it completes the bare xs
kernel. In the direct term K\ :sup:`d`, where xs carries the screened W of the kernel, every atom pair of the
correction is screened with the ratio that xs applies to the same atoms,

.. math::

   W_{1c} = W_{\rm xs} + S_{AB}\,\Delta_{AB},\qquad
   S_{AB} = \frac{\sum_{\mu\in A,\nu\in B} W_{\mu\nu}}{\sum_{\mu\in A,\nu\in B}(\mu\mu|\nu\nu)} .

For the Resta-type kernels S\ :sub:`AA` = 1: the on-site blocks stay bare and the dipole tail between atoms is
screened like the xs charges. Without this, the bare tail overbinds the triplet of Cd\ :sub:`16` by 0.4 eV
with the strongly screened Resta W. The BSE products never form pair densities:

.. math::

   (K^x x)_{ia} \mathrel{+}= \Big(T^{\mathsf T}U\bar V\Big)_{ia},\quad U_A=\sum_B\Delta_{AB}:\bar T_B X V_B^{\mathsf T};\qquad
   (K^d x)_{ia} \mathrel{+}= \sum_{AB}\Big(T_A^{\mathsf T}\,\big[\Delta_{AB}:\bar T_A X V_B^{\mathsf T}\big]\,\bar V_B\Big)_{ia},

with T and V the Löwdin coefficients of the active occupied and virtual orbitals (one pair per spin component
for spinors) and X the amplitudes. The spin-free (singlet and triplet), SOC spinor and UKS solvers include it,
in the Davidson products, the dense solver (``full_diag``) and the diagonal BSE; the GPU path adds the
correction on the CPU. Not covered: SOC with a UKS reference (refused), the perturbative selection (its
estimates use the xs couplings; a warning is printed) and the NAMD kernels.

**Accuracy** (exact four-centre integrals over a 30 × 30 window, screened TDA with an RPA W in each
representation, virtual levels scissored to the QP gap; errors in meV):

.. list-table::
   :header-rows: 1

   * - Representation
     - Cd\ :sub:`16`\ Se\ :sub:`13`\ Cl\ :sub:`6`: S₁ / T₁ / ΔST
     - Cs\ :sub:`20`\ Pb\ :sub:`8`\ Br\ :sub:`36`: S₁ / T₁ / ΔST
   * - mnok (Mulliken)
     - −1 / +225 / −226
     - +1089 / +1488 / −398
   * - xs
     - −156 / +9 / −165
     - −269 / −17 / −252
   * - xs + one-centre, on-site
     - −166 / −155 / −11
     - −53 / −85 / +32
   * - xs + one-centre, R_c = 3.2 Å
     - −174 / −155 / −20
     - −8 / −49 / +41
   * - xs + one-centre, R_c = 8 Å
     - −79 / −50 / −29
     - +24 / −45 / +69
   * - **xs + one-centre, on-site + dipole tail (default)**
     - **−41 / −2 / −40**
     - **+65 / −5 / +70**
   * - all one-centre blocks (NDDO)
     - −51 / −21 / −30
     - +68 / +2 / +66

The on-site blocks restore the exchange (ΔST within 10–40 meV instead of 165–250 meV too small). The
absolute S₁ and T₁ also depend on how the on-site densities interact with their neighbours up to ~8 Å; the
dipole tail recovers that (W within 1.8 % of the reference in both clusters, as with all one-centre blocks),
so exact blocks beyond the atom itself are not needed.

**Cost.** One block of (n\ :sub:`A`\ ²)² numbers per atom (Cd 625², 3 MB, stored in two layouts; Se, Cl, Br
169²) and one per atom pair within R\ :sub:`c`, batched by element pair; each BSE product costs
O(n\ :sub:`ao`\ n\ :sub:`occ`\ n\ :sub:`virt`) plus the blocks. The dipole tail adds a 4N\ :sub:`at` × 4N\ :sub:`at`
coupling and, for the Davidson products, the hole and electron moment matrices (4N\ :sub:`at`\ ·(n\ :sub:`occ`\ ² +
n\ :sub:`virt`\ ²) numbers, four times those of the MNOK direct term); the diagonal BSE needs only their
diagonals. The xs kernel is unchanged.

Which to use
------------

* **Large dots (≥ 2 nm):** ``mnok``. The two agree within 0.05 eV and MNOK is much cheaper.
* **Small clusters and molecules,** where the short range dominates: ``xs``. It removes the hardness
  parameters and keeps the AO shape of the densities.
