Part 3: Quasiparticle Corrections & The Scaled GW Models
=========================================================

Standard semi-local Kohn-Sham Density Functional Theory (DFT) using functionals like PBE severely underestimates the fundamental band gap of semiconductor nanostructures. For instance, PBE predicts a band gap of :math:`\approx 1.5\text{ eV}` for bulk :math:`\text{CsPbBr}_3`, whereas the experimental quasiparticle gap is :math:`\approx 2.35\text{ eV}`.

In ``QDEX``, single-particle excitation energies are corrected via an analytical and physically grounded hierarchy of **Scaled GW Quasiparticle Models** that incorporate quantum confinement, microscopic electrostatic screening, dielectric solvation, dynamic renormalization, and orbital relaxation.

---

1. Quasiparticle Theory and Hedin's GW Approximation
-----------------------------------------------------

In many-body perturbation theory, the true single-particle electron addition (electron affinity, EA) and removal (ionization potential, IP) energies are determined by the **quasiparticle (QP) equation**:

.. math::

   \left( \hat{T} + \hat{V}_{\mathrm{ext}}(\mathbf{r}) + \hat{V}_H(\mathbf{r}) \right) \psi_k(\mathbf{r}) + \int \Sigma(\mathbf{r}, \mathbf{r}'; \varepsilon_k^{\mathrm{QP}}) \psi_k(\mathbf{r}') \, d\mathbf{r}' = \varepsilon_k^{\mathrm{QP}} \psi_k(\mathbf{r})

where:

* :math:`\hat{V}_H(\mathbf{r})` is the classical Hartree potential.
* :math:`\Sigma(\mathbf{r}, \mathbf{r}'; \omega)` is the non-local, energy-dependent **electron self-energy** operator.

Hedin's Equations & The :math:`G_0W_0` Approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Hedin's pentagon of many-body equations, the self-energy is expanded to first order in the dynamically screened Coulomb interaction :math:`W`:

.. math::

   \Sigma(\mathbf{r}, \mathbf{r}'; \omega) = \frac{i}{2\pi} \int e^{i \omega' \eta} G_0(\mathbf{r}, \mathbf{r}'; \omega + \omega') W_0(\mathbf{r}, \mathbf{r}'; \omega') \, d\omega'

where:

* :math:`G_0` is the non-interacting single-particle Green's function.
* :math:`W_0 = \epsilon^{-1} v` is the screened Coulomb interaction, mediated by the dynamic microscopic dielectric function :math:`\epsilon(\mathbf{r}, \mathbf{r}'; \omega)`.

The quasiparticle energy shift :math:`\Delta \varepsilon_k^{\mathrm{QP}} = \varepsilon_k^{\mathrm{QP}} - \varepsilon_k^{\mathrm{DFT}}` is given to first order by:

.. math::

   \Delta \varepsilon_k^{\mathrm{QP}} = Z_k \, \operatorname{Re}\langle \psi_k | \Sigma(\varepsilon_k^{\mathrm{DFT}}) - v_{xc}^{\mathrm{DFT}} | \psi_k \rangle

where :math:`Z_k = \left( 1 - \left. \frac{\partial \operatorname{Re}\Sigma}{\partial \omega} \right|_{\varepsilon_k} \right)^{-1}` is the quasiparticle renormalization weight.

---

2. Unified Two-Body Interaction Engine (``2e-integrals`` vs. ``kernel``)
------------------------------------------------------------------------

In many-body perturbation theory, both the single-particle electron self-energy :math:`\Sigma = iGW` (this Part) and the two-particle electron-hole interaction kernel :math:`\Xi = \delta \Sigma / \delta G \approx v - W` (:doc:`../part4_excited_states/index`) are fundamentally constructed from the **exact same pair of two-body operators**:

1. **The Bare Coulomb Repulsion Operator**:

   .. math::

      v(\mathbf{r}_1, \mathbf{r}_2) = \frac{1}{|\mathbf{r}_1 - \mathbf{r}_2|}

2. **The Dynamically Screened Coulomb Interaction Operator**:

   .. math::

      W(\mathbf{r}_1, \mathbf{r}_2; \omega) = \int \epsilon^{-1}(\mathbf{r}_1, \mathbf{r}_3; \omega) \, v(\mathbf{r}_3, \mathbf{r}_2) \, d\mathbf{r}_3

In ``QDEX``, rather than coupling these operators into an opaque, monolithic choice, they are decoupled into two orthogonal, physically transparent configuration keywords:

.. list-table::
   :widths: 22 28 50
   :header-rows: 1

   * - Configuration Keyword
     - Options
     - Mathematical & Physical Meaning
   * - ``2e-integrals``
     - ``mnok``, ``xs`` (or ``xs-qdex``)
     - **Bare Coulomb Representation**: Chooses how :math:`v(\mathbf{r}_1, \mathbf{r}_2)` is evaluated (semi-empirical atom-centered point charges vs. exact analytical 4-center Gaussian integrals).
   * - ``kernel``
     - ``resta``, ``dim``, ``rpa``, ``sbse``, ``bse``
     - **Dielectric Screening Profile**: Chooses how :math:`\epsilon^{-1}` and :math:`W(\mathbf{r}_1, \mathbf{r}_2)` are computed across the nanocluster and its surrounding environment.

Shared Two-Body Operator Architecture
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The diagram below illustrates how this unified engine simultaneously powers Quasiparticle Theory and Optical Excitations:

.. code-block:: text

   +-----------------------------------------------------------------------------------------+
   |                          Unified Two-Body Interaction Engine                            |
   |                                                                                         |
   |   Bare Coulomb Operator v                        Screened Interaction W = ε⁻¹ v         |
   |   (keyword: '2e-integrals')                      (keyword: 'kernel')                    |
   |   • mnok (Ohno damped point charges)             • resta (Thomas-Fermi + Penn scaling)  |
   |   • xs   (Libint2 4-center Gaussian integrals)   • dim   (Atomistic polarizable dipoles)|
   |                                                  • rpa   (Microscopic transition polar.)|
   |                                                  • sbse  (Simplified polarizability)    |
   |                                                  • bse   (Uniform bulk ε_∞ scaling)     |
   +----------------------------+-----------------------------------+------------------------+
                                |                                   |
                                v                                   v
   +--------------------------------------------+   +----------------------------------------+
   |   Part 3: Quasiparticle (GW) Corrections   |   |  Part 4: Optical Excitations (BSE)    |
   |                                            |   |                                        |
   |   • Self-energy: Σ = i G W                 |   |   • Bethe-Salpeter Kernel: Ξ = v - W   |
   |   • Shift: ΔΣ = Σ_COH[ΔW] + Σ_SEX[ΔW]      |   |   • Bare exchange: K_x(v)              |
   |   • Approach B: Microscopic Asymmetry      |   |   • Direct attraction: K_d(W)          |
   |     HOMO/LUMO split (f_H, f_L) from MOs    |   |   • Wannier-Mott bulk limit as R -> ∞  |
   +--------------------------------------------+   +----------------------------------------+

Role in Quasiparticle Theory vs. Excited States
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 35 40
   :header-rows: 1

   * - Theoretical Dimension
     - Role in Part 3 (Quasiparticles: :math:`\Sigma = iGW`)
     - Role in Part 4 (Excitations: :math:`\Xi = v - W`)
   * - **Bare Repulsion (``2e-integrals``)**
     - Evaluates the baseline atomic interaction matrix :math:`\mathbf{\Gamma}` used in screening polarization and on-site charging.
     - Evaluates the bare exchange matrix :math:`K^x`, determining singlet-triplet exchange splitting.
   * - **Screened Interaction (``kernel``)**
     - Computes the confinement shift :math:`\Delta W = W^{\mathrm{QD}} - W^{\mathrm{bulk}} + W^{\mathrm{solv}}` for COHSEX self-energies and :math:`qsGW` Hamiltonian updates.
     - Computes the direct screened attraction :math:`K^d`, binding the electron and hole into an exciton.
   * - **Short-Range Limit (:math:`r \to 0`)**
     - :math:`W \to v = 1/r` (:math:`\epsilon \to 1`), preventing unphysical damping of on-site self-energy shifts.
     - Preserves full on-site Coulomb repulsion and atomic exchange splitting.
   * - **Bulk Limit (:math:`R_{\mathrm{QD}} \to \infty`)**
     - :math:`\Delta W \to 0`, naturally recovering the bulk quasiparticle band gap :math:`\Delta_{\mathrm{bulk}}`.
     - :math:`-K^d \to -e^2/(\epsilon_\infty r)`, recovering the bulk Wannier-Mott exciton binding energy :math:`E_b^{\mathrm{bulk}} = R_y^*`.

---

3. The DFT Band Gap Problem
---------------------------

The failure of Kohn-Sham DFT with semi-local functionals (LDA, PBE) to predict fundamental band gaps arises from two fundamental physical sources:

1. **Missing Derivative Discontinuity (:math:`\Delta_{xc}`)**:
   The true fundamental band gap is:

   .. math::

      E_g = \mathrm{IP} - \mathrm{EA} = E_g^{\mathrm{KS}} + \Delta_{xc}

   where :math:`\Delta_{xc} = \left. \frac{\delta E_{xc}}{\delta \rho} \right|_{N+\delta} - \left. \frac{\delta E_{xc}}{\delta \rho} \right|_{N-\delta}` is the integer discontinuity of the exchange-correlation functional. In standard LDA/GGA functionals, :math:`E_{xc}[\rho]` is a continuous function of density, identically yielding :math:`\Delta_{xc} = 0`. Consequently, :math:`E_g^{\mathrm{PBE}}` underestimates the true quasiparticle gap by 30–50%.

2. **Self-Interaction Error (SIE)**:
   In semi-local functionals, an electron spuriously interacts with its own charge density through the Hartree term. This unphysical electrostatic repulsion artificially destabilizes occupied states and over-delocalizes frontier wavefunctions.

---

4. The Nanocrystal Scaling Bottleneck
-------------------------------------

Although full :math:`G_0W_0` calculations resolve the band gap problem for small molecules or bulk crystals, computing :math:`G_0W_0` explicitly for colloidal quantum dots is computationally prohibitive:

* **Computational Complexity**: Standard :math:`G_0W_0` scales as :math:`O(N^4)` with real-space basis sets and :math:`O(N^5)` with plane-wave expansions.
* **Memory Footprint**: Calculating the polarizability matrix :math:`\chi_0(\mathbf{r}, \mathbf{r}'; \omega)` requires summing over thousands of unoccupied conduction states, demanding terabytes of RAM for clusters with :math:`> 500` atoms.
* **Nanocrystal Realities**: Chemically realistic colloidal quantum dots comprise 1,000 to 10,000 atoms, including surface passivation ligands, rendering direct *ab initio* :math:`G_0W_0` impossible for high-throughput screening or non-adiabatic molecular dynamics trajectories.

---

5. Hierarchy of Quasiparticle Models in QDEX
--------------------------------------------

To bypass this bottleneck across diverse computational regimes, ``QDEX`` provides three complementary avenues of quasiparticle models:

.. list-table::
   :widths: 22 23 55
   :header-rows: 1

   * - Avenue
     - Method Identifiers
     - Theoretical Characterization
   * - **Two-Anchor Scaled GW**
     - ``sgw-anchor`` (alias ``gw``)
     - Analytic size-interpolation between vacuum cluster anchor (:math:`R_0`) and bulk limit (:math:`R \to \infty`) with classical image-charge solvation. Sub-millisecond evaluation.
   * - **Microscopic Dielectric Shift (:math:`\Delta W`)**
     - ``sgw-dim``, ``sgw-resta``, ``sgw``
     - Microscopic screened Coulomb difference :math:`\Delta W = W^{\mathrm{QD}} - W^{\mathrm{bulk}} + W^{\mathrm{solv}}`. Local :math:`v_{xc}` cancels identically; parameter-free.
   * - **Self-Consistent Extensions**
     - ``evgw-dim``, ``evgw-resta``, ``qsgw-dim``, ``qsgw-resta``
     - Iterative eigenvalue self-consistency (:math:`evGW`) and full AO-basis orbital relaxation (:math:`qsGW`) with dynamic plasmon-pole renormalization :math:`Z_p`.

---

6. Avenue 1: Two-Anchor Scaled GW (``sgw-anchor``)
--------------------------------------------------

Instead of guessing empirical parameters, ``sgw-anchor`` (historically called ``gw``) anchors the quasiparticle gap between two rigorously computed physical limits:

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Anchor Limit
     - Physical System
     - Theoretical Characterization
   * - **Anchor 1: Smallest Vacuum Anchor (:math:`R_0`)**
     - Monomer or smallest stoichiometric Wulff cluster
     - Relaxed cluster computed in vacuum with hybrid DFT (:math:`\text{PBE0}`) and eigenvalue-self-consistent GW (``EV_GW_ITER 4``).
   * - **Anchor 2: Bulk Limit (:math:`R \to \infty`)**
     - Periodic crystal
     - High-accuracy bulk :math:`G_0W_0` quasiparticle gap (:math:`E_g^{\mathrm{GW, bulk}}`), calibrated against experimental ARPES.

The Confinement Interpolation Formula
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The radius :math:`R` is the equivalent-volume radius of the inorganic core, :math:`R = (3V/4\pi)^{1/3}`. The scissor added to a PBE eigenvalue difference is:

.. math::

   \Delta_{\mathrm{GW}}(R)
   = \Delta_{\mathrm{bulk}}
   + \frac{\kappa_{\mathrm{out}}}{R+\ell}
   + A\left(\frac{R_0}{R}\right)^{p}

with the bulk shift defined by:

.. math::

   \Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW,bulk}} - E_g^{\mathrm{PBE,bulk}}.

The dielectric prefactor is:

.. math::

   \kappa_{\mathrm{out}} = 11.52\,\mathrm{eV\,\AA}\left(\frac{1}{\epsilon_{\mathrm{out}}} - \frac{1}{\epsilon_\infty}\right).

In vacuum (:math:`\epsilon_{\mathrm{out}}=1`), :math:`\kappa_{\mathrm{vac}} = 11.52\,\mathrm{eV\,\AA}\,(1 - 1/\epsilon_\infty)`. The standard defaults are regularization length :math:`\ell = 1.0\,\text{Å}` and power :math:`p = 2`.

The anchor amplitude :math:`A` is fixed by the finite cluster opening:

.. math::

   A = \Delta(R_0) - \Delta_{\mathrm{bulk}} - \frac{\kappa_{\mathrm{vac}}}{R_0+\ell}

where :math:`\Delta(R_0) = E_g^{\mathrm{GW,cluster}} - E_g^{\mathrm{PBE,cluster}}`.

Implementation in QDEX (``sgw-anchor``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-anchor scaled GW model is implemented in :func:`qdex.hardness.estimate_anchor_scaled_qp_gap`:

1. **Nanocrystal Core Sizing**: Evaluates the equivalent spherical core radius :math:`R = (3V/4\pi)^{1/3}` or gyro-radius from the 3D atomic coordinates :math:`\mathbf{R}_A` of the inorganic core.
2. **Database Query**: Queries tabulated bulk and monomer anchor parameters from ``MATERIAL_DB`` (:math:`E_g^{\mathrm{PBE, bulk}}, E_g^{\mathrm{GW, bulk}}, R_0, \Delta(R_0), f_{\mathrm{homo}}, f_{\mathrm{lumo}}`).
3. **Electrostatic Image Charges**: Computes dielectric prefactors :math:`\kappa_{\mathrm{vac}} = 11.52 (1 - 1/\epsilon_\infty)` and :math:`\kappa_{\mathrm{out}} = 11.52 (1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty)` in :math:`\mathrm{eV\cdot\mathring{A}}`.
4. **Analytic Scissor Calculation**: Evaluates the total opening :math:`\Delta_{\mathrm{GW}}(R) = \Delta_{\mathrm{bulk}} + \frac{\kappa_{\mathrm{out}}}{R+\ell} + A (R_0/R)^p`.
5. **Level Alignment & Provenance**: Partitions :math:`\Delta_{\mathrm{GW}}` across occupied and virtual manifolds using database fractions :math:`f_{\mathrm{homo}}, f_{\mathrm{lumo}}` (Approach A) and stores detailed diagnostics in the returned ``details`` dictionary.

*CLI & YAML Invocation*:

.. code-block:: bash

   qdex --mos ground_state.mos --material CSPBBR3 --qp_gap sgw-anchor --eps-out 2.25

---

7. Avenue 2: Microscopic Dielectric Shift Model (:math:`\Delta W`)
------------------------------------------------------------------

When finite anchor clusters are unavailable or when state-specific orbital self-energies are required, the **Microscopic :math:`\Delta W` Model** computes the quasiparticle shift directly from the difference between the confined nanocrystal screened interaction :math:`W^{\mathrm{QD}}` and the bulk crystal screened interaction :math:`W^{\mathrm{bulk}}`.

The Physical Rationale: Cancellation of :math:`v_{xc}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Hedin's :math:`G_0W_0` quasiparticle equation, the energy shift relative to DFT is:

.. math::

   \Delta \varepsilon_p^{\mathrm{QP}} = Z_p \, \langle \psi_p | \Sigma - v_{xc}^{\mathrm{DFT}} | \psi_p \rangle.

In semiconductors, the short-range exchange-correlation potential :math:`v_{xc}^{\mathrm{DFT}}(\mathbf{r})` is primarily determined by local atomic density and core-valence overlap, which is nearly identical between the quantum dot interior and the bulk crystal. By referencing the quantum dot self-energy to the known bulk quasiparticle correction :math:`\Delta \varepsilon_p^{\mathrm{bulk}} = Z^{\mathrm{bulk}} \langle \psi_p | \Sigma^{\mathrm{bulk}} - v_{xc}^{\mathrm{DFT}} | \psi_p \rangle`, the local potential :math:`v_{xc}^{\mathrm{DFT}}` cancels out identically:

.. math::

   \Delta \Sigma(\mathbf{r}, \mathbf{r}') = \Sigma^{\mathrm{QD}}(\mathbf{r}, \mathbf{r}') - \Sigma^{\mathrm{bulk}}(\mathbf{r}, \mathbf{r}') \approx \Delta \Sigma^{\mathrm{COH}}(\mathbf{r}, \mathbf{r}') + \Delta \Sigma^{\mathrm{SEX}}(\mathbf{r}, \mathbf{r}').

Screened COHSEX Operator
~~~~~~~~~~~~~~~~~~~~~~~~

In the static Coulomb-Hole plus Screened-Exchange (COHSEX) approximation:

.. math::

   \Delta W(\mathbf{r}, \mathbf{r}') = W^{\mathrm{QD}}(\mathbf{r}, \mathbf{r}') - W^{\mathrm{bulk}}(\mathbf{r}, \mathbf{r}') + W^{\mathrm{solv}}(\mathbf{r}, \mathbf{r}')

where:

1. **Screened Exchange (:math:`\Delta \Sigma^{\mathrm{SEX}}`)**:

   .. math::

      \Delta \Sigma^{\mathrm{SEX}}(\mathbf{r}, \mathbf{r}') = -\rho(\mathbf{r}, \mathbf{r}') \, \Delta W(\mathbf{r}, \mathbf{r}')

2. **Coulomb Hole (:math:`\Delta \Sigma^{\mathrm{COH}}`)**:

   .. math::

      \Delta \Sigma^{\mathrm{COH}}(\mathbf{r}, \mathbf{r}') = \frac{1}{2} \, \delta(\mathbf{r} - \mathbf{r}') \, \Delta W(\mathbf{r}, \mathbf{r}).

Projected onto state :math:`p`:

.. math::

   \Delta \varepsilon_p^{\mathrm{conf}} = \frac{Z_p}{2} \sum_{A=1}^{N_{\mathrm{atoms}}} q_A^p \, \Delta W_{AA} - Z_p \sum_{A,B} q_A^p \, \Delta W_{AB} \, q_B^{\mathrm{occ}}

where :math:`q_A^p = \sum_{\mu \in A, \nu} C_{\mu p} S_{\mu \nu} C_{\nu p}` is the Mulliken or Lowdin atomic orbital population of state :math:`p`.

The total quasiparticle scissor is:

.. math::

   \Delta_{\mathrm{sGW}} = \Delta_{\mathrm{bulk}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}} - \Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}}.

Screening Formulations for :math:`W^{\mathrm{QD}}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``QDEX`` implements three distinct physical methods for evaluating :math:`W^{\mathrm{QD}}`:

1. **Discrete Dipole Interaction Model (``sgw-dim``)**:
   Models the semiconductor lattice as an array of atom-centered polarizable point dipoles :math:`\boldsymbol{\mu}_A = \alpha_A \mathbf{E}^{\mathrm{tot}}(\mathbf{R}_A)`.
   The total local field includes external fields and self-consistent dipolar interactions:

   .. math::

      \boldsymbol{\mu}_A = \alpha_A \left[ \mathbf{E}_0(\mathbf{R}_A) - \sum_{B \neq A} \mathbf{T}_{AB} \boldsymbol{\mu}_B \right]

   where :math:`\mathbf{T}_{AB} = \frac{3 \mathbf{R}_{AB} \otimes \mathbf{R}_{AB} - R_{AB}^2 \mathbf{I}}{R_{AB}^5} f_{\mathrm{Thole}}(R_{AB})` is the Thole-damped dipole interaction tensor.
   Solving the linear system :math:`(\mathbf{I} + \boldsymbol{\alpha}\mathbf{T}) \boldsymbol{\mu} = \boldsymbol{\alpha}\mathbf{E}_0` yields the microscopic screened interaction :math:`W_{AB}^{\mathrm{QD}}`.

2. **Resta-Penn Screened Dielectric Model (``sgw-resta``)**:
   Evaluates distance-dependent electronic screening based on the Thomas-Fermi model of the valence electron gas:

   .. math::

      W_{AB}^{\mathrm{QD}} = \frac{1}{\epsilon_{\mathrm{in}}(R) r_{AB}} + \frac{1 - \epsilon_{\mathrm{in}}(R)^{-1}}{r_{AB}} \exp\left( -\frac{r_{AB}}{\lambda_{\mathrm{TF}}} \right)

   where quantum confinement suppresses the core dielectric constant according to the Penn model:

   .. math::

      \epsilon_{\mathrm{in}}(R) = 1 + (\epsilon_\infty - 1) \frac{1}{1 + (R_p / R)^2}, \quad R_p = \frac{\pi}{2} \left( \frac{13.6\text{ eV}}{E_g^{\mathrm{bulk}}} \right) a_0.

3. **Simplified BSE Screened Interaction (``sgw``)**:
   Implements the Cho, Bintrim, and Berkelbach [JCTC 18, 3054 (2022)] polarizability kernel:

   .. math::

      W = (\mathbf{I} + \mathbf{J}_{\mathrm{solv}} \boldsymbol{\Pi}^0)^{-1} \mathbf{J}_{\mathrm{solv}}

   where :math:`\boldsymbol{\Pi}^0` is the non-interacting transition polarizability matrix.

Spatial Asymptotics of the Dielectric Kernel: Why Screening Fits Nanocrystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The interaction kernel :math:`W(\mathbf{r}, \mathbf{r}')` exhibits three essential spatial regimes that make it uniquely suited for semiconductor nanocrystals:

1. **Short-Range Limit (:math:`r \to 0`, intra-atomic / on-site)**:
   At sub-nanometer distances, the valence electron gas cannot displace fast enough to screen charge fluctuations without violating quantum kinetic energy constraints (Pauli exclusion / Thomas-Fermi cutoff). Therefore, the screening function satisfies:

   .. math::

      \lim_{r \to 0} W(r) = v(r) = \frac{1}{r}, \quad (\epsilon \to 1).

   Standard macroscopic dielectric approximations (:math:`W = v / \epsilon_\infty`) artificially divide on-site Coulomb repulsion by :math:`\epsilon_\infty \approx 6 - 10`, severely underestimating atomic charging energies and exciton binding energies. The Resta and DIM kernels preserve the true unscreened Coulomb interaction at short range.

2. **Intermediate Confinement & Solvent Boundary (:math:`r \sim R_{\mathrm{QD}}`)**:
   Within the nanocrystal core, quantum confinement suppresses electronic polarizability, reducing the effective internal permittivity :math:`\epsilon_{\mathrm{in}}(R) < \epsilon_\infty`. At the nanocrystal-solvent interface, the dielectric mismatch with the surrounding medium (:math:`\epsilon_{\mathrm{out}}`) induces image charges that screen carriers through the reaction field :math:`W^{\mathrm{solv}}`.

3. **Asymptotic Bulk Limit (:math:`r \to \infty` or :math:`R_{\mathrm{QD}} \to \infty`)**:
   Across large distances inside the crystal, dielectric polarization fully develops:

   .. math::

      \lim_{r \to \infty} W(r) = \frac{1}{\epsilon_\infty r}.

   Crucially, as the nanocrystal size grows (:math:`R_{\mathrm{QD}} \to \infty`), the cluster screening profile converges to the periodic crystal limit:

   .. math::

      \lim_{R_{\mathrm{QD}} \to \infty} W^{\mathrm{QD}}(\mathbf{r}, \mathbf{r}') = W^{\mathrm{bulk}}(\mathbf{r}, \mathbf{r}').

   Consequently, the net confinement shift :math:`\Delta W = W^{\mathrm{QD}} - W^{\mathrm{bulk}} + W^{\mathrm{solv}} \to 0` vanishes identically for large crystals, naturally recovering the bulk quasiparticle band gap and bulk exciton spectrum.

Approach B: Microscopic Wavefunction Asymmetry (:math:`f_H^{\mathrm{micro}}, f_L^{\mathrm{micro}}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A major advantage of the microscopic :math:`\Delta W` models (``sgw-dim``, ``sgw-resta``, ``evgw``, ``qsgw``) over the two-anchor model is that they **completely skip the need for calibrated monomer anchor clusters**.

Instead of inheriting empirical fractions from a database, the state-specific self-energies are computed directly from the actual frontier wavefunctions:

.. math::

   \Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} = Z_H \, \sigma_H^{\mathrm{stat}} = \frac{1}{2} Z_H \sum_{A, B=1}^{N_{\mathrm{atoms}}} q_A^{\mathrm{HOMO}} \, \Delta W_{AB} \, q_B^{\mathrm{HOMO}}

.. math::

   \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}} = Z_L \, \sigma_L^{\mathrm{stat}} = \frac{1}{2} Z_L \sum_{A, B=1}^{N_{\mathrm{atoms}}} q_A^{\mathrm{LUMO}} \, \Delta W_{AB} \, q_B^{\mathrm{LUMO}}.

Because the HOMO (typically anion :math:`p` orbitals) and LUMO (typically cation :math:`s` orbitals) possess distinct spatial delocalizations and chemical hardnesses, their dielectric self-energies are naturally asymmetric. ``QDEX`` defines the **microscopic asymmetry fractions**:

.. math::

   f_H^{\mathrm{micro}} = \frac{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}}}{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}, \quad
   f_L^{\mathrm{micro}} = \frac{\Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}

subject to :math:`f_H^{\mathrm{micro}} + f_L^{\mathrm{micro}} = 1.0`.

The bulk shift :math:`\Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW, bulk}} - E_g^{\mathrm{PBE, bulk}}` is then partitioned using these exact same microscopic fractions:

.. math::

   \varepsilon_i^{\mathrm{QP}} = \varepsilon_i^{\mathrm{DFT}} - f_H^{\mathrm{micro}} \Delta_{\mathrm{total}} \quad (i \in \mathrm{occ})

.. math::

   \varepsilon_a^{\mathrm{QP}} = \varepsilon_a^{\mathrm{DFT}} + f_L^{\mathrm{micro}} \Delta_{\mathrm{total}} \quad (a \in \mathrm{virt})

where :math:`\Delta_{\mathrm{total}} = \Delta_{\mathrm{bulk}} + \Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}`.

This provides:

1. **Anchor-Free Evaluation**: Can be applied to any material, core/shell geometry, or surface ligand shell without requiring pre-computed vacuum Wulff clusters.
2. **Consistent Level Alignment**: Accurately predicts absolute Ionization Potentials (IP) and Electron Affinities (EA) tailored to the specific nanocrystal shape and surface termination.
3. **Rigorous qsGW Background**: Partitions the bulk reference Hamiltonian :math:`\mathbf{H}_{\mathrm{bulk}} = -f_H^{\mathrm{micro}} \Delta_{\mathrm{bulk}} (0.5 \mathbf{P}_{\mathrm{occ}}) + f_L^{\mathrm{micro}} \Delta_{\mathrm{bulk}} \mathbf{Q}_{\mathrm{virt}}`.

Implementation in QDEX (``sgw-dim``, ``sgw-resta``, ``sgw``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The microscopic :math:`\Delta W` algorithms are implemented in :mod:`qdex.hardness`:

1. **Discrete Dipole Interaction Model (``sgw-dim``)**:
   Implemented in :func:`qdex.hardness.estimate_sgw_dim_qp_gap`:
   - Builds the :math:`3N \times 3N` Thole dipole interaction matrix :math:`\mathbf{M} = \boldsymbol{\alpha}^{-1} + \mathbf{T}` from the polarizability table ``POLARIZABILITY_TABLE_AU``.
   - Inverts :math:`\mathbf{M}` to solve for induced dipoles under 3D probe fields, obtaining atom-specific effective polarizabilities :math:`p_A^{\mathrm{eff}}`.
   - Computes pairwise screening factors :math:`S_{AB}` and assembles the screened potential :math:`W_{AB}^{\mathrm{DIM}} = S_{AB} \gamma_{AB}^{\mathrm{Ohno}}`.
   - Evaluates the confinement contrast matrix :math:`\Delta W_{AB} = \max(0, W_{AB}^{\mathrm{DIM}} - W_{AB}^{\mathrm{bulk}}) + W_{AB}^{\mathrm{solv}}`.
   - Contracts with Löwdin atomic populations :math:`\mathbf{q}_H, \mathbf{q}_L` to evaluate static shifts :math:`\sigma_H^{\mathrm{stat}}, \sigma_L^{\mathrm{stat}}`.
   - Evaluates :math:`f_H^{\mathrm{micro}}, f_L^{\mathrm{micro}}` (Approach B) and exports them in the ``provenance`` dictionary.

2. **Resta Screened Dielectric Model (``sgw-resta``)**:
   Implemented in :func:`qdex.hardness.estimate_sgw_resta_qp_gap`:
   - Determines the median core nearest-neighbor bond distance :math:`d_{\mathrm{NN}}`.
   - Evaluates the Penn-scaled core dielectric constant :math:`\epsilon_{\mathrm{eff}}(R)`.
   - Computes the Thomas-Fermi screening wavevector :math:`k_s = \sqrt{\epsilon_{\mathrm{eff}} - 1} / d_{\mathrm{NN}}`.
   - Assembles the Resta screened potential :math:`W_{AB}^{\mathrm{Resta}} = [c_\infty + (1 - c_\infty) e^{-k_s R_{AB}}] \gamma_{AB}^{\mathrm{Ohno}}`.
   - Computes :math:`\Delta W` and evaluates state shifts and Approach B splitting fractions.

3. **Simplified BSE Polarizability Model (``sgw``)**:
   Implemented in :func:`qdex.hardness.estimate_sgw_qp_gap`:
   - Evaluates the non-interacting transition polarizability matrix :math:`\boldsymbol{\Pi}^0` from active Kohn-Sham orbital pairs.
   - Evaluates the dynamic inversion :math:`\mathbf{W} = (\mathbf{I} + \mathbf{J}_{\mathrm{solv}} \boldsymbol{\Pi}^0)^{-1} \mathbf{J}_{\mathrm{solv}}`.

*CLI & YAML Invocation*:

.. code-block:: bash

   # DIM polarizable dipoles (recommended default)
   qdex --mos ground_state.mos --material CDSE --qp_gap sgw-dim --dynamic_z --eps-out 2.40

   # Resta Thomas-Fermi screening
   qdex --mos ground_state.mos --material CDSE --qp_gap sgw-resta --dynamic_z --eps-out 2.40

---

8. Avenue 3: Dynamic Renormalization & Self-Consistency
-------------------------------------------------------

Dynamic Renormalization Factor :math:`Z_p` (``--dynamic_z``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than assuming a fixed empirical renormalization factor (such as :math:`Z \approx 0.8`), ``QDEX`` can compute state-dependent dynamic weights :math:`Z_p` from the **Plasmon-Pole Model (PPM) constrained by the :math:`f`-sum rule**:

.. math::

   Z_p = \left( 1 - \left. \frac{\partial \operatorname{Re}\Sigma_p}{\partial \omega} \right|_{\varepsilon_p} \right)^{-1} = \left( 1 + \frac{\Delta \Sigma_p^{\mathrm{stat}}}{\tilde{\Omega}_p} \right)^{-1}

where :math:`\Delta \Sigma_p^{\mathrm{stat}}` is the static COHSEX self-energy shift and :math:`\tilde{\Omega}_p` is the screened plasmon frequency:

.. math::

   \tilde{\Omega}_p = \sqrt{ \frac{\Omega_p^2}{\max(1.0, \epsilon_{\mathrm{eff}} - 1.0)} + E_g^2 }.

Here :math:`\Omega_p = \sqrt{4\pi n_v e^2 / m_e}` is the valence electron plasmon energy (:math:`\approx 15 - 20\text{ eV}`). In large quantum dots with strong screening (:math:`\epsilon_{\mathrm{eff}} \gg 1`), :math:`Z_p \to 0.80 - 0.85`; in ultra-small clusters with suppressed screening, :math:`Z_p \to 0.90 - 0.95`.

Implementation in QDEX (``--dynamic_z``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dynamic state-dependent renormalization factors are evaluated by :func:`qdex.hardness.compute_dynamic_z`:

- Retrieves material-specific valence plasmon energy :math:`\Omega_p \approx 15 - 20\text{ eV}` from ``MATERIAL_PLASMON_ENERGIES``.
- Evaluates the screened plasmon frequency :math:`\tilde{\Omega}_p = \sqrt{\frac{\Omega_p^2}{\max(1.0, \epsilon_{\mathrm{eff}} - 1.0)} + E_g^2}` from the current gap and effective screening constant :math:`\epsilon_{\mathrm{eff}}`.
- Computes state renormalization factors :math:`Z_H = (1 + \sigma_H^{\mathrm{stat}} / \tilde{\Omega}_p)^{-1}` and :math:`Z_L = (1 + \sigma_L^{\mathrm{stat}} / \tilde{\Omega}_p)^{-1}`.
- Clamps values to the physically robust interval :math:`Z \in [0.50, 0.98]` to safeguard against numerical unphysicalities in ultra-small clusters.

Eigenvalue Self-Consistent GW (``evgw-dim``, ``evgw-resta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard one-shot :math:`G_0W_0` depends on the starting DFT functional. The **Eigenvalue Self-Consistent GW (:math:`evGW`)** model eliminates this bias by updating the quasiparticle eigenvalues in the Green's function iteratively:

.. math::

   \varepsilon_p^{(n+1)} = \varepsilon_p^{\mathrm{DFT}} + \Delta \varepsilon_p^{\mathrm{QP}}\big( \{\varepsilon_q^{(n)}\} \big)

until :math:`|\varepsilon_p^{(n+1)} - \varepsilon_p^{(n)}| < 10^{-4}\text{ eV}` (typically 3–5 iterations).

Implementation in QDEX (``evgw-dim``, ``evgw-resta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eigenvalue self-consistency is implemented in :func:`qdex.hardness.estimate_evgw_dim_qp_gap` and :func:`qdex.hardness.estimate_sgw_resta_qp_gap` (with ``self_consistent=True``):

- Runs an iterative Dyson loop up to ``max_iter=25`` with a default energy convergence tolerance of ``tol=1e-4`` eV.
- In each cycle :math:`n`, updates the Penn core dielectric permittivity :math:`\epsilon_{\mathrm{eff}}(E_g^{(n)})`, recomputes the screening contrast :math:`\Delta W^{(n)}`, and updates dynamic renormalization factors :math:`Z_H^{(n)}, Z_L^{(n)}`.
- Stabilizes iteration dynamics via linear damping: :math:`E_g^{(n+1)} = (1 - \alpha) E_g^{(n)} + \alpha E_g^{\mathrm{target}}` with mixing factor :math:`\alpha = 0.5`.
- Successfully converges in 3–5 iterations with negligible computational overhead (:math:`< 0.5\text{ s}`).

Quasiparticle Self-Consistent GW (``qsgw-dim``, ``qsgw-resta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When orbital wavefunctions undergo significant polarization (e.g. in core/shell nanocrystals, type-II heterostructures, or surface defect states), the single-particle wavefunctions must be relaxed.

``QDEX`` implements full AO-basis **Quasiparticle Self-Consistent GW (:math:`qsGW`)**:

1. Constructs the effective non-local Hamiltonian in the full :math:`N_{\mathrm{AO}} \times N_{\mathrm{AO}}` atomic orbital basis:

   .. math::

      \mathbf{H}_{\mathrm{eff}} = \mathbf{H}_{\mathrm{DFT}} + \Delta \mathbf{H}_{\mathrm{bulk}} + Z \left( \boldsymbol{\Sigma}^{\mathrm{SEX}}[\mathbf{P}] + \boldsymbol{\Sigma}^{\mathrm{COH}} \right)

   where:

   * :math:`\mathbf{P} = \mathbf{C}_{\mathrm{occ}} \mathbf{C}_{\mathrm{occ}}^T` is the single-particle density matrix.
   * :math:`\Sigma_{\mu \nu}^{\mathrm{SEX}} = -P_{\mu \nu} \, \Delta W_{AB}` for :math:`\mu \in A, \nu \in B`.
   * :math:`\Sigma_{\mu \nu}^{\mathrm{COH}} = \frac{1}{2} \, \Delta W_{AA} \, S_{\mu \nu}`.

2. Solves the generalized Hermitian eigenvalue problem:

   .. math::

      \mathbf{H}_{\mathrm{eff}} \, \mathbf{C}_{\mathrm{QP}} = \mathbf{S} \, \mathbf{C}_{\mathrm{QP}} \, \mathbf{E}_{\mathrm{QP}}.

3. Evaluates orbital fidelity to monitor wavefunction reconstruction:

   .. math::

      \mathcal{F}_p = |\langle \psi_p^{(0)} | \psi_p^{(n)} \rangle| = | (\mathbf{C}_p^{(0)})^\dagger \mathbf{S} \, \mathbf{C}_p^{(n)} |.

Implementation in QDEX (``qsgw-dim``, ``qsgw-resta``, ``--update_orbitals``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full AO-basis orbital relaxation is implemented in :func:`qdex.hardness.estimate_qsgw_dim_qp_gap` and :func:`qdex.hardness.estimate_qsgw_resta_qp_gap`:

1. **Löwdin Symmetric Orthogonalization**: Constructs the transformation :math:`\mathbf{S}^{-1/2} = \mathbf{U}_S \boldsymbol{\Lambda}_S^{-1/2} \mathbf{U}_S^T` from the overlap matrix :math:`\mathbf{S}` and transforms initial DFT orbitals into the orthogonal Löwdin basis: :math:`\mathbf{C}_{\mathrm{Low}}^{(0)} = \mathbf{S}^{1/2} \mathbf{C}_{\mathrm{DFT}}`.
2. **Hamiltonian in Löwdin Basis**: Evaluates :math:`\mathbf{H}_{\mathrm{DFT}}^{\mathrm{Low}} = \mathbf{C}_{\mathrm{Low}}^{(0)} \operatorname{diag}(\boldsymbol{\varepsilon}^{\mathrm{DFT}}) (\mathbf{C}_{\mathrm{Low}}^{(0)})^T`.
3. **AO Screening Contrast**: Expands atom-by-atom screening contrast :math:`\Delta W_{AB}` to atomic orbital block pairs :math:`\Delta W_{\mu \nu} = \Delta W_{AB}` for :math:`\mu \in A, \nu \in B`.
4. **Iterative Self-Consistent Loop** (up to ``max_iter=30``):
   - Constructs single-particle density matrix :math:`\mathbf{P}^{\mathrm{Low}} = 2 \mathbf{C}_{\mathrm{occ}}^{\mathrm{Low}} (\mathbf{C}_{\mathrm{occ}}^{\mathrm{Low}})^T` and projector :math:`\mathbf{Q}^{\mathrm{Low}} = \mathbf{I} - 0.5 \mathbf{P}^{\mathrm{Low}}`.
   - Computes frontier charges :math:`\mathbf{q}_H, \mathbf{q}_L` and dynamic asymmetry fractions :math:`f_H^{\mathrm{micro}}, f_L^{\mathrm{micro}}` (Approach B).
   - Assembles bulk reference operator: :math:`\mathbf{H}_{\mathrm{bulk}}^{\mathrm{Low}} = -f_H^{\mathrm{micro}} \Delta_{\mathrm{bulk}} (0.5 \mathbf{P}^{\mathrm{Low}}) + f_L^{\mathrm{micro}} \Delta_{\mathrm{bulk}} \mathbf{Q}^{\mathrm{Low}}`.
   - Builds non-local self-energies: :math:`\boldsymbol{\Sigma}^{\mathrm{SEX}} = -0.5 \mathbf{P}^{\mathrm{Low}} \odot \Delta \mathbf{W}_{\mathrm{AO}}` and :math:`\boldsymbol{\Sigma}^{\mathrm{COH}} = 0.5 \operatorname{diag}(\Delta \mathbf{W}_{\mathrm{AO}})`.
   - Updates target effective Hamiltonian: :math:`\mathbf{H}_{\mathrm{eff}}^{\mathrm{target}} = \mathbf{H}_{\mathrm{DFT}}^{\mathrm{Low}} + \mathbf{H}_{\mathrm{bulk}}^{\mathrm{Low}} + Z (\boldsymbol{\Sigma}^{\mathrm{SEX}} + \boldsymbol{\Sigma}^{\mathrm{COH}})`.
   - Mixes with previous Hamiltonian using linear damping :math:`\mathbf{H}_{\mathrm{eff}} = (1 - \alpha) \mathbf{H}_{\mathrm{eff}}^{\mathrm{prev}} + \alpha \mathbf{H}_{\mathrm{eff}}^{\mathrm{target}}` (:math:`\alpha = 0.5`).
   - Diagonalizes :math:`\mathbf{H}_{\mathrm{eff}} \mathbf{C}^{\mathrm{Low}} = \mathbf{C}^{\mathrm{Low}} \boldsymbol{\varepsilon}^{\mathrm{QP}}` via LAPACK/NumPy ``eigh``.
   - Aligns eigenvector phases: :math:`C_{\mu p} \leftarrow C_{\mu p} \times \operatorname{sign}(\langle \psi_p^{\mathrm{prev}} | \psi_p \rangle)` to eliminate arbitrary sign flips.
   - Monitors wavefunction reconstruction fidelity: :math:`\mathcal{F}_p = |(\mathbf{C}_p^{(0)})^\dagger \mathbf{C}_p^{(n)}|^2`.
5. **Back-Transformation to AO Basis**: Converged eigenvectors are transformed back to the non-orthogonal AO basis: :math:`\mathbf{C}_{\mathrm{QP}} = \mathbf{S}^{-1/2} \mathbf{C}^{\mathrm{Low}}`.
6. **Integration with BSE**: In :mod:`qdex.cli`, the relaxed wavefunctions :math:`\mathbf{C}_{\mathrm{QP}}` and quasiparticle eigenvalues :math:`\boldsymbol{\varepsilon}^{\mathrm{QP}}` directly replace the DFT starting point for subsequent Bethe-Salpeter excited-state calculations!

*CLI & YAML Invocation*:

.. code-block:: bash

   # Full qsGW orbital update with DIM polarizable dipoles
   qdex --mos ground_state.mos --material CDSE --qp_gap qsgw-dim --update_orbitals --dynamic_z --eps-out 2.40

---

9. Environmental Dielectric Polarization
----------------------------------------

In practical applications, colloidal quantum dots are dispersed in liquid solvents (toluene :math:`\epsilon_{\mathrm{out}} = 2.38`, hexane :math:`\epsilon_{\mathrm{out}} = 1.88`) or embedded in polymer or dielectric matrices.

Classical Image Charge Solvation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a carrier charge :math:`q` is placed inside a sphere of core permittivity :math:`\epsilon_\infty` surrounded by solvent permittivity :math:`\epsilon_{\mathrm{out}}`, image charges are induced at the boundary:

.. math::

   \Delta_{\mathrm{pol}}(R, \epsilon_{\mathrm{out}}) = -\frac{e^2}{R} \left( \frac{\epsilon_\infty - \epsilon_{\mathrm{out}}}{\epsilon_\infty + \epsilon_{\mathrm{out}}} \right).

* **Vacuum (:math:`\epsilon_{\mathrm{out}} = 1.0`)**: Image potential is repulsive; gap opens due to lack of external dielectric screening.
* **Dielectric Matching (:math:`\epsilon_{\mathrm{out}} = \epsilon_\infty`)**: Boundary polarization vanishes identically.
* **High-Dielectric Matrix (:math:`\epsilon_{\mathrm{out}} > \epsilon_\infty`)**: External screening reduces carrier charging energies, compressing the fundamental gap.

---

10. Rigid Scissor Operator & Frontier Splitting Strategies (Approach A vs. Approach B)
--------------------------------------------------------------------------------------

To predict absolute valence and conduction band edge alignments, the total quasiparticle gap scissor shift :math:`\Delta_{\mathrm{total}}` must be partitioned between occupied and virtual manifolds:

.. math::

   \varepsilon_i^{\mathrm{QP}} = \varepsilon_i^{\mathrm{DFT}} - f_{\mathrm{homo}} \, \Delta_{\mathrm{total}} \quad (i \in \mathrm{occ}), \qquad
   \varepsilon_a^{\mathrm{QP}} = \varepsilon_a^{\mathrm{DFT}} + f_{\mathrm{lumo}} \, \Delta_{\mathrm{total}} \quad (a \in \mathrm{virt})

where :math:`f_{\mathrm{homo}} + f_{\mathrm{lumo}} = 1.0`. By referencing eigenvalues to the vacuum zero, ``QDEX`` yields the true absolute **Ionization Potential (IP)** and **Electron Affinity (EA)**:

.. math::

   \mathrm{IP} = -\varepsilon_{\mathrm{HOMO}}^{\mathrm{QP}}, \qquad \mathrm{EA} = -\varepsilon_{\mathrm{LUMO}}^{\mathrm{QP}}.

``QDEX`` provides two distinct physical strategies for determining :math:`f_{\mathrm{homo}}` and :math:`f_{\mathrm{lumo}}`:

Approach A: Database Monomer Anchor Frontier Splitting (Used by ``sgw-anchor``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the two-anchor model (``sgw-anchor``), individual self-energy corrections for the HOMO and LUMO levels are extracted from high-level vacuum GW calculations on the monomer anchor cluster stored in ``MATERIAL_DB``:

.. math::

   \delta_h = \varepsilon_{\mathrm{HOMO}}^{\mathrm{GW, mono}} - \varepsilon_{\mathrm{HOMO}}^{\mathrm{PBE, mono}} \le 0, \quad
   \delta_l = \varepsilon_{\mathrm{LUMO}}^{\mathrm{GW, mono}} - \varepsilon_{\mathrm{LUMO}}^{\mathrm{PBE, mono}} \ge 0.

The asymmetry fractions are:

.. math::

   f_{\mathrm{homo}}^{\mathrm{anchor}} = \frac{-\delta_h}{\delta_l - \delta_h}, \quad
   f_{\mathrm{lumo}}^{\mathrm{anchor}} = \frac{\delta_l}{\delta_l - \delta_h}.

*Applicability*: Fast and reliable for pristine, stoichiometric nanocrystals with tabulated material data.

Approach B: Microscopic Wavefunction Asymmetry (Used by ``sgw-dim``, ``sgw-resta``, ``evgw``, ``qsgw``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the microscopic :math:`\Delta W` models, ``QDEX`` **completely skips the need for monomer anchor clusters** by calculating the self-energy shifts directly on the actual frontier wavefunctions of the quantum dot:

.. math::

   \Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} = \frac{1}{2} Z_H \sum_{A, B=1}^{N_{\mathrm{atoms}}} q_A^{\mathrm{HOMO}} \, \Delta W_{AB} \, q_B^{\mathrm{HOMO}}, \quad
   \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}} = \frac{1}{2} Z_L \sum_{A, B=1}^{N_{\mathrm{atoms}}} q_A^{\mathrm{LUMO}} \, \Delta W_{AB} \, q_B^{\mathrm{LUMO}}.

The state-specific asymmetry fractions are computed as:

.. math::

   f_H^{\mathrm{micro}} = \frac{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}}}{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}, \quad
   f_L^{\mathrm{micro}} = \frac{\Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}.

Because the HOMO (typically anion :math:`p` orbitals) and LUMO (typically cation :math:`s` orbitals) possess distinct spatial delocalizations and chemical hardnesses, :math:`f_H^{\mathrm{micro}}` and :math:`f_L^{\mathrm{micro}}` capture the true physical asymmetry of the dot.

Furthermore, in full AO-basis quasiparticle self-consistent GW (:math:`qsGW`), this exact same ratio partitions the bulk reference Hamiltonian:

.. math::

   \mathbf{H}_{\mathrm{bulk}} = -f_H^{\mathrm{micro}} \Delta_{\mathrm{bulk}} (0.5 \mathbf{P}_{\mathrm{occ}}) + f_L^{\mathrm{micro}} \Delta_{\mathrm{bulk}} \mathbf{Q}_{\mathrm{virt}}.

*Why Approach B is the Recommended Default*:

1. **No Anchor Database Required**: Works for any chemical composition, core/shell hetero-interface, facet termination, or organic ligand shell without requiring pre-computed monomer cluster data.
2. **First-Principles Consistency**: Dynamically evolves as nanocrystal size, shape, and aspect ratio change.
3. **Physical Soundness**: Reflects the actual orbital localization of the frontier states rather than an idealized small-molecule surrogate.

Implementation in QDEX (Frontier Alignment & CLI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frontier level alignment is orchestrated between :mod:`qdex.hardness` and the main command-line entry point :mod:`qdex.cli`:

1. **Provenance Inspection**:
   When any microscopic model (``sgw-dim``, ``sgw-resta``, ``evgw-dim``, ``evgw-resta``, ``qsgw-dim``, ``qsgw-resta``) is selected, :mod:`qdex.hardness` returns ``f_homo_micro`` and ``f_lumo_micro`` in the ``provenance`` dictionary.
   
2. **Automatic Dispatching in CLI**:
   - If ``f_homo_micro`` is present, :mod:`qdex.cli` automatically activates **Approach B** and logs:

     .. code-block:: text

        [Absolute Band Edges (IP & EA - Microscopic Wavefunction Asymmetry)]
          -> Shift Split: HOMO takes 55.2%, LUMO takes 44.8% (SGW_DIM)

   - If ``sgw-anchor`` (or ``gw``) is selected, the CLI falls back to **Approach A**, reading :math:`f_{\mathrm{homo}}, f_{\mathrm{lumo}}` from ``MATERIAL_DB``.

3. **Rigid Scissor Application**:
   Occupied Kohn-Sham orbital energies are lowered by :math:`-f_{\mathrm{homo}} \Delta_{\mathrm{total}}` and virtual energies are raised by :math:`+f_{\mathrm{lumo}} \Delta_{\mathrm{total}}`:

   .. math::

      \varepsilon_i^{\mathrm{QP}} = \varepsilon_i^{\mathrm{DFT}} - f_{\mathrm{homo}} \Delta_{\mathrm{total}} \quad (i \le \mathrm{HOMO}), \qquad
      \varepsilon_a^{\mathrm{QP}} = \varepsilon_a^{\mathrm{DFT}} + f_{\mathrm{lumo}} \Delta_{\mathrm{total}} \quad (a \ge \mathrm{LUMO}).

4. **Absolute IP & EA Output**:
   When ``--qp_energy_reference vacuum`` is set (the default), the CLI reports the absolute Ionization Potential (:math:`\mathrm{IP} = -\varepsilon_{\mathrm{HOMO}}^{\mathrm{QP}}`) and Electron Affinity (:math:`\mathrm{EA} = -\varepsilon_{\mathrm{LUMO}}^{\mathrm{QP}}`), establishing direct contact with ultraviolet photoelectron spectroscopy (UPS) and cyclic voltammetry experiments.

---

11. Material Database Reference
-------------------------------

Every energy in the tables below is spin-free. The two bulk columns are a scalar PBE fundamental gap and a spin-free GW gap. Their difference is :math:`\Delta_{\mathrm{bulk}}`. The cluster columns are the same difference on the finite anchor cluster, plus the separate HOMO and LUMO shifts that fix :math:`f_{\mathrm{homo}}` and :math:`f_{\mathrm{lumo}}`.

Perovskites
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0` (Å)
     - PBE bulk
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - Cs\ :sub:`3`\ Bi\ :sub:`2`\ Br\ :sub:`9`
     - 7.85
     - 3.33
     - 3.65
     - +0.32
     - 3.65
     - 7.35
     - +3.70
     - 0.46
     - 0.54
   * - CsPbCl\ :sub:`3`
     - 7.60
     - 2.45
     - 3.65
     - +1.20
     - 3.24
     - 6.90
     - +3.66
     - 0.50
     - 0.50
   * - CsPbBr\ :sub:`3`
     - 7.91
     - 1.85
     - 2.85
     - +1.00
     - 3.08
     - 6.83
     - +3.75
     - 0.43
     - 0.57
   * - CsPbI\ :sub:`3`
     - 8.37
     - 1.55
     - 2.45
     - +0.90
     - 2.69
     - 6.02
     - +3.33
     - 0.34
     - 0.66
   * - MAPbI\ :sub:`3`
     - 
     - 1.55
     - 2.50
     - +0.95
     - 
     - 
     - 
     - 
     - 
   * - FAPbI\ :sub:`3`
     - 
     - 1.45
     - 2.35
     - +0.90
     - 
     - 
     - 
     - 
     - 

II–VI Semiconductors
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0` (Å)
     - PBE bulk
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - ZnS
     - 4.82
     - 2.11
     - 3.73
     - +1.62
     - 3.45
     - 7.26
     - +3.81
     - 0.43
     - 0.57
   * - ZnSe
     - 4.94
     - 1.25
     - 2.76
     - +1.51
     - 3.29
     - 6.84
     - +3.55
     - 0.40
     - 0.60
   * - ZnTe
     - 5.19
     - 1.18
     - 2.35
     - +1.17
     - 2.98
     - 6.10
     - +3.12
     - 0.39
     - 0.61
   * - CdS
     - 5.15
     - 1.14
     - 2.55
     - +1.41
     - 2.72
     - 6.36
     - +3.65
     - 0.44
     - 0.56
   * - CdSe
     - 5.26
     - 0.64
     - 1.91
     - +1.27
     - 2.63
     - 6.03
     - +3.39
     - 0.41
     - 0.59
   * - CdTe
     - 5.44
     - 0.61
     - 1.62
     - +1.01
     - 2.82
     - 6.08
     - +3.26
     - 0.37
     - 0.63
   * - HgS
     - 5.08
     - 0.00
     - 0.50
     - +0.50
     - 2.23
     - 5.47
     - +3.25
     - 0.38
     - 0.62
   * - HgSe
     - 5.20
     - −0.42
     - 0.12
     - +0.54
     - 2.16
     - 5.16
     - +3.00
     - 0.35
     - 0.65
   * - HgTe
     - 5.41
     - −0.61
     - 0.10
     - +0.71
     - 2.18
     - 4.90
     - +2.72
     - 0.34
     - 0.66

III–V Semiconductors
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0` (Å)
     - PBE bulk
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - AlP
     - 5.18
     - 1.59
     - 2.44
     - +0.85
     - 1.94
     - 5.84
     - +3.90
     - 0.47
     - 0.53
   * - AlAs
     - 5.23
     - 1.39
     - 2.15
     - +0.76
     - 2.10
     - 5.73
     - +3.63
     - 0.45
     - 0.55
   * - AlSb
     - 5.56
     - 1.18
     - 1.83
     - +0.65
     - 1.93
     - 5.17
     - +3.24
     - 0.43
     - 0.57
   * - GaP
     - 5.32
     - 1.61
     - 2.36
     - +0.75
     - 1.98
     - 5.40
     - +3.41
     - 0.48
     - 0.52
   * - GaAs
     - 5.45
     - 0.49
     - 1.42
     - +0.93
     - 1.89
     - 5.03
     - +3.13
     - 0.45
     - 0.55
   * - GaSb
     - 5.72
     - 0.11
     - 0.77
     - +0.66
     - 1.48
     - 4.31
     - +2.83
     - 0.44
     - 0.56
   * - InP
     - 5.63
     - 0.46
     - 1.44
     - +0.98
     - 1.75
     - 4.99
     - +3.24
     - 0.45
     - 0.55
   * - InAs
     - 5.70
     - −0.42
     - 0.35
     - +0.77
     - 1.60
     - 4.61
     - +3.01
     - 0.43
     - 0.57
   * - InSb
     - 5.95
     - −0.61
     - 0.23
     - +0.84
     - 1.56
     - 4.30
     - +2.74
     - 0.42
     - 0.58

---

12. Comprehensive Quasiparticle Options Matrix
----------------------------------------------

The following table summarizes all quasiparticle modes supported by ``QDEX``:

.. list-table::
   :widths: 16 14 12 12 46
   :header-rows: 1

   * - ``qp_gap`` Choice
     - Self-Consistency
     - Orbitals
     - Execution Time
     - Recommended Application
   * - ``sgw-anchor`` (alias ``gw``)
     - One-shot
     - DFT (Fixed)
     - :math:`< 1\text{ ms}`
     - High-throughput screening & multi-thousand frame NAMD trajectories.
   * - ``sgw-dim``
     - One-shot
     - DFT (Fixed)
     - :math:`\sim 0.1\text{ s}`
     - Microscopic atomistic polarizable dipole screening without monomer anchors.
   * - ``sgw-resta``
     - One-shot
     - DFT (Fixed)
     - :math:`\sim 0.05\text{ s}`
     - Resta Thomas-Fermi screening with Penn size-dependent dielectric scaling.
   * - ``evgw-dim``
     - Eigenvalue
     - DFT (Fixed)
     - :math:`\sim 0.5\text{ s}`
     - Removes DFT starting-point eigenvalue bias iteratively via Dyson equation.
   * - ``evgw-resta``
     - Eigenvalue
     - DFT (Fixed)
     - :math:`\sim 0.2\text{ s}`
     - Fast eigenvalue self-consistency with Resta screening profile.
   * - ``qsgw-dim``
     - Full Quasiparticle
     - **Updated (AO)**
     - :math:`\sim 5\text{ s}`
     - Full AO orbital relaxation in polarizable dielectric environment.
   * - ``qsgw-resta``
     - Full Quasiparticle
     - **Updated (AO)**
     - :math:`\sim 3\text{ s}`
     - Full AO orbital relaxation with Resta screening profile.
   * - ``brus``
     - Empirical
     - DFT (Fixed)
     - :math:`< 1\text{ ms}`
     - Standard effective-mass confinement approximation.
   * - ``pbe``
     - None
     - DFT (Fixed)
     - :math:`0\text{ s}`
     - Uncorrected Kohn-Sham DFT eigenvalues.

---

13. Recommended Workflow Presets
--------------------------------

Preset 1: Ultra-Fast NAMD Trajectory Dynamics (``sgw-anchor``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For non-adiabatic dynamics across 1,000–10,000 time steps, ``sgw-anchor`` provides sub-millisecond evaluation while reproducing experimental band gaps and solvatochromic shifts:

.. code-block:: yaml

   physics:
     qp_gap: "sgw-anchor"      # Two-anchor scaled GW (backward-compatible alias: "gw")
     material: "CSPBBR3"
     eps_out: 2.25
     excitation_mode: "diagonal_bse"
     kernel: "resta"
     2e-integrals: "mnok"

Preset 2: First-Principles Nanocrystal Screening (``sgw-dim``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For accurate single-point spectra without relying on calibrated monomer anchors, use atomistic polarizable dipoles with dynamic :math:`Z_p` from the plasmon-pole :math:`f`-sum rule:

.. code-block:: yaml

   physics:
     qp_gap: "sgw-dim"         # Microscopic Delta-W with Atomistic Polarizable Dipoles
     dynamic_z: true           # State-dependent Z_p via Plasmon-Pole f-sum rule
     material: "CDSE"
     eps_out: 2.40
     excitation_mode: "bse"
     kernel: "dim"
     2e-integrals: "mnok"

Preset 3: Eigenvalue Self-Consistency (``evgw-dim``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To remove starting-point DFT eigenvalue bias without the computational overhead of full orbital updates:

.. code-block:: yaml

   physics:
     qp_gap: "evgw-dim"        # Eigenvalue self-consistent GW with DIM screening
     material: "CDSE"
     eps_out: 2.40
     excitation_mode: "bse"
     kernel: "resta"
     2e-integrals: "mnok"

Preset 4: Full Quasiparticle Self-Consistency (``qsgw-dim``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For core/shell quantum dots, type-II heterojunctions, or surface-passivated dots where orbital polarization is critical:

.. code-block:: yaml

   physics:
     qp_gap: "qsgw-dim"        # Full AO-basis Quasiparticle Self-Consistent GW
     update_orbitals: true     # Diagonalizes H_eff to relax wavefunctions
     dynamic_z: true
     material: "CDSE"
     eps_out: 2.40
     excitation_mode: "bse"
     kernel: "dim"
     2e-integrals: "mnok"

---

14. Combining Electronic Structure Analysis with QP Shifts
----------------------------------------------------------

``QDEX`` seamlessly integrates QP corrections into ground-state electronic structure analyses:

* **QP-Corrected PDOS**: Convolves the projected density of states on the corrected quasiparticle energy axis :math:`\varepsilon^{\mathrm{QP}}`, opening the gap to experimental values while preserving Mulliken orbital weights.
* **QP-Corrected Fuzzy Bands**: Unfolds nanocrystal orbitals onto bulk :math:`k`-paths using the quasiparticle dispersion:

  .. math::

     P_{\mathbf{k}}(\varepsilon) = \sum_m |\langle e^{i \mathbf{k} \cdot \mathbf{r}} | \phi_m \rangle|^2 \, \delta(\varepsilon - \varepsilon_m^{\mathrm{QP}}).

Controlled via ``--dashboard_energy_mode {dft, qp, both}``.

---

15. CLI Flags & YAML Configuration Reference
--------------------------------------------

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--qp_gap <choice>``
     - ``brus``
     - Quasiparticle model: ``sgw-anchor`` (alias ``gw``), ``sgw-dim``, ``sgw-resta``, ``evgw-dim``, ``evgw-resta``, ``qsgw-dim``, ``qsgw-resta``, ``brus``, ``pbe``, or numeric gap in eV.
   * - ``--dynamic_z``
     - ``False``
     - Compute state-dependent dynamic renormalization :math:`Z_p` from the plasmon-pole :math:`f`-sum rule.
   * - ``--update_orbitals``
     - ``False``
     - Perform full AO-basis Quasiparticle Self-Consistent GW (:math:`qsGW`) orbital update.
   * - ``--material <name>``
     - ``DEFAULT``
     - Material key in ``MATERIAL_DB`` (e.g. ``CSPBBR3``, ``CDSE``, ``INAS``).
   * - ``--eps-out <float>``
     - ``2.0``
     - Optical dielectric constant :math:`\epsilon_{\mathrm{out}}` of surrounding solvent or matrix.
   * - ``--qp-regularization-length <float>``
     - ``1.0``
     - Short-range regularization length :math:`\ell` (in Å) in the anchor confinement formula.
   * - ``--qp-residual-power <float>``
     - ``2.0``
     - Exponent :math:`p` for the quantum confinement power-law decay.
   * - ``--dashboard_energy_mode <choice>``
     - ``dft``
     - Energy axis for Fuzzy Band and PDOS dashboards: ``dft``, ``qp``, or ``both``.
   * - ``--qp_energy_reference <choice>``
     - ``vacuum``
     - Reference zero for QP spectra: ``vacuum`` (absolute IP/EA) or ``fermi`` (:math:`E_F = 0`).

YAML Configuration Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   physics:
     qp_gap: "sgw-anchor"      # "sgw-anchor", "sgw-dim", "evgw-dim", "qsgw-dim", "brus", "pbe"
     dynamic_z: true           # Dynamically compute Z_p via PPM f-sum rule
     update_orbitals: false    # Full AO orbital relaxation (qsGW)
     material: "CSPBBR3"
     eps_out: 2.25
     qp_regularization_length: 1.0
     qp_residual_power: 2.0

   fuzzy:
     run: true
     dashboard_energy_mode: "both"
     qp_energy_reference: "vacuum"
