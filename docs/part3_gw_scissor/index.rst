Part 3: Quasiparticle Corrections & The Scaled GW Model
========================================================

Standard semi-local Kohn-Sham Density Functional Theory (DFT) using functionals like PBE severely underestimates the fundamental band gap of semiconductor nanostructures. For instance, PBE predicts a band gap of :math:`\approx 1.5\text{ eV}` for bulk :math:`\text{CsPbBr}_3`, whereas the experimental quasiparticle gap is :math:`\approx 2.35\text{ eV}`.

In ``QDEX``, single-particle excitation energies are corrected via an analytical, physically grounded **Scaled GW Quasiparticle Model** that incorporates quantum confinement, dielectric solvation screening, and absolute energy level alignment.

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

2. The DFT Band Gap Problem
---------------------------

The failure of Kohn-Sham DFT with semi-local functionals (LDA, PBE) to predict fundamental band gaps arises from two fundamental physical sources:

1. **Missing Derivative Discontinuity (:math:`\Delta_{xc}`)**:
   The true fundamental band gap is:

   .. math::

      E_g = \mathrm{IP} - \mathrm{EA} = E_g^{\mathrm{KS}} + \Delta_{xc}

   where :math:`\Delta_{xc} = \left. \frac{\delta E_{xc}}{\delta \rho} \right|_{N+\delta} - \left. \frac{\delta E_{xc}}{\delta \rho} \right|_{N-\delta}` is the integer discontinuity of the exchange-correlation functional. In standard LDA/GGA functionals, :math:`E_{xc}[\rho]` is a continuous function of density, identically yielding :math:`\Delta_{xc} = 0`. Consequently, :math:`E_g^{\mathrm{PBE}}` underestimates the true quasiparticle gap by 30–50%.

2. **Self-Interaction Error (SIE)**:
   In semi-local functionals, an electron spurious interacts with its own charge density through the Hartree term. This unphysical electrostatic repulsion artificially destabilizes occupied states and over-delocalizes frontier wavefunctions.

---

3. The Nanocrystal Scaling Bottleneck
-------------------------------------

Although full :math:`G_0W_0` calculations resolve the band gap problem for small molecules or bulk crystals, computing :math:`G_0W_0` explicitly for colloidal quantum dots is computationally prohibitive:

* **Computational Complexity**: Standard :math:`G_0W_0` scales as :math:`O(N^4)` with real-space basis sets and :math:`O(N^5)` with plane-wave expansions.
* **Memory Footprint**: Calculating the polarizability matrix :math:`\chi_0(\mathbf{r}, \mathbf{r}'; \omega)` requires summing over thousands of unoccupied conduction states, demanding terabytes of RAM for clusters with :math:`> 500` atoms.
* **Nanocrystal Realities**: Chemically realistic colloidal quantum dots comprise 1,000 to 10,000 atoms, including surface passivation ligands, rendering direct *ab initio* :math:`G_0W_0` impossible for high-throughput screening or molecular dynamics trajectories.

---

4. Rationale for the Two-Anchor Model
-------------------------------------

To overcome this bottleneck while maintaining quantitative accuracy, ``QDEX`` introduces the **Two-Anchor Scaled GW Model**. 

Instead of guessing empirical parameters or running prohibitively expensive calculations for every nanocrystal size, the model anchors the quasiparticle gap between two rigorously computed physical limits:

.. list-table::
   :widths: 20 25 55
   :header-rows: 1

   * - Anchor Limit
     - Physical System
     - Theoretical Characterization
   * - **Anchor 1: Smallest Vacuum Anchor (:math:`R_0`)**
     - Monomer or smallest Wulff-constructed cluster
     - Fully relaxed, stoichiometric cluster computed in vacuum with high-level hybrid DFT (:math:`\text{PBE0}`) and full :math:`G_0W_0`. Captures maximum quantum confinement and atomic-scale self-energy shifts.
   * - **Anchor 2: Bulk Limit (:math:`R \to \infty`)**
     - Periodic crystal
     - High-accuracy bulk :math:`G_0W_0` quasiparticle gap (:math:`E_g^{\mathrm{GW, bulk}}`), calibrated against experimental angle-resolved photoemission spectroscopy (ARPES).

Quantum Confinement Scaling Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For an intermediate quantum dot of effective radius :math:`R = \left( \frac{3 V_{\mathrm{cluster}}}{4\pi} \right)^{1/3}`, the quantum confinement contribution to the GW scissor shift interpolates smoothly between the monomer anchor and the bulk limit:

.. math::

   \Delta_{\mathrm{conf}}(R) = \frac{A}{(R + \ell)^p}

where:
* :math:`A` is the confinement amplitude calibrated to match the vacuum anchor at :math:`R = R_0`.
* :math:`\ell` is a short-range regularization length (default 1.0 Å) that prevents unphysical divergences as :math:`R \to 0`.
* :math:`p` is the confinement power (typically :math:`p \approx 1.5 - 2.0`, reflecting the effective mass kinetic energy scaling).

---

5. Environmental Dielectric Polarization
----------------------------------------

In practical applications, semiconductor quantum dots are dispersed in liquid solvents (e.g., toluene :math:`\epsilon_{\mathrm{out}} = 2.38`, hexane :math:`\epsilon_{\mathrm{out}} = 1.88`) or embedded in polymer or oxide matrices.

Classical Image Charge Formalism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When an electron or hole is added to a nanocrystal of core optical dielectric constant :math:`\epsilon_\infty` embedded in a surrounding dielectric medium of permittivity :math:`\epsilon_{\mathrm{out}}`, image charges are induced at the spherical dielectric boundary:

.. math::

   \Delta_{\mathrm{pol}}(R, \epsilon_{\mathrm{out}}) = -\frac{e^2}{R} \left( \frac{\epsilon_\infty - \epsilon_{\mathrm{out}}}{\epsilon_\infty + \epsilon_{\mathrm{out}}} \right)

Physical Behavior:
* **Nanocrystal in Vacuum (:math:`\epsilon_{\mathrm{out}} = 1.0`)**: Since :math:`\epsilon_\infty \gg 1`, :math:`\Delta_{\mathrm{pol}} < 0`, reflecting the absence of external screening (large charging energy).
* **Dielectric Matching (:math:`\epsilon_{\mathrm{out}} = \epsilon_\infty`)**: The boundary image charge vanishes: :math:`\Delta_{\mathrm{pol}} = 0`.
* **High-Dielectric Matrix (:math:`\epsilon_{\mathrm{out}} > \epsilon_\infty`)**: The surrounding medium screens the carrier charge more effectively than the internal lattice, lowering the fundamental gap.

The total scaled GW scissor shift is therefore:

.. math::

   \Delta_{\mathrm{GW}}(R, \epsilon_{\mathrm{out}}) = \Delta_{\mathrm{bulk}} + \frac{A}{(R + \ell)^p} - \frac{e^2}{R} \left( \frac{\epsilon_\infty - \epsilon_{\mathrm{out}}}{\epsilon_\infty + \epsilon_{\mathrm{out}}} \right)

where :math:`\Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW, bulk}} - E_g^{\mathrm{PBE, bulk}}`.

---

6. Rigid Scissor Operator & Frontier Splitting
----------------------------------------------

In standard electronic structure calculations, a constant scissor operator simply shifts the entire conduction band relative to the valence band. However, to correctly predict absolute energy levels, the shift must be partitioned between occupied and virtual manifolds.

Frontier Splitting Fractions (:math:`f_{\mathrm{homo}}, f_{\mathrm{lumo}}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From the high-level vacuum GW calculation on the monomer anchor, ``QDEX`` extracts the individual self-energy corrections for the HOMO and LUMO levels:

.. math::

   \delta_h = \varepsilon_{\mathrm{HOMO}}^{\mathrm{GW, mono}} - \varepsilon_{\mathrm{HOMO}}^{\mathrm{PBE, mono}} \le 0

.. math::

   \delta_l = \varepsilon_{\mathrm{LUMO}}^{\mathrm{GW, mono}} - \varepsilon_{\mathrm{LUMO}}^{\mathrm{PBE, mono}} \ge 0

The asymmetry fractions are defined by the ratio of individual level shifts to the total gap opening:

.. math::

   f_{\mathrm{homo}} = \frac{-\delta_h}{\delta_l - \delta_h}, \quad
   f_{\mathrm{lumo}} = \frac{\delta_l}{\delta_l - \delta_h}

subject to the conservation condition:

.. math::

   f_{\mathrm{homo}} + f_{\mathrm{lumo}} = 1.0

In lead halide perovskites (:math:`\text{CsPbBr}_3`), the occupied valence band (Pb :math:`6s` / Br :math:`4p`) and unoccupied conduction band (Pb :math:`6p`) respond asymmetrically to correlation, resulting in :math:`f_{\mathrm{homo}} \approx 0.43` and :math:`f_{\mathrm{lumo}} \approx 0.57`.

Quasiparticle Energy Shifts
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The quasiparticle energies are computed by applying the partitioned scissor operator:

.. math::

   \varepsilon_i^{\mathrm{QP}} = \varepsilon_i^{\mathrm{DFT}} - f_{\mathrm{homo}} \, \Delta_{\mathrm{GW}} \quad (i \in \mathrm{occ})

.. math::

   \varepsilon_a^{\mathrm{QP}} = \varepsilon_a^{\mathrm{DFT}} + f_{\mathrm{lumo}} \, \Delta_{\mathrm{GW}} \quad (a \in \mathrm{virt})

Dynamic Prediction of Absolute IP & EA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In standard CP2K cluster calculations, the raw eigenvalues float with respect to an arbitrary electrostatic reference potential. By anchoring the levels to the vacuum-calibrated monomer data stored in ``MATERIAL_DB``, ``QDEX`` reconstructs the true absolute **Ionization Potential (IP)** and **Electron Affinity (EA)**:

.. math::

   \mathrm{IP} = -\varepsilon_{\mathrm{HOMO}}^{\mathrm{QP}}, \quad
   \mathrm{EA} = -\varepsilon_{\mathrm{LUMO}}^{\mathrm{QP}}

This dynamic reconstruction bypasses the arbitrary floating vacuum level of standard periodic or supercell DFT codes.

---

7. Critical Assessment: Pros and Cons
-------------------------------------

Advantages of the Scaled GW Scissor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Extreme Computational Efficiency**: The scissor shift :math:`\Delta_{\mathrm{GW}}` is evaluated in :math:`< 1\text{ ms}`, making it possible to apply quasiparticle corrections across thousands of MD trajectory frames.
2. **Asymptotic Physical Rigor**: It guarantees exact behavior at both extremes: recovering the full GW monomer gap at :math:`R_0` and the experimental macroscopic gap as :math:`R \to \infty`.
3. **Dielectric Environment Awareness**: Accurately reproduces solvatochromic and dielectric screening shifts without requiring explicit solvent molecules.

Limitations
~~~~~~~~~~~

1. **Rigid Band Shift**: Assumes that all occupied orbitals experience a uniform shift :math:`-f_{\mathrm{homo}}\Delta_{\mathrm{GW}}` and all virtual orbitals experience :math:`+f_{\mathrm{lumo}}\Delta_{\mathrm{GW}}`, neglecting state-specific orbital self-energy variations deep within the bands.
2. **Database Dependency**: Requires calibrated parameters (:math:`\epsilon_\infty`, :math:`R_0`, :math:`\Delta_{\mathrm{bulk}}`) in ``MATERIAL_DB`` for the target semiconductor.

---

8. Combining Electronic Structure Analysis with QP Shifts
---------------------------------------------------------

``QDEX`` seamlessly integrates QP corrections into ground-state electronic structure analyses:

* **QP-Corrected PDOS**: Convolves the projected density of states on the corrected quasiparticle energy axis :math:`\varepsilon^{\mathrm{QP}}`, opening the gap to experimental values while preserving Mulliken orbital weights.
* **QP-Corrected Fuzzy Bands**: Unfolds nanocrystal orbitals onto bulk :math:`k`-paths using the quasiparticle dispersion:

  .. math::

     P_{\mathbf{k}}(\varepsilon) = \sum_m |\langle e^{i \mathbf{k} \cdot \mathbf{r}} | \phi_m \rangle|^2 \, \delta(\varepsilon - \varepsilon_m^{\mathrm{QP}})

Controlled via ``--dashboard_energy_mode {dft, qp, both}``.

---

9. CLI Flags & YAML Configuration Reference
-------------------------------------------

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--qp_gap <choice>``
     - ``gw``
     - Quasiparticle model: ``gw`` (Scaled GW model), ``brus`` (Brus effective mass model), ``pbe`` (uncorrected DFT), or explicit numeric gap in eV.
   * - ``--material <name>``
     - ``DEFAULT``
     - Material key in ``MATERIAL_DB`` (e.g. ``CSPBBR3``, ``CSPBI3``, ``INAS``, ``CDSE``).
   * - ``--eps-out <float>``
     - ``2.0``
     - Optical dielectric constant :math:`\epsilon_{\mathrm{out}}` of the surrounding solvent or matrix.
   * - ``--qp-regularization-length <float>``
     - ``1.0``
     - Short-range regularization length :math:`\ell` (in Angstroms) in the confinement formula.
   * - ``--qp-residual-power <float>``
     - ``2.0``
     - Exponent :math:`p` for the quantum confinement power-law decay.
   * - ``--dashboard_energy_mode <choice>``
     - ``dft``
     - Energy axis for Fuzzy Band and PDOS dashboards: ``dft``, ``qp``, or ``both``.
   * - ``--qp_energy_reference <choice>``
     - ``vacuum``
     - Reference zero for QP spectra: ``vacuum`` (absolute IP/EA) or ``fermi`` (:math:`E_F = 0`).

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   physics:
     qp_gap: "gw"
     material: "CSPBBR3"
     eps_out: 2.25
     qp_regularization_length: 1.0
     qp_residual_power: 2.0

   fuzzy:
     run: true
     dashboard_energy_mode: "both"
     qp_energy_reference: "vacuum"
