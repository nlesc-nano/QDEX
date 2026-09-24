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

What the interpolation actually evaluates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The radius :math:`R` is the equivalent-volume radius of the inorganic core, :math:`R = (3V/4\pi)^{1/3}`. The scissor added to a PBE eigenvalue difference is

.. math::

   \Delta_{\mathrm{GW}}(R)
   = \Delta_{\mathrm{bulk}}
   + \frac{\kappa_{\mathrm{out}}}{R+\ell}
   + A\left(\frac{R_0}{R}\right)^{p}

with the bulk piece taken from the two bulk gaps in the table below,

.. math::

   \Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW,bulk}} - E_g^{\mathrm{PBE,bulk}}.

The dielectric prefactor is :math:`\kappa_{\mathrm{out}} = 11.52\,\mathrm{eV\,\AA}\,(1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty)`. In vacuum, :math:`\epsilon_{\mathrm{out}}=1` and :math:`\kappa_{\mathrm{vac}} = 11.52\,\mathrm{eV\,\AA}\,(1 - 1/\epsilon_\infty)`. The defaults are :math:`\ell = 1.0\,\text{Å}` and :math:`p = 2`.

:math:`A` is not a free fit. It is fixed by the finite cluster so that, in vacuum and at :math:`R=R_0`, the whole expression equals the cluster gap opening :math:`\Delta(R_0) = E_g^{\mathrm{GW,cluster}} - E_g^{\mathrm{PBE,cluster}}`:

.. math::

   A = \Delta(R_0) - \Delta_{\mathrm{bulk}} - \frac{\kappa_{\mathrm{vac}}}{R_0+\ell}.

Both table columns are therefore required. :math:`\Delta_{\mathrm{bulk}}` is what survives for a large crystal. :math:`\Delta(R_0)` is what the small cluster actually opened, and :math:`A` carries the difference. At :math:`R_0` the bulk term cancels and the vacuum scissor is exactly the cluster opening. As :math:`R` grows, :math:`(R_0/R)^p` and :math:`1/(R+\ell)` die, and the scissor falls to :math:`\Delta_{\mathrm{bulk}}`. A PBE calculation of the bulk is then shifted onto the spin-free GW gap stored in the table. Materials without a finite anchor (MAPbI\ :sub:`3`, FAPbI\ :sub:`3`) keep only :math:`\Delta_{\mathrm{bulk}}` and the dielectric term.

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

In the implemented model this environment dependence is the term :math:`\kappa_{\mathrm{out}}/(R+\ell)` written above, with :math:`\kappa_{\mathrm{out}} \propto (1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty)`. It vanishes when the surroundings have the same optical dielectric constant as the crystal, it is positive in vacuum, and it changes sign when the matrix screens more strongly than the crystal.

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

7. Where the anchor numbers come from
-------------------------------------

Every energy in the tables below is spin-free. CP2K, which produces the finite clusters, cannot apply spin–orbit coupling with the GTH pseudopotentials used here, so a bulk number that already contains spin–orbit coupling is not the same correction and is not stored. The two bulk columns are a spin-free PBE gap and a spin-free GW gap. Their difference is :math:`\Delta_{\mathrm{bulk}}`. The cluster columns are the same difference on the finite anchor, plus the separate HOMO and LUMO shifts that fix :math:`f_{\mathrm{homo}}` and :math:`f_{\mathrm{lumo}}`.

Finite cluster
~~~~~~~~~~~~~~

The II–VI anchors are not diatomic monomers. They are stoichiometric Wulff clusters :math:`\mathrm{M}_{16}\mathrm{X}_{13}\mathrm{Cl}_6` with :math:`\mathrm{M} = \mathrm{Zn,\,Cd,\,Hg}` and :math:`\mathrm{X} = \mathrm{S,\,Se,\,Te}`, chlorine-terminated, relaxed with spin-free PBE. :math:`R_0` is the convex-hull radius of that relaxed geometry. The III–V and perovskite anchors are the spin-free clusters whose eigenvalues are stored; :math:`R_0` is obtained the same way. MAPbI\ :sub:`3` and FAPbI\ :sub:`3` have bulk gaps only.

On that fixed geometry the electronic protocol is:

1. Spin-free PBE. The HOMO and LUMO are the reference that the scissor corrects.
2. Spin-free PBE0 on the same density path: 75% PBE exchange and 25% exact exchange.
3. Eigenvalue-self-consistent GW on those PBE0 orbitals (``EV_GW_ITER 4`` in CP2K), with the frontier occupied and virtual states corrected. The orbitals stay the PBE0 orbitals. This is not one-shot :math:`G_0W_0` on the PBE density, and it is not quasiparticle-self-consistent GW.

One-shot :math:`G_0W_0` depends on the functional that builds :math:`G_0` and :math:`W_0`. Iterating the eigenvalues removes most of that dependence from the gap. The anchor opening is therefore

.. math::

   \Delta(R_0) = \big(E_{\mathrm{LUMO}}-E_{\mathrm{HOMO}}\big)_{\mathrm{evGW@PBE0}}
   - \big(E_{\mathrm{LUMO}}-E_{\mathrm{HOMO}}\big)_{\mathrm{PBE}}.

That opening is several eV on these clusters because confinement and the self-energy both contribute. It is the number :math:`A` is built from. The separate HOMO and LUMO shifts, not the gap alone, set how the scissor is split between occupied and virtual states.

Bulk crystal
~~~~~~~~~~~~

The bulk PBE column is a scalar-relativistic PBE fundamental gap: highest occupied band to lowest empty band, no spin–orbit coupling. For an ordinary zincblende or perovskite crystal this is a positive number and is taken from scalar PBE literature (and, for the mercury salts, checked by an AMS BAND calculation with scalar ZORA and PBE).

The bulk GW column is the spin-free quasiparticle gap at the same definition. A published GW gap that includes spin–orbit coupling is converted before it is stored, and only when the band order is the normal one. Spin–orbit coupling raises the :math:`\Gamma_8` valence maximum by :math:`\Delta_0/3` above the scalar :math:`\Gamma_{15}` barycenter, so

.. math::

   E_g^{\mathrm{spin\text{-}free}} \approx E_0 + \frac{\Delta_0}{3},

where :math:`E_0` is the spin–orbit gap and :math:`\Delta_0 = E(\Gamma_8)-E(\Gamma_7)` is the valence splitting. The whole splitting :math:`\Delta_0` is not added. This correction is a few hundredths of an eV for ZnS and CdS and about :math:`0.3\,\mathrm{eV}` for CdTe, ZnTe and GaSb. It does not apply to the lead salts, whose gap is at L.

HgTe is the check that this conversion was done on the right states. A scalar ZORA PBE band structure (AMS, DZP, :math:`a = 6.43\,\text{Å}`) has no fundamental gap: the Te :math:`p` states touch at the Fermi level. The Hg :math:`s` state lies :math:`0.620\,\mathrm{eV}` below them, so the signed gap that matches a normal II–VI gap is :math:`E(\Gamma_1)-E(\Gamma_{15}) = -0.620\,\mathrm{eV}`. The table stores :math:`-0.61\,\mathrm{eV}`. The same cell with spin–orbit ZORA gives :math:`\Delta_0 = 0.768\,\mathrm{eV}` and :math:`E_0 = E(\Gamma_6)-E(\Gamma_8) = -0.883\,\mathrm{eV}`. Rebuilding the scalar separation,

.. math::

   E_0 + \frac{\Delta_0}{3} = -0.627\,\mathrm{eV},

agrees with the scalar result to :math:`7\,\mathrm{meV}`. The GW entry :math:`+0.10\,\mathrm{eV}` is the same separation after the quasiparticle correction. It matches the hybrid quasiparticle-self-consistent GW result of Svane et al., Phys. Rev. B **84**, 205205 (2011), :math:`E_0 = -0.18\,\mathrm{eV}` and :math:`\Delta_0 = 0.80\,\mathrm{eV}`, converted by the same formula to :math:`+0.09\,\mathrm{eV}`. One-shot :math:`G_0W_0` from that paper stays near :math:`-0.08\,\mathrm{eV}` and is not the value stored. Full quasiparticle-self-consistent GW un-inverts HgTe and is not used. HgSe and InAs, InSb are stored the same way: a negative PBE number is a spin-free :math:`s`–:math:`p` separation, not a spin–orbit :math:`E_0` copied in unchanged. HgS has a closed scalar PBE fundamental gap, stored as 0.

Perovskites
~~~~~~~~~~~

Gaps are in eV and :math:`R_0` is in Å. :math:`\Delta_{\mathrm{bulk}}` and :math:`\Delta(R_0)` are the two openings the interpolator combines. :math:`f_H` and :math:`f_L` are the occupied and virtual fractions of the cluster opening.

.. list-table::
   :header-rows: 1
   :widths: 14 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0`
     - PBE bulk
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE cluster
     - GW cluster
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

II–VI
~~~~~

The cluster is :math:`\mathrm{M}_{16}\mathrm{X}_{13}\mathrm{Cl}_6`. For HgTe the bulk columns are the signed :math:`s`–:math:`p` separation discussed above, not the fundamental gap of zero.

.. list-table::
   :header-rows: 1
   :widths: 10 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0`
     - PBE bulk
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE cluster
     - GW cluster
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

III–V
~~~~~

AlP, AlAs, AlSb and GaP are indirect. The valence maximum is still at :math:`\Gamma`, so :math:`\Delta_0/3` remains the leading spin–orbit correction to a literature gap. InAs and InSb keep a negative spin-free PBE separation; their GW column is the positive spin-free quasiparticle gap.

.. list-table::
   :header-rows: 1
   :widths: 10 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0`
     - PBE bulk
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE cluster
     - GW cluster
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

Read a row as follows. For CdSe, :math:`\Delta_{\mathrm{bulk}} = 1.27\,\mathrm{eV}` and :math:`\Delta(R_0) = 3.39\,\mathrm{eV}`. A dot the size of the anchor gets a vacuum scissor of :math:`3.39\,\mathrm{eV}`, of which 41% lowers the occupied states and 59% raises the virtual states. A very large CdSe crystal gets :math:`1.27\,\mathrm{eV}`. Everything in between is the formula above. For HgTe the same reading gives :math:`2.72\,\mathrm{eV}` on the Wulff cluster and :math:`0.71\,\mathrm{eV}` in the bulk, which moves the signed :math:`s`–:math:`p` separation from :math:`-0.61\,\mathrm{eV}` to :math:`+0.10\,\mathrm{eV}`.

8. Critical Assessment: Pros and Cons
-------------------------------------

Advantages of the Scaled GW Scissor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Extreme Computational Efficiency**: The scissor shift :math:`\Delta_{\mathrm{GW}}` is evaluated in :math:`< 1\text{ ms}`, making it possible to apply quasiparticle corrections across thousands of MD trajectory frames.
2. **Asymptotic Physical Rigor**: In vacuum the scissor equals the finite-cluster opening at :math:`R_0`. As :math:`R \to \infty` it equals :math:`E_g^{\mathrm{GW,bulk}} - E_g^{\mathrm{PBE,bulk}}`, so a bulk PBE gap is shifted onto the spin-free GW gap in the table.
3. **Dielectric Environment Awareness**: Accurately reproduces solvatochromic and dielectric screening shifts without requiring explicit solvent molecules.

Limitations
~~~~~~~~~~~

1. **Rigid Band Shift**: Assumes that all occupied orbitals experience a uniform shift :math:`-f_{\mathrm{homo}}\Delta_{\mathrm{GW}}` and all virtual orbitals experience :math:`+f_{\mathrm{lumo}}\Delta_{\mathrm{GW}}`, neglecting state-specific orbital self-energy variations deep within the bands.
2. **Database Dependency**: Requires calibrated parameters (:math:`\epsilon_\infty`, :math:`R_0`, :math:`\Delta_{\mathrm{bulk}}`) in ``MATERIAL_DB`` for the target semiconductor.

---

9. Combining Electronic Structure Analysis with QP Shifts
---------------------------------------------------------

``QDEX`` seamlessly integrates QP corrections into ground-state electronic structure analyses:

* **QP-Corrected PDOS**: Convolves the projected density of states on the corrected quasiparticle energy axis :math:`\varepsilon^{\mathrm{QP}}`, opening the gap to experimental values while preserving Mulliken orbital weights.
* **QP-Corrected Fuzzy Bands**: Unfolds nanocrystal orbitals onto bulk :math:`k`-paths using the quasiparticle dispersion:

  .. math::

     P_{\mathbf{k}}(\varepsilon) = \sum_m |\langle e^{i \mathbf{k} \cdot \mathbf{r}} | \phi_m \rangle|^2 \, \delta(\varepsilon - \varepsilon_m^{\mathrm{QP}})

Controlled via ``--dashboard_energy_mode {dft, qp, both}``.

---

10. CLI Flags & YAML Configuration Reference
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
