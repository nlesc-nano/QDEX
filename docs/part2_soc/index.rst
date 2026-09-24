Part 2: Relativistic Spin-Orbit Coupling (SOC)
==============================================

In materials containing heavy chemical elements—such as lead halide perovskites (:math:`\text{CsPb}X_3`), bismuth double perovskites (:math:`\text{Cs}_2\text{AgBiBr}_6`), or heavy chalcogenide quantum dots (:math:`\text{InAs}`, :math:`\text{PbS}`, :math:`\text{CdSe}`)—non-relativistic Schrödinger quantum mechanics breaks down. 

Relativistic **Spin-Orbit Coupling (SOC)** fundamentally alters the electronic band structure, lifting orbital degeneracies, reshaping the density of states, and opening or contracting the optical band gap.

---

1. Physical Origin of Spin-Orbit Coupling
-----------------------------------------

The spin-orbit interaction arises naturally from the relativistic **Dirac equation** for a spin-1/2 fermion in an electrostatic potential :math:`V(\mathbf{r})`. Under the low-velocity (Pauli) reduction to order :math:`(v/c)^2`, an electron moving with velocity :math:`\mathbf{v}` experiences the electric field of the nucleus :math:`\mathbf{E} = -\boldsymbol{\nabla} V(\mathbf{r})` transformed in its rest frame into an effective magnetic field :math:`\mathbf{B}_{\mathrm{eff}} = -\frac{1}{c} \mathbf{v} \times \mathbf{E}`.

The Zeeman coupling of the electron's intrinsic magnetic dipole moment :math:`\boldsymbol{\mu}_s = -g_s \frac{e}{2m_e} \mathbf{S}` to this relativistic field yields the spin-orbit Hamiltonian:

.. math::

   \hat{H}_{\mathrm{SO}} = -\boldsymbol{\mu}_s \cdot \mathbf{B}_{\mathrm{eff}} = \frac{\hbar}{4 m_e^2 c^2} \left( \boldsymbol{\nabla} V \times \mathbf{p} \right) \cdot \boldsymbol{\sigma}

For a central spherically symmetric potential :math:`V(r)` where :math:`\boldsymbol{\nabla} V(r) = \frac{\mathbf{r}}{r} \frac{dV}{dr}`, recognizing the orbital angular momentum operator :math:`\mathbf{L} = \mathbf{r} \times \mathbf{p}` and the spin operator :math:`\mathbf{S} = \frac{\hbar}{2} \boldsymbol{\sigma}` yields:

.. math::

   \hat{H}_{\mathrm{SO}} = \frac{1}{2 m_e^2 c^2} \frac{1}{r} \frac{dV(r)}{dr} \, \mathbf{L} \cdot \mathbf{S} = \xi(r) \, \mathbf{L} \cdot \mathbf{S}

The :math:`Z^4` Scaling Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Near the atomic core, the Coulomb potential behaves as :math:`V(r) \approx -Z e^2 / r`, leading to :math:`\frac{1}{r} \frac{dV}{dr} \approx Z e^2 / r^3`. Evaluating the expectation value :math:`\langle r^{-3} \rangle` with hydrogenic wavefunctions produces a strong dependence on the atomic number :math:`Z`:

.. math::

   \langle \xi(r) \rangle \propto Z^4

This steep scaling explains why SOC is negligible for carbon (:math:`Z = 6`) or oxygen (:math:`Z = 8`), moderate for sulfur (:math:`Z = 16`), but overwhelmingly dominant in heavy elements such as iodine (:math:`Z = 53`), cesium (:math:`Z = 55`), lead (:math:`Z = 82`), and bismuth (:math:`Z = 83`).

---

2. Impact of SOC on Semiconductor Nanocrystals
----------------------------------------------

In lead halide perovskites, the frontier conduction band arises from the hybridization of empty lead :math:`6p` orbitals with halogen :math:`np` states. 

Band Inversion & Giant Gap Contraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the absence of SOC, the :math:`p`-like conduction band is triply degenerate (six-fold degenerate including spin). The operator :math:`\mathbf{L} \cdot \mathbf{S}` acts on total angular momentum :math:`\mathbf{J} = \mathbf{L} + \mathbf{S}`:

.. math::

   \mathbf{L} \cdot \mathbf{S} = \frac{1}{2} \left( \mathbf{J}^2 - \mathbf{L}^2 - \mathbf{S}^2 \right) = \frac{\hbar^2}{2} \left[ j(j+1) - l(l+1) - s(s+1) \right]

For :math:`p`-orbitals (:math:`l = 1, s = 1/2`):
* **Lower :math:`j = 1/2` doublet**: Eigenvalue :math:`-\hbar^2`. The conduction band minimum shifts downward in energy by :math:`\approx 0.65\text{ eV}` in :math:`\text{CsPbBr}_3`.
* **Upper :math:`j = 3/2` quartet**: Eigenvalue :math:`+\frac{1}{2}\hbar^2`.

As a direct consequence, standard non-relativistic DFT severely misidentifies the nature of the conduction band minimum. When SOC is activated:
1. The fundamental band gap contracts dramatically.
2. The effective mass of conduction band electrons decreases.
3. Carrier cooling dynamics within the conduction band manifold accelerate significantly due to dense non-adiabatic couplings between spinor levels.

Rashba-Dresselhaus Splitting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the surface of a quantum dot or in the presence of asymmetric ligand termination, the local inversion symmetry is broken. The resulting macroscopic gradient :math:`\boldsymbol{\nabla} V \neq \mathbf{0}` couples with electron momentum via the Rashba mechanism:

.. math::

   \hat{H}_{\mathrm{Rashba}} = \alpha_R \left( \boldsymbol{\sigma} \times \mathbf{k} \right) \cdot \hat{\mathbf{z}}

This splits spin-degenerate states in :math:`k`-space, protecting carriers against non-radiative electron-hole recombination and dramatically prolonging photoluminescence lifetimes.

---

3. Survey of Relativistic Methods & Rationale for GTH Pseudopotentials
----------------------------------------------------------------------

Multiple methodologies exist in modern quantum chemistry and solid-state physics for treating spin-orbit coupling:

.. list-table::
   :widths: 20 25 55
   :header-rows: 1

   * - Approach
     - Scaling / Complexity
     - Physical Characteristics & Limitations
   * - **4-Component Dirac-Kohn-Sham (DKS)**
     - Extremely High (:math:`O(N^4) - O(N^5)`)
     - Solves the fully coupled Dirac 4-spinor problem (large and small components). Computationally intractable for nanoclusters containing thousands of atoms.
   * - **2-Component ZORA / DKH**
     - High (:math:`O(N^3) - O(N^4)`)
     - Eliminates small components through Foldy-Wouthuysen or Zero-Order Regular Approximation. Requires all-electron contracted basis sets with severe picture-change errors and high memory overhead.
   * - **Empirical Tight-Binding (TB)**
     - Low (:math:`O(N^2)`)
     - Introduces empirical atom-centered atomic spin-orbit constants :math:`\lambda_{\mathrm{SO}}`. Lacks *ab initio* wavefunctions and orbital overlap consistency.
   * - **Separable GTH Pseudopotentials (QDEX)**
     - Highly Scalable & Fast (:math:`O(N_{\mathrm{act}}^3)`)
     - Fully *ab initio* relativistic core representation using Goedecker-Teter-Hutter (GTH) separable angular momentum projectors. Perfectly compatible with CP2K GPW calculations.

Why QDEX Chooses GTH SOC Pseudopotentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``QDEX`` is tailored for post-processing CP2K calculations. In CP2K Quickstep, the core electrons are integrated out using norm-conserving GTH pseudopotentials. By evaluating SOC within the exact same **separable GTH projector framework**:
1. **Rigorous Consistency**: The relativistic core potentials match the underlying Kohn-Sham ground state without parameter re-fitting.
2. **Semi-Local Projector Form**: The SOC Hamiltonian is completely factorized into separable inner products with atom-centered Gaussian projector functions.
3. **Analytic C++ Evaluation**: Matrix elements of the Gaussian projectors are evaluated analytically via Libint2, eliminating numerical integration grids.

---

4. Mathematical Formulation: Separable GTH Pseudopotential
----------------------------------------------------------

In the GTH relativistic pseudopotential formalism, the spin-orbit potential :math:`\hat{V}_{\mathrm{SO}}` is expressed as a sum of separable semi-local projectors centered on each atom :math:`I`:

.. math::

   \hat{V}_{\mathrm{SO}} = \sum_I \sum_{l=1}^{l_{\max}} \sum_{m, m'=-l}^{l} \sum_{i=1}^{N_{\mathrm{proj}}} \sum_{j=1}^{N_{\mathrm{proj}}} |p_i^{Ilm}\rangle \, k_{ij}^{Il} \, \left( \mathbf{L} \cdot \mathbf{S} \right)_{mm'} \, \langle p_j^{Ilm'}|

where:
* :math:`|p_i^{Ilm}\rangle` are atom-centered Gaussian-type projector functions with angular momentum :math:`(l, m)` and radial index :math:`i`.
* :math:`k_{ij}^{Il}` are material-specific relativistic coupling coefficients tabulated in ``GTH_SOC_POTENTIALS.txt``.
* :math:`\mathbf{L} = (L_x, L_y, L_z)` are orbital angular momentum operators in the complex spherical harmonic basis.
* :math:`\mathbf{S} = \frac{1}{2} (\sigma_x, \sigma_y, \sigma_z)` are the Pauli spin matrices:

.. math::

   \sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
   \sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad
   \sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}

Angular Momentum Matrices in Spherical Harmonics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The matrix elements of :math:`L_z`, :math:`L_+`, and :math:`L_-` in the standard complex spherical harmonic basis :math:`|l, m\rangle` are:

.. math::

   \langle l, m | L_z | l, m' \rangle = m \, \delta_{m, m'}

.. math::

   \langle l, m \pm 1 | L_\pm | l, m \rangle = \sqrt{l(l+1) - m(m \pm 1)}

The Cartesian components are obtained via:

.. math::

   L_x = \frac{1}{2} (L_+ + L_-), \quad L_y = \frac{1}{2i} (L_+ - L_-)

Two-Component Spinor Hamiltonian Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Expanding the total Hamiltonian in the two-component spinor basis :math:`(\alpha, \beta)`:

.. math::

   |\psi_k^{\mathrm{spinor}}\rangle = \sum_{m=1}^{N_{\mathrm{mo}}} \left[ U_{mk}^\alpha |\phi_m\rangle \otimes |\alpha\rangle + U_{mk}^\beta |\phi_m\rangle \otimes |\beta\rangle \right]

The complete Hamiltonian takes the :math:`2 N_{\mathrm{mo}} \times 2 N_{\mathrm{mo}}` block form:

.. math::

   \mathbf{H}_{\mathrm{total}} = \begin{pmatrix}
     \mathbf{H}_0 - \frac{1}{2} \mathbf{H}_z & -\frac{1}{2} (\mathbf{H}_x - i \mathbf{H}_y) \\
     -\frac{1}{2} (\mathbf{H}_x + i \mathbf{H}_y) & \mathbf{H}_0 + \frac{1}{2} \mathbf{H}_z
   \end{pmatrix}

where :math:`\mathbf{H}_0 = \operatorname{diag}(\varepsilon_1^{\mathrm{DFT}}, \dots, \varepsilon_{N_{\mathrm{mo}}}^{\mathrm{DFT}})` is the diagonal matrix of spin-free Kohn-Sham orbital energies, and the Cartesian spin-orbit blocks in the molecular orbital active space are:

.. math::

   \mathbf{H}_x = \frac{1}{2} \mathbf{B}_{\mathrm{mo}} \, \mathbf{K}_x \, \mathbf{B}_{\mathrm{mo}}^T

.. math::

   \mathbf{H}_y = \frac{i}{2} \mathbf{B}_{\mathrm{mo}} \, \tilde{\mathbf{K}}_y \, \mathbf{B}_{\mathrm{mo}}^T

.. math::

   \mathbf{H}_z = \frac{1}{2} \mathbf{B}_{\mathrm{mo}} \, \mathbf{K}_z \, \mathbf{B}_{\mathrm{mo}}^T

Here, :math:`\mathbf{B}_{\mathrm{mo}} = \mathbf{C}_{\mathrm{act}}^T \mathbf{B}_{\mathrm{raw}}` represents the projection of the atomic orbital-to-projector overlap matrix into the active space.

---

5. High-Performance Sparse Assembly & DGEMM Optimization
--------------------------------------------------------

For large nanocrystals containing :math:`> 10,000` AOs and thousands of molecular orbitals, naively constructing the Cartesian SOC blocks via triple nested loops over all 1,750 atomic projector blocks requires minutes per geometry frame.

``QDEX`` achieves **sub-second spinor diagonalization** through three architectural optimizations:

1. **Global Sparse CSR Projector Representation**:
   All atom-centered angular projector matrices :math:`k_{ij}^{Il} (L_\kappa)_{mm'}` are pre-assembled once into global sparse Compressed Sparse Row (CSR) matrices :math:`\mathbf{K}_x, \tilde{\mathbf{K}}_y, \mathbf{K}_z`. These sparse structures remain invariant across MD steps and are cached in memory.

2. **Purely Real BLAS DGEMM**:
   Although the spin-orbit Hamiltonian is complex Hermitian, the Cartesian building blocks :math:`\mathbf{B}_{\mathrm{mo}}`, :math:`\mathbf{K}_x`, :math:`\tilde{\mathbf{K}}_y`, and :math:`\mathbf{K}_z` are **strictly real-valued** (`float64`). ``QDEX`` executes all intermediate tensor contractions using highly optimized real BLAS Level-3 DGEMM routines, completely bypassing expensive complex matrix multiplications.

3. **In-Place Contiguous Memory Allocation**:
   The four blocks of :math:`\mathbf{H}_{\mathrm{total}}` are populated directly into a single contiguous :math:`(2N_{\mathrm{act}} \times 2N_{\mathrm{act}})` complex array, eliminating auxiliary memory copies before LAPACK ``zheevd`` eigensolving.

---

6. Spinor Representation & Unrestricted Kohn-Sham (UKS)
-------------------------------------------------------

Diagonalization of :math:`\mathbf{H}_{\mathrm{total}}` yields the relativistic spinor eigenvalues :math:`\varepsilon_k^{\mathrm{spinor}}` and the unitary expansion matrix :math:`\mathbf{U}`:

.. math::

   \mathbf{H}_{\mathrm{total}} \mathbf{U} = \mathbf{U} \operatorname{diag}\left( \varepsilon_1^{\mathrm{spinor}}, \dots, \varepsilon_{2N_{\mathrm{act}}}^{\mathrm{spinor}} \right)

The :math:`k`-th spinor wavefunction is partitioned into its :math:`\alpha` and :math:`\beta` spin components:

.. math::

   \psi_k^\alpha(\mathbf{r}) = \sum_m U_{mk}^\alpha \phi_m(\mathbf{r}), \quad
   \psi_k^\beta(\mathbf{r}) = \sum_m U_{mk}^\beta \phi_m(\mathbf{r})

The local spinor probability density is given by:

.. math::

   \rho_k^{\mathrm{spinor}}(\mathbf{r}) = |\psi_k^\alpha(\mathbf{r})|^2 + |\psi_k^\beta(\mathbf{r})|^2

UKS Spin-Preserving Framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When starting from an Unrestricted Kohn-Sham (UKS) calculation with different molecular orbitals for alpha and beta spins (:math:`\mathbf{C}_\alpha \neq \mathbf{C}_\beta`), ``QDEX`` maps the active spaces independently:

.. math::

   \mathbf{H}_{\mathrm{total}}^{\mathrm{UKS}} = \begin{pmatrix}
     \mathbf{H}_{0, \alpha} - \frac{1}{2} \mathbf{H}_{z, \alpha \alpha} & -\frac{1}{2} (\mathbf{H}_{x, \alpha \beta} - i \mathbf{H}_{y, \alpha \beta}) \\
     -\frac{1}{2} (\mathbf{H}_{x, \beta \alpha} + i \mathbf{H}_{y, \beta \alpha}) & \mathbf{H}_{0, \beta} + \frac{1}{2} \mathbf{H}_{z, \beta \beta}
   \end{pmatrix}

This allows studying doped quantum dots, open-shell radicals, or spin-polarized nanocrystals without loss of relativistic accuracy.

---

7. CLI Flags & YAML Configuration Reference
-------------------------------------------

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--soc_flag``
     - ``False``
     - Enable relativistic Spin-Orbit Coupling across the calculation.
   * - ``--gth_file <path>``
     - Built-in table
     - Path to custom GTH SOC pseudopotential parameter file (overriding internal database).
   * - ``--soc_window <float>``
     - ``10.0``
     - Energy window (in eV) centered at the Fermi level for selecting MOs in the SOC active space.
   * - ``--device <choice>``
     - ``auto``
     - Compute device: ``cpu``, ``cuda``, ``mps``, or ``numpy``.
   * - ``--namd-soc``
     - ``False``
     - Enable SOC specifically for Non-Adiabatic Molecular Dynamics trajectory precomputation.

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   physics:
     soc: true
     soc_window_ev: 8.0
     gth_potentials: "GTH_SOC_POTENTIALS.txt"

   system:
     mo_file: "CsPbI3_MOs.mbse"
     xyz: "CsPbI3_QD.xyz"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-SR-GTH"
