Part 4: Optical Excitations & Four Excited-State Frameworks
============================================================

The description of neutral optical excitations in semiconductor nanostructures requires treating the two-particle correlated motion of an electron promoted to the conduction band and the hole left behind in the valence band.

``QDEX`` provides four distinct theoretical frameworks for computing excited states, ranging from non-interacting single-particle transitions to the fully coupled **Bethe-Salpeter Equation (BSE)** under the Tamm-Dancoff Approximation (TDA), coupled with a clear separation between **Two-Electron Integral Representations** and **Dielectric Screening Kernels**.

---

1. The Two-Particle Excitation Problem
--------------------------------------

In a neutral optical excitation, an electron is removed from an occupied valence orbital :math:`i` and placed into an unoccupied conduction orbital :math:`a`, creating an electron-hole pair configuration :math:`|ia\rangle = a_a^\dagger a_i |0\rangle`.

The exact correlated excited-state wavefunction :math:`|\Psi_S\rangle` is expressed as a linear superposition of electron-hole configurations:

.. math::

   |\Psi_S\rangle = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, |ia\rangle

where :math:`X_{ia}^S` are the configuration interaction amplitudes and :math:`\Omega_S` is the corresponding optical excitation energy.

Tamm-Dancoff Approximation (TDA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In linear-response many-body perturbation theory, the full Bethe-Salpeter Equation contains both resonant excitations (:math:`\mathbf{A}`) and anti-resonant de-excitations (:math:`\mathbf{B}`):

.. math::

   \begin{pmatrix}
     \mathbf{A} & \mathbf{B} \\
     -\mathbf{B}^* & -\mathbf{A}^*
   \end{pmatrix}
   \begin{pmatrix}
     \mathbf{X} \\
     \mathbf{Y}
   \end{pmatrix}
   = \Omega
   \begin{pmatrix}
     \mathbf{X} \\
     \mathbf{Y}
   \end{pmatrix}.

Under the **Tamm-Dancoff Approximation (TDA)**, coupling to ground-state de-excitations is neglected (:math:`\mathbf{B} \approx \mathbf{0}`). This reduces the problem to a standard Hermitian eigenvalue problem:

.. math::

   \mathbf{A} \mathbf{X}_S = \Omega_S \mathbf{X}_S.

The TDA is exceptionally robust for semiconductor nanostructures: it eliminates triplet instabilities, guarantees purely real excitation energies, and reduces computational complexity by a factor of 4 with negligible loss of accuracy for optical transitions well below the plasma frequency.

Singlet vs. Triplet Matrix Elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The elements of the resonant matrix :math:`A_{ia, jb}` depend on the spin multiplicity:

.. math::

   A_{ia, jb}^{\mathrm{singlet}} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} + 2 K_{ia, jb}^x - K_{ia, jb}^d

.. math::

   A_{ia, jb}^{\mathrm{triplet}} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} - K_{ia, jb}^d.

The bare exchange term :math:`K_{ia, jb}^x` is strictly absent in triplet states because electrons with parallel spins experience identical spatial exchange cancellation. The factor of :math:`2 K_{ia, jb}^x` in singlets is responsible for the singlet-triplet exchange splitting.

Relativistic 2-Component Spinor BSE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When spin-orbit coupling is enabled (``soc: true``), the single-particle spatial orbitals are expanded into 2-component Kramers spinors :math:`\psi_p(\mathbf{r}) = \begin{pmatrix} \phi_{p\alpha}(\mathbf{r}) \\ \phi_{p\beta}(\mathbf{r}) \end{pmatrix}`. The BSE Hamiltonian is transformed into the spinor electron-hole basis :math:`|IA\rangle = a_A^\dagger a_I |0\rangle`:

.. math::

   A_{IA, JB}^{\mathrm{spinor}} = \left( \varepsilon_A^{\mathrm{QP}} - \varepsilon_I^{\mathrm{QP}} \right) \delta_{IJ} \delta_{AB} + K_{IA, JB}^x - K_{IA, JB}^d

capturing fine-structure splittings, bright-dark exciton order inversion, and Rashba effects without empirical parameters.

---

2. Key Distinction: ``2e-integrals`` vs. ``kernel``
---------------------------------------------------

A critical conceptual distinction in ``QDEX`` is the clean separation between the **Two-Electron Integral Representation** and the **Dielectric Screening Kernel**.

As introduced in :doc:`../part3_gw_scissor/index`, both the single-particle self-energy :math:`\Sigma = iGW` and the two-particle Bethe-Salpeter kernel :math:`\Xi = \delta \Sigma / \delta G \approx v - W` share the **exact same two-body interaction engine**:

.. math::

   \Xi(\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3, \mathbf{r}_4) \approx \underbrace{v(\mathbf{r}_1, \mathbf{r}_3) \, \delta(\mathbf{r}_1 - \mathbf{r}_2) \, \delta(\mathbf{r}_3 - \mathbf{r}_4)}_{\text{Bare Exchange } K^x} - \underbrace{W(\mathbf{r}_1, \mathbf{r}_2) \, \delta(\mathbf{r}_1 - \mathbf{r}_3) \, \delta(\mathbf{r}_2 - \mathbf{r}_4)}_{\text{Direct Screened Attraction } K^d}.

In ``QDEX``, these two operators are configured independently:

.. list-table::
   :widths: 22 28 50
   :header-rows: 1

   * - Keyword
     - Options
     - Physical Role in BSE
   * - ``2e-integrals``
     - ``mnok``, ``xs`` (or ``xs-qdex``)
     - Defines the **bare Coulomb representation** :math:`v(\mathbf{r}_1, \mathbf{r}_2)` used to evaluate the electron-hole exchange matrix :math:`K^x` (governing singlet-triplet splitting).
   * - ``kernel``
     - ``resta``, ``dim``, ``rpa``, ``sbse``, ``bse``
     - Defines the **physical dielectric screening profile** :math:`W(\mathbf{r}, \mathbf{r}')` used to evaluate the direct electron-hole attraction matrix :math:`K^d` (binding the exciton).

.. code-block:: text

   +-----------------------------------------------------------------------------------------+
   |                          Unified Two-Body Interaction Engine                            |
   |                                                                                         |
   |   Bare Coulomb Operator v                        Screened Interaction W = ε⁻¹ v         |
   |   (keyword: '2e-integrals')                      (keyword: 'kernel')                    |
   |   • mnok (Ohno point charges)                    • resta (Thomas-Fermi + Penn scaling)  |
   |   • xs   (Libint2 4-center Gaussians)            • dim   (Atomistic polarizable dipoles)|
   |                                                  • rpa   (Microscopic ZDO polarizability|
   |                                                  • sbse  (Simplified BSE)               |
   |                                                  • bse   (Uniform bulk ε_∞ scaling)     |
   +----------------------------+-----------------------------------+------------------------+
                                |                                   |
                                | (Exchange Splitting)              | (Direct Attraction)
                                v                                   v
                +-------------------------------+   +-------------------------------+
                |    Bare Exchange Matrix K_x   |   |   Screened Direct Matrix K_d  |
                +---------------+---------------+   +---------------+---------------+
                                |                                   |
                                +-----------------+-----------------+
                                                  |
                                                  v
                                +-----------------------------------+
                                |     Bethe-Salpeter Matrix A       |
                                |                                   |
                                | Singlet: (ε_a - ε_i) + 2 K_x - K_d|
                                | Triplet: (ε_a - ε_i) - K_d        |
                                | Spinor:  (ε_A - ε_I) + K_x - K_d  |
                                +-----------------------------------+

.. note::
   The keyword ``2e-integrals`` is also accepted via YAML as ``two_electron_integrals`` or ``2e_integrals``, and on the CLI as ``--2e-integrals`` or ``--two-electron-integrals``. The legacy flag ``kernel_type`` is retained for full backward compatibility.

---

3. Two-Electron Integral Representations (``2e-integrals``)
-----------------------------------------------------------

Semi-Empirical Atom-Centered Representation (``2e-integrals: mnok``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To achieve high efficiency for nanocrystals containing up to 10,000 atoms, ``QDEX`` contracts transition densities into atom-centered point charges using Mulliken or Lowdin population analysis:

.. math::

   q_A^{ia} = \sum_{\mu \in A} \sum_{\nu=1}^{N_{\mathrm{ao}}} C_{\mu i} S_{\mu \nu} C_{\nu a}.

For diagonal electron and hole charge densities:

* Hole density on atom :math:`A`: :math:`q_A^{ii} = \sum_{\mu \in A} \sum_\nu C_{\mu i} S_{\mu \nu} C_{\nu i}`
* Electron density on atom :math:`B`: :math:`q_B^{aa} = \sum_{\mu \in B} \sum_\nu C_{\mu a} S_{\mu \nu} C_{\nu a}`

The four-center integrals are replaced by pairwise contractions over atomic sites:

.. math::

   K_{ia, jb}^x = \sum_{A=1}^{N_{\mathrm{atoms}}} \sum_{B=1}^{N_{\mathrm{atoms}}} q_A^{ia} \, \gamma_{AB}^{\mathrm{Ohno}} \, q_B^{jb}

.. math::

   K_{ia, jb}^d = \sum_{A=1}^{N_{\mathrm{atoms}}} \sum_{B=1}^{N_{\mathrm{atoms}}} q_A^{ij} \, W_{AB} \, q_B^{ab}

where :math:`\gamma_{AB}^{\mathrm{Ohno}}` is the Mataga-Nishimoto-Ohno-Klopman (MNOK) damped Coulomb potential:

.. math::

   \gamma_{AB}^{\mathrm{Ohno}} = \frac{1}{\sqrt{R_{AB}^2 + a_{AB}^2}}

and :math:`a_{AB} = 2 / (\eta_A + \eta_B)` is the Ohno-Klopman damping parameter derived from atomic chemical hardnesses :math:`\eta_A` and :math:`\eta_B`.

**Advantages**: Scales strictly as :math:`O(N_{\mathrm{atoms}}^2)`. Memory requirements are minimal (:math:`< 10\text{ MB}`).

Exact Analytical Gaussian Representation (``2e-integrals: xs``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For smaller clusters (:math:`\le 500` atoms) or benchmark comparisons, ``QDEX`` evaluates the **exact four-center two-electron Gaussian repulsion integrals** analytically via the C++ Libint2 library:

.. math::

   (\mu \nu | \lambda \sigma) = \iint \chi_\mu(\mathbf{r}_1) \chi_\nu(\mathbf{r}_1) \, \frac{1}{|\mathbf{r}_1 - \mathbf{r}_2|} \, \chi_\lambda(\mathbf{r}_2) \chi_\sigma(\mathbf{r}_2) \, d\mathbf{r}_1 \, d\mathbf{r}_2.

The two-electron atomic repulsion matrix is computed over AO density pairs:

.. math::

   \Gamma_{\mu \nu} = (\mu \mu | \nu \nu).

The molecular orbital matrix elements are constructed by exact AO-to-MO contraction:

.. math::

   K_{ia, jb}^x = \sum_{\mu, \nu=1}^{N_{\mathrm{ao}}} \left( C_{\mu i} C_{\nu a} \right) \Gamma_{\mu \nu} \left( C_{\mu j} C_{\nu b} \right)

.. math::

   K_{ia, jb}^d = \sum_{\mu, \nu=1}^{N_{\mathrm{ao}}} \left( C_{\mu i} C_{\nu j} \right) W_{\mu \nu} \left( C_{\mu a} C_{\nu b} \right).

**Advantages**: Preserves non-spherical orbital angular momentum components (e.g. anisotropic :math:`p` and :math:`d` orbital bonding interactions), eliminating any reliance on spherical atomic charge partitioning.

Implementation in QDEX (``2e-integrals``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-electron integral representations are managed across :mod:`qdex.hardness`, :mod:`qdex.integrals`, and :mod:`qdex.exciton_hamiltonian`:

1. **Semi-Empirical Atom-Centered Points (``2e-integrals: mnok``)**:
   - Implemented via :func:`qdex.hardness.build_damped_mnok_matrix` and :class:`qdex.exciton_hamiltonian.ExcitonHamiltonian`.
   - Projects molecular orbitals onto atomic centers using Löwdin or Mulliken populations: :math:`\mathbf{q}_A^{ia} = \sum_{\mu \in A, \nu} C_{\mu i} S_{\mu \nu} C_{\nu a}` via optimized BLAS ``DGEMM`` routines.
   - Pairs transition charges with the damped Ohno matrix :math:`\mathbf{\Gamma}_{AB}^{\mathrm{Ohno}}` for exchange and :math:`\mathbf{W}_{AB}` for direct screening.
   - Fast, memory-lean (< 10 MB RAM), and scales easily to nanocrystals containing 1,000 to 10,000 atoms.

2. **Analytical 4-Center GTO Integrals (``2e-integrals: xs``)**:
   - Evaluated via the C++ extension :mod:`libint_cpp` and wrapped in :func:`qdex.integrals.compute_two_electron_ao`.
   - Calculates exact two-electron Gaussian repulsion integrals :math:`(\mu \mu | \nu \nu)` over contracted GTO basis shells.
   - Preserves complete angular orbital anisotropy without spherical approximations. Recommended for molecular benchmarks and small nanoclusters (:math:`\le 500` atoms).

---

4. Dielectric Screening Kernels (``kernel``)
--------------------------------------------

The direct electron-hole attraction :math:`K^d` is mediated by the screened interaction :math:`W`. ``QDEX`` provides five distinct dielectric screening kernels:

1. Resta Screened Dielectric Kernel (``kernel: resta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In semiconductor nanoclusters, dielectric screening varies continuously from the bulk optical dielectric constant :math:`\epsilon_\infty` at large distances down to unscreened vacuum interaction (:math:`\epsilon = 1`) at short, intra-atomic distances.

``QDEX`` implements the **Resta model of electronic screening**:

.. math::

   W(r) = \frac{1}{\epsilon_\infty r} + \frac{1 - \epsilon_\infty^{-1}}{r} \exp\left( -\frac{r}{\lambda_s} \right)

where :math:`\lambda_s` is the **Thomas-Fermi screening length** of the valence electron gas:

.. math::

   \lambda_s = \sqrt{ \frac{\pi}{4 k_F} } = \left( \frac{\pi}{4 (3\pi^2 n_v)^{1/3}} \right)^{1/2}.

Damped over atomic centers with Ohno-Klopman hardness parameters, the discrete kernel is:

.. math::

   W_{AB}^{\mathrm{Resta}} = \frac{1}{\epsilon_\infty \sqrt{R_{AB}^2 + a_{AB}^2}} + \frac{1 - \epsilon_\infty^{-1}}{\sqrt{R_{AB}^2 + a_{AB}^2}} \exp\left( -\frac{\sqrt{R_{AB}^2 + a_{AB}^2}}{\lambda_s} \right).

2. Atomistic Discrete Dipole Interaction Kernel (``kernel: dim``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluates screening via the Discrete Dipole Interaction Model (DIM / Thole model). Each atom :math:`A` responds with an induced dipole :math:`\boldsymbol{\mu}_A = \alpha_A \mathbf{E}^{\mathrm{tot}}(\mathbf{R}_A)` to the electron-hole charge distribution:

.. math::

   (\mathbf{I} + \boldsymbol{\alpha}\mathbf{T}) \boldsymbol{\mu} = \boldsymbol{\alpha}\mathbf{E}_0

yielding an atom-specific screened potential :math:`W_{AB}^{\mathrm{DIM}} = S_{AB} \, \gamma_{AB}`.

3. Parameter-Free Microscopic ZDO-RPA Kernel (``kernel: rpa`` / ``xs-rpa``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Computes microscopic dielectric screening directly from the non-interacting transition polarizability of the quantum dot without any empirical parameters:

.. math::

   \Pi_{\mu \nu}^0 = 4 \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} \frac{C_{\mu i} C_{\mu a} C_{\nu i} C_{\nu a}}{\varepsilon_a - \varepsilon_i}

.. math::

   \boldsymbol{\epsilon} = \mathbf{I} + \boldsymbol{\Gamma} \boldsymbol{\Pi}^0, \quad \mathbf{W} = \boldsymbol{\epsilon}^{-1} \boldsymbol{\Gamma}.

4. Simplified BSE Kernel (``kernel: sbse``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Constructs the screened interaction according to Cho, Bintrim, and Berkelbach:

.. math::

   \mathbf{W} = (\mathbf{I} + \mathbf{J}_{\mathrm{solv}} \boldsymbol{\Pi}^0)^{-1} \mathbf{J}_{\mathrm{solv}}.

5. Uniform Dielectric Kernel (``kernel: bse``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scales the interaction uniformly by :math:`1/\epsilon_\infty` with an empirical scaling factor :math:`\alpha` (default :math:`\alpha = 1.0`).

Spatial Asymptotics & Wannier-Mott Bulk Limit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As established in the quasiparticle theory (:doc:`../part3_gw_scissor/index`), the screened interaction :math:`W(r)` connects short-range atomic scales to the macroscopic crystal:

1. **Short-Range Limit (:math:`r \to 0`)**: :math:`W(r) \to v(r) = 1/r` (:math:`\epsilon \to 1`). Electronic screening ceases at sub-nanometer distances because the valence electrons cannot instantaneously compress inside an atomic core. This prevents the catastrophic underestimation of singlet-triplet exchange splitting and on-site Coulomb repulsion.
2. **Nanocrystal Boundary (:math:`r \sim R_{\mathrm{QD}}`)**: Dielectric mismatch between the dot (:math:`\epsilon_\infty`) and the solvent (:math:`\epsilon_{\mathrm{out}}`) generates an image-charge reaction field :math:`W^{\mathrm{solv}}`. In low-permittivity solvents, dielectric confinement strongly enhances the direct electron-hole attraction :math:`K^d`.
3. **Asymptotic Bulk Limit (:math:`r \to \infty`, :math:`R_{\mathrm{QD}} \to \infty`)**: :math:`W(r) \to \frac{1}{\epsilon_\infty r}`. In this limit, the Bethe-Salpeter equation continuously reduces to the hydrogenic Wannier-Mott exciton equation:

   .. math::

      \left( -\frac{\hbar^2 \nabla_{\mathbf{r}}^2}{2\mu} - \frac{e^2}{\epsilon_\infty r} \right) \phi_{\mathrm{exc}}(\mathbf{r}) = -E_b^{\mathrm{bulk}} \phi_{\mathrm{exc}}(\mathbf{r})

   recovering the bulk Rydberg binding energy:

   .. math::

      E_b^{\mathrm{bulk}} = \frac{\mu e^4}{2 \hbar^2 \epsilon_\infty^2} = R_y^*

   with bulk exciton Bohr radius :math:`a_{\mathrm{exc}} = a_0 \epsilon_\infty (m_0 / \mu)`. Consequently, as the nanocrystal diameter surpasses the Bohr radius (:math:`R_{\mathrm{QD}} \gg a_{\mathrm{exc}}`), the BSE optical transition energy :math:`\Omega_1` smoothly converges to the bulk band edge minus the Wannier-Mott binding energy: :math:`\Omega_1 \to E_g^{\mathrm{bulk}} - E_b^{\mathrm{bulk}}`.

---

5. The Four Excitation Frameworks (``excitation_mode``)
-------------------------------------------------------

``QDEX`` provides four progressive levels of physical theory via ``--excitation-mode``:

.. list-table::
   :widths: 22 28 50
   :header-rows: 1

   * - Framework Mode
     - Energy Expression
     - Physical Characteristics
   * - **(A) `independent_dft`**
     - :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{DFT}} - \varepsilon_i^{\mathrm{DFT}}`
     - Bare Kohn-Sham transitions. Completely ignores quasiparticle self-energy corrections and electron-hole Coulomb interactions.
   * - **(B) `independent_qp`**
     - :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}`
     - Non-interacting quasiparticles. Applies the scaled GW scissor shift :math:`\Delta_{\mathrm{GW}}`, opening the gap to experimental values, but neglects electron-hole binding (:math:`E_b = 0`).
   * - **(C) `diagonal_bse`**
     - :math:`\Omega_{ia} = (\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}) + 2 K_{ia,ia}^x - K_{ia,ia}^d`
     - Diagonal BSE. Accounts for both quasiparticle self-energy and diagonal electron-hole Coulomb binding, but omits off-diagonal configuration mixing. Extremely fast.
   * - **(D) `bse` (sTDA)**
     - Full diagonalization of :math:`A_{ia, jb}`
     - Fully coupled configuration interaction. Solves the complete resonant matrix, capturing spatial exciton delocalization, state mixing, and oscillator strength redistribution.

Why Diagonal BSE Works in Nanocrystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In standard BSE, diagonalizing the full :math:`A_{ia, jb}` matrix requires :math:`O(N_{\mathrm{pairs}}^3)` operations. In a quantum dot with 100 occupied and 100 virtual orbitals, :math:`N_{\mathrm{pairs}} = 10,000`, requiring gigabytes of memory and long compute times.

The **Diagonal BSE** framework omits off-diagonal configuration interaction (:math:`ia \neq jb`), evaluating:

.. math::

   \Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} + 2 K_{ia, ia}^x - K_{ia, ia}^d.

* **Dominance of Diagonal Coulomb Attraction**: In quantum dots, the diagonal term :math:`K_{ia, ia}^d` represents the direct electrostatic attraction between the electron distribution :math:`|\phi_a|^2` and hole distribution :math:`|\phi_i|^2`, accounting for **over 90% of the total exciton binding energy** :math:`E_b`.
* **Essential for NAMD**: In non-adiabatic molecular dynamics simulations where excited states must be evaluated at every time step (e.g. 5,000 steps), full BSE diagonalization is computationally prohibitive. Diagonal BSE provides an accurate, energy-conserving potential energy surface at a fraction of the cost.

Implementation in QDEX (``--excitation-mode``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The excitation frameworks are implemented across :mod:`qdex.cli`, :mod:`qdex.exciton_hamiltonian`, and :mod:`qdex.solver`:

1. **Bare DFT Transitions (``independent_dft``)**:
   Returns raw differences :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{DFT}} - \varepsilon_i^{\mathrm{DFT}}` directly from Kohn-Sham eigenvalues without building two-body Coulomb matrices.
2. **Quasiparticle Transitions (``independent_qp``)**:
   Applies the QP scissor shifts :math:`\varepsilon^{\mathrm{QP}}` from Part 3: :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}`. Evaluates single-particle transition dipoles :math:`\boldsymbol{\mu}_{ia}` and oscillator strengths without electron-hole binding.
3. **Diagonal BSE (``diagonal_bse``)**:
   Implemented via :func:`qdex.solver.solve_diagonal_bse`. Builds only the diagonal elements :math:`A_{ia, ia} = (\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}) + 2 K_{ia,ia}^x - K_{ia,ia}^d` in :math:`O(N_{\mathrm{pairs}})` time. Completely bypasses matrix diagonalization, making it the ideal engine for multi-thousand step non-adiabatic dynamics.
4. **Coupled Bethe-Salpeter Equation (``bse``)**:
   Implemented via :func:`qdex.solver.davidson`. Constructs the active space transition basis, applies energy truncation thresholds (``--ci_threshold``), and solves for the lowest :math:`N_{\mathrm{roots}}` exciton eigenvectors using a block-Davidson iterative subspace algorithm.

---

6. Transition Dipoles, Oscillator Strengths & Superradiance
-----------------------------------------------------------

The single-particle transition dipole moment between occupied orbital :math:`i` and virtual orbital :math:`a` is:

.. math::

   \boldsymbol{\mu}_{ia} = \langle \phi_i | \mathbf{r} | \phi_a \rangle = \sum_{\mu \nu} C_{\mu i} \, \mathbf{D}_{\mu \nu} \, C_{\nu a}

where :math:`\mathbf{D}_{\mu \nu} = \langle \chi_\mu | \mathbf{r} | \chi_\nu \rangle` is the AO dipole matrix evaluated analytically via Libint2.

Under full configuration interaction, the exciton transition dipole moment is a coherent linear superposition:

.. math::

   \boldsymbol{\mu}_S = \sum_{ia} X_{ia}^S \, \boldsymbol{\mu}_{ia}.

The corresponding dimensionless oscillator strength is:

.. math::

   f_S = \frac{2}{3} \, \Omega_S \, |\boldsymbol{\mu}_S|^2.

This coherent summation describes **superradiance** and intensity borrowing, where optical strength from high-energy transitions is transferred into the lowest bright exciton.

---

7. Full BSE & The Davidson Iterative Solver
-------------------------------------------

To solve for the lowest :math:`k` roots without dense :math:`O(N_{\mathrm{pairs}}^3)` matrix diagonalization, ``QDEX`` implements an optimized **Davidson iterative subspace solver**:

1. Projects :math:`\mathbf{A}` into a small trial subspace :math:`\mathbf{V} = [\mathbf{v}_1, \dots, \mathbf{v}_m]`.
2. Computes the matrix-vector product :math:`\mathbf{w}_j = \mathbf{A} \mathbf{v}_j` on-the-fly.
3. Solves the projected eigenvalue problem :math:`\mathbf{V}^\dagger \mathbf{A} \mathbf{V} \mathbf{y} = \omega \mathbf{y}`.
4. Computes the residual :math:`\mathbf{r} = \mathbf{A} \mathbf{x} - \omega \mathbf{x}` and preconditioner :math:`\boldsymbol{\delta} = (\operatorname{diag}(\mathbf{A}) - \omega)^{-1} \mathbf{r}`.
5. Expands the subspace until the norm :math:`||\mathbf{r}|| < \text{tol}` (default :math:`10^{-5}`).

This reduces memory requirements from :math:`O(N_{\mathrm{pairs}}^2)` to :math:`O(N_{\mathrm{pairs}} \times k)`.

---

8. Comprehensive Exciton Calculation Matrix
-------------------------------------------

.. list-table::
   :widths: 18 15 15 14 38
   :header-rows: 1

   * - ``excitation_mode``
     - ``2e-integrals``
     - ``kernel``
     - Complexity
     - Recommended Purpose
   * - ``bse``
     - ``mnok``
     - ``resta``
     - :math:`O(N_{\mathrm{atoms}}^2)`
     - **Production Standard**: Absorption spectra for large colloidal quantum dots (up to 10,000 atoms).
   * - ``bse``
     - ``mnok``
     - ``dim``
     - :math:`O(N_{\mathrm{atoms}}^2)`
     - Core/shell or anisotropic nanocrystals with polarizable dielectric boundaries.
   * - ``diagonal_bse``
     - ``mnok``
     - ``resta``
     - :math:`O(N_{\mathrm{pairs}})`
     - **NAMD Production**: Ultrafast excited-state dynamics and non-adiabatic trajectories.
   * - ``bse``
     - ``xs``
     - ``rpa``
     - :math:`O(N_{\mathrm{ao}}^4)`
     - **Ab Initio Benchmark**: Parameter-free exact Gaussian integrals and microscopic RPA screening.
   * - ``bse``
     - ``xs``
     - ``dim``
     - :math:`O(N_{\mathrm{ao}}^4)`
     - Exact Gaussian integrals with atomistic polarizable dipole screening.
   * - ``independent_qp``
     - None
     - None
     - :math:`O(N_{\mathrm{pairs}})`
     - Single-particle joint density of states with quasiparticle gap correction.
   * - ``independent_dft``
     - None
     - None
     - :math:`O(1)`
     - Uncorrected baseline Kohn-Sham single-particle transitions.

---

9. Recommended Workflow Presets
-------------------------------

Preset 1: Standard Colloidal QD Absorption Spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended default for production calculations on colloidal nanocrystals:

.. code-block:: yaml

   physics:
     excitation_mode: "bse"
     2e-integrals: "mnok"              # Fast semi-empirical atom-centered representation
     kernel: "resta"                   # Thomas-Fermi electronic screening
     charge_type: "mulliken"           # "mulliken" or "lowdin"
     qp_gap: "sgw-dim"                 # Microscopic polarizable dipole QP gap
     dynamic_z: true                   # State-dependent Z_p via PPM f-sum rule
     nhomos: 50
     nlumos: 50

   bse:
     nroots: 15
     full_diag: false
     tol: 1.0e-5

Preset 2: Ultrafast Non-Adiabatic Molecular Dynamics (NAMD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optimized for high-speed trajectory propagation across thousands of frames:

.. code-block:: yaml

   physics:
     excitation_mode: "diagonal_bse"   # Diagonal e-h attraction without CI mixing
     2e-integrals: "mnok"
     kernel: "resta"
     qp_gap: "sgw-anchor"              # Sub-millisecond two-anchor scaled GW
     nhomos: 30
     nlumos: 30

Preset 3: Benchmark First-Principles Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exact analytical Gaussian two-electron integrals with microscopic RPA screening:

.. code-block:: yaml

   physics:
     excitation_mode: "bse"
     2e-integrals: "xs"                # Exact 4-center Gaussian integrals via Libint2
     kernel: "xs-rpa"                  # Parameter-free microscopic ZDO-RPA screening
     qp_gap: "qsgw-dim"                # Full AO orbital relaxation (qsGW)
     update_orbitals: true
     nhomos: 25
     nlumos: 25

   bse:
     nroots: 10
     full_diag: false

Preset 4: Spin-Orbit Coupling & Dark Excitons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Relativistic 2-component spinor Bethe-Salpeter calculation:

.. code-block:: yaml

   physics:
     excitation_mode: "bse"
     2e-integrals: "mnok"
     kernel: "resta"
     soc: true                         # Enable 2-component spinor Hamiltonian
     soc_window_ev: 8.0                # Active window around Fermi level in eV
     qp_gap: "sgw-anchor"
     nhomos: 40
     nlumos: 40

   bse:
     nroots: 20
     full_diag: false

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
   * - ``--excitation-mode <choice>``
     - ``bse``
     - Excitation framework: ``bse`` (full sTDA), ``diagonal_bse``, ``independent_qp``, or ``independent_dft``.
   * - ``--2e-integrals <choice>``
     - ``mnok``
     - Two-electron integral representation: ``mnok`` (semi-empirical atom-centered) or ``xs`` (exact analytical Gaussian integrals).
   * - ``--kernel <choice>``
     - ``bse``
     - Dielectric screening kernel: ``resta``, ``dim``, ``rpa``, ``sbse``, ``xs-resta``, ``xs-dim``, ``xs-rpa``, or ``bse``.
   * - ``--nhomos <int>``
     - ``25``
     - Number of occupied frontier orbitals to include in the BSE active space.
   * - ``--nlumos <int>``
     - ``25``
     - Number of virtual frontier orbitals to include in the BSE active space.
   * - ``--e_thresh <float>``
     - ``None``
     - Energy threshold (in eV) to automatically select active pairs with :math:`\varepsilon_a - \varepsilon_i \le E_{\mathrm{thresh}}`.
   * - ``--f_thresh <float>``
     - ``0.0``
     - Minimum oscillator strength threshold to print and log excited states.
   * - ``--nroots <int>``
     - ``10``
     - Number of lowest exciton roots to compute via Davidson diagonalization.
   * - ``--full-diag``
     - ``False``
     - Force full dense LAPACK diagonalization instead of the iterative Davidson solver.
   * - ``--tol <float>``
     - ``1e-5``
     - Convergence tolerance for the Davidson solver residual norm.
   * - ``--triplet``
     - ``False``
     - Perform triplet excited-state BSE calculation (omitting :math:`2K^x`).
   * - ``--charge_type <choice>``
     - ``mulliken``
     - Transition charge partitioning: ``mulliken`` or ``lowdin``.
