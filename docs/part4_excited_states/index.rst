Part 4: Optical Excitations & Four Excited-State Frameworks
============================================================

The description of neutral optical excitations in semiconductor nanostructures requires treating the two-particle correlated motion of an electron promoted to the conduction band and the hole left behind in the valence band.

``QDEX`` provides four distinct theoretical frameworks for computing excited states, ranging from non-interacting single-particle transitions to the fully coupled **Bethe-Salpeter Equation (BSE)** under the Tamm-Dancoff Approximation (TDA).

---

1. The Two-Particle Excitation Problem
--------------------------------------

In a neutral excitation, an electron is removed from an occupied valence orbital :math:`i` and placed into an unoccupied conduction orbital :math:`a`, creating an electron-hole configuration :math:`|ia\rangle = a_a^\dagger a_i |0\rangle`.

The exact correlated excited-state wavefunction :math:`|\Psi_S\rangle` is expressed as a linear combination of electron-hole configurations:

.. math::

   |\Psi_S\rangle = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, |ia\rangle

where :math:`X_{ia}^S` are the configuration interaction amplitudes and :math:`\Omega_S` is the corresponding optical excitation energy.

Tamm-Dancoff Approximation (TDA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In linear-response many-body theory, the full Bethe-Salpeter Equation contains both resonant excitations (:math:`\mathbf{A}`) and anti-resonant de-excitations (:math:`\mathbf{B}`):

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
   \end{pmatrix}

Under the **Tamm-Dancoff Approximation (TDA)**, coupling to ground-state de-excitations is neglected (:math:`\mathbf{B} \approx \mathbf{0}`). This reduces the problem to a standard Hermitian eigenvalue problem:

.. math::

   \mathbf{A} \mathbf{X}_S = \Omega_S \mathbf{X}_S

The TDA is exceptionally robust for semiconductor nanostructures: it eliminates triplet instabilities, guarantees purely real excitation energies, and reduces the computational complexity by a factor of 4 with negligible loss of accuracy for optical transitions well below the plasma frequency.

Singlet vs. Triplet Matrix Elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The elements of the resonant matrix :math:`A_{ia, jb}` depend on the spin state:

.. math::

   A_{ia, jb}^{\mathrm{singlet}} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} + 2 K_{ia, jb}^x - K_{ia, jb}^d

.. math::

   A_{ia, jb}^{\mathrm{triplet}} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} - K_{ia, jb}^d

The bare exchange term :math:`K_{ia, jb}^x` is strictly absent in triplet states because electrons with parallel spins experience identical spatial exchange cancellation. The factor of :math:`2 K_{ia, jb}^x` in singlets is responsible for the singlet-triplet exchange splitting.

---

2. In-Depth Examination of the Four Frameworks
----------------------------------------------

``QDEX`` enables selecting among four progressive levels of physical theory via ``--excitation-mode``:

.. list-table::
   :widths: 20 25 55
   :header-rows: 1

   * - Framework Mode
     - Energy Formula
     - Physical Characteristics
   * - **(A) `independent_dft`**
     - :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{DFT}} - \varepsilon_i^{\mathrm{DFT}}`
     - Bare Kohn-Sham transitions. Completely ignores quasiparticle self-energy corrections and electron-hole Coulomb interactions. Severely underestimates optical gaps.
   * - **(B) `independent_qp`**
     - :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}`
     - Non-interacting quasiparticles. Applies the scaled GW scissor shift :math:`\Delta_{\mathrm{GW}}`, correctly opening the band gap, but neglects electron-hole binding energy (:math:`E_b = 0`).
   * - **(C) `diagonal_bse`**
     - :math:`\Omega_{ia} = (\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}) + 2 K_{ia,ia}^x - K_{ia,ia}^d`
     - Diagonal BSE. Accounts for both quasiparticle self-energy and diagonal electron-hole Coulomb binding, but omits off-diagonal configuration mixing. Extremely fast.
   * - **(D) `bse` (sTDA)**
     - Full diagonalization of :math:`A_{ia, jb}`
     - Fully coupled configuration interaction. Solves the complete resonant matrix, capturing spatial exciton delocalization, state mixing, and oscillator strength redistribution.

---

3. Continuous 4-Center BSE vs. Atom-Centered MNOK / Resta Approximation
-----------------------------------------------------------------------

Standard Full 4-Center BSE
~~~~~~~~~~~~~~~~~~~~~~~~~~

In conventional quantum chemistry and solid-state codes, the bare exchange :math:`K^x` and screened direct attraction :math:`K^d` are defined as continuous four-center electron-hole integrals:

.. math::

   K_{ia, jb}^x = \iint \phi_i^*(\mathbf{r}) \phi_a(\mathbf{r}) \, \frac{1}{|\mathbf{r} - \mathbf{r}'|} \, \phi_j(\mathbf{r}') \phi_b^*(\mathbf{r}') \, d\mathbf{r} \, d\mathbf{r}'

.. math::

   K_{ia, jb}^d = \iint \phi_i^*(\mathbf{r}) \phi_j(\mathbf{r}) \, W(\mathbf{r}, \mathbf{r}') \, \phi_a^*(\mathbf{r}') \phi_b(\mathbf{r}') \, d\mathbf{r} \, d\mathbf{r}'

Evaluating these integrals requires transforming four-center two-electron atomic orbital integrals :math:`(\mu \nu | \lambda \sigma)` into the molecular orbital basis, an :math:`O(N_{\mathrm{ao}}^5)` operation that requires petabytes of storage for nanocrystals containing thousands of atoms.

The Atom-Centered MNOK / Resta Approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To achieve high efficiency without sacrificing dielectric physics, ``QDEX`` contracts transition densities into atom-centered point charges using Mulliken population analysis:

.. math::

   q_A^{ia} = \sum_{\mu \in A} \sum_{\nu=1}^{N_{\mathrm{ao}}} C_{\mu i} S_{\mu \nu} C_{\nu a}

For diagonal electron-hole densities:
* Hole density on atom :math:`A`: :math:`q_A^{ii} = \sum_{\mu \in A} \sum_\nu C_{\mu i} S_{\mu \nu} C_{\nu i}`
* Electron density on atom :math:`B`: :math:`q_B^{aa} = \sum_{\mu \in B} \sum_\nu C_{\mu a} S_{\mu \nu} C_{\nu a}`

The four-center integrals are replaced by pairwise contractions over atomic sites:

.. math::

   K_{ia, jb}^x = \sum_{A=1}^{N_{\mathrm{atoms}}} \sum_{B=1}^{N_{\mathrm{atoms}}} q_A^{ia} \, \gamma_{AB}^{\mathrm{Ohno}} \, q_B^{jb}

.. math::

   K_{ia, jb}^d = \sum_{A=1}^{N_{\mathrm{atoms}}} \sum_{B=1}^{N_{\mathrm{atoms}}} q_A^{ij} \, W_{AB}^{\mathrm{Resta}} \, q_B^{ab}

where :math:`\gamma_{AB}^{\mathrm{Ohno}}` is the Mataga-Nishimoto-Ohno-Klopman (MNOK) damped Coulomb potential:

.. math::

   \gamma_{AB}^{\mathrm{Ohno}} = \frac{1}{\sqrt{R_{AB}^2 + a_{AB}^2}}

and :math:`a_{AB} = 2 / (\eta_A + \eta_B)` is the Ohno-Klopman damping parameter derived from the atomic chemical hardness :math:`\eta_A` and :math:`\eta_B`.

The Resta Screened Dielectric Kernel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In semiconductor nanoclusters, dielectric screening varies continuously from the bulk optical dielectric constant :math:`\epsilon_\infty` at large distances down to unscreened vacuum interaction (:math:`\epsilon = 1`) at short, intra-atomic distances.

``QDEX`` implements the **Resta model of electronic screening**:

.. math::

   W(r) = \frac{1}{\epsilon_\infty r} + \frac{1 - \epsilon_\infty^{-1}}{r} \exp\left( -\frac{r}{\lambda_s} \right)

where :math:`\lambda_s` is the **Thomas-Fermi screening length** of the valence electron gas:

.. math::

   \lambda_s = \sqrt{ \frac{\pi}{4 k_F} } = \left( \frac{\pi}{4 (3\pi^2 n_v)^{1/3}} \right)^{1/2}

Here, :math:`n_v` is the valence electron density of the semiconductor lattice.

Damped over atomic centers with Ohno-Klopman hardness parameters, the discrete kernel is:

.. math::

   W_{AB}^{\mathrm{Resta}} = \frac{1}{\epsilon_\infty \sqrt{R_{AB}^2 + a_{AB}^2}} + \frac{1 - \epsilon_\infty^{-1}}{\sqrt{R_{AB}^2 + a_{AB}^2}} \exp\left( -\frac{\sqrt{R_{AB}^2 + a_{AB}^2}}{\lambda_s} \right)

Why This Approximation Works for Quantum Dots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Avoids Over-Screening at Atomic Distances**: At short range (:math:`r \ll \lambda_s`), :math:`W(r) \to 1/r`. Conventional constant-dielectric models (:math:`W(r) = 1 / (\epsilon_\infty r)`) artificially divide the Coulomb attraction by :math:`\epsilon_\infty \approx 5 - 10` even when the electron and hole reside on the same atom, underestimating the exciton binding energy by hundreds of meV.
2. **Captures Macroscopic Screening**: At large distances (:math:`r \gg \lambda_s`), :math:`W(r) \to 1 / (\epsilon_\infty r)`, correctly reproducing the bulk dielectric polarization of the nanocrystal core.
3. **Atomic Charge Centroids**: In semiconductors, valence and conduction band states are predominantly composed of atomic :math:`s` and :math:`p` orbitals. Mulliken transition charges accurately preserve the spatial monopole centroids of the electron-hole charge distribution.

---

4. Physical Rationale for Diagonal BSE
--------------------------------------

In standard BSE, diagonalizing the full :math:`A_{ia, jb}` matrix requires :math:`O(N_{\mathrm{pairs}}^3)` operations. In a quantum dot with 50 occupied and 50 virtual orbitals, :math:`N_{\mathrm{pairs}} = 2,500`, which is manageable. But in larger clusters with 200 occupied and 200 virtual states, :math:`N_{\mathrm{pairs}} = 40,000`, requiring gigabytes of memory and prohibitive diagonalization times.

The **Diagonal BSE** framework omits off-diagonal configuration interaction (:math:`ia \neq jb`), evaluating:

.. math::

   \Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} + 2 K_{ia, ia}^x - K_{ia, ia}^d

Why Diagonal BSE Works in Nanocrystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Dominance of the Diagonal Coulomb Attraction**: The diagonal kernel element :math:`K_{ia, ia}^d` represents the direct electrostatic attraction between the electron distribution :math:`|\phi_a|^2` and the hole distribution :math:`|\phi_i|^2`. In quantum dots, this diagonal term accounts for **over 90% of the total exciton binding energy** :math:`E_b`.
* **High Density of States**: In dense excitonic manifolds, off-diagonal mixing redistributes intensity among neighboring quasi-degenerate pairs without substantially altering the global optical absorption envelope.
* **Essential for Molecular Dynamics (NAMD)**: In non-adiabatic molecular dynamics simulations where excited states must be evaluated at every time step (e.g. 5,000 steps), full BSE diagonalization is computationally impossible. Diagonal BSE provides an accurate, energy-conserving potential energy surface at a fraction of the cost.

---

5. Full BSE (sTDA) & The Davidson Iterative Solver
--------------------------------------------------

When accurate oscillator strength borrowing, state mixing, or fine excitonic splittings are required, ``QDEX`` employs the full **sTDA** framework.

Transition Dipole Redistribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The single-particle transition dipole moment between occupied orbital :math:`i` and virtual orbital :math:`a` is:

.. math::

   \boldsymbol{\mu}_{ia} = \langle \phi_i | \mathbf{r} | \phi_a \rangle = \sum_{\mu \nu} C_{\mu i} \mathbf{D}_{\mu \nu} C_{\nu a}

where :math:`\mathbf{D}_{\mu \nu} = \langle \chi_\mu | \mathbf{r} | \chi_\nu \rangle` is the AO dipole matrix evaluated via Libint2.

Under full configuration interaction, the exciton transition dipole moment is a coherent linear superposition:

.. math::

   \boldsymbol{\mu}_S = \sum_{ia} X_{ia}^S \, \boldsymbol{\mu}_{ia}

The corresponding dimensionless oscillator strength is:

.. math::

   f_S = \frac{2}{3} \, \Omega_S \, |\boldsymbol{\mu}_S|^2

This coherent summation describes **superradiance** and intensity borrowing, where optical strength from high-energy transitions is transferred into the lowest bright exciton.

The Davidson Eigensolver
~~~~~~~~~~~~~~~~~~~~~~~~

To solve for the lowest :math:`k` roots without dense :math:`O(N_{\mathrm{pairs}}^3)` matrix diagonalization, ``QDEX`` implements an optimized **Davidson iterative subspace solver**:

1. Projects :math:`\mathbf{A}` into a small trial subspace :math:`\mathbf{V} = [\mathbf{v}_1, \dots, \mathbf{v}_m]`.
2. Computes the matrix-vector product :math:`\mathbf{w}_j = \mathbf{A} \mathbf{v}_j` on-the-fly.
3. Solves the projected eigenvalue problem :math:`\mathbf{V}^\dagger \mathbf{A} \mathbf{V} \mathbf{y} = \omega \mathbf{y}`.
4. Computes the residual :math:`\mathbf{r} = \mathbf{A} \mathbf{x} - \omega \mathbf{x}` and preconditioner :math:`\mathbf{\delta} = (\operatorname{diag}(\mathbf{A}) - \omega)^{-1} \mathbf{r}`.
5. Expands the subspace until the norm :math:`||\mathbf{r}|| < \text{tol}` (default :math:`10^{-5}`).

This reduces memory requirements from :math:`O(N_{\mathrm{pairs}}^2)` to :math:`O(N_{\mathrm{pairs}} \times k)`.

---

6. CLI Flags & YAML Configuration Reference
-------------------------------------------

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
     - Excitation model: ``bse`` (full sTDA), ``diagonal_bse``, ``independent_qp``, or ``independent_dft``.
   * - ``--kernel <choice>``
     - ``bse``
     - Dielectric screening kernel: ``bse`` (legacy) or ``resta`` (electronic screening model).
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

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   physics:
     excitation_mode: "diagonal_bse"
     kernel: "resta"
     charge_type: "mulliken"
     nhomos: 50
     nlumos: 50
     triplet: false

   bse:
     nroots: 20
     full_diag: false
     tol: 1.0e-5
