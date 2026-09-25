Architecture
============

Part of :doc:`/interactions/index`.

.. important::

   The ``--2e-integrals`` axis chooses the representation of the bare interaction; ``--kernel`` chooses a screening builder. Some kernel names force AO or atom resolution in ``qdex.solver.ExcitonSolver``. Inspect that dispatcher for supported combinations.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.build_gamma``
* CLI: ``--2e-integrals, --kernel, --eps-out``
* YAML: ``physics.2e-integrals, physics.kernel, physics.eps_out``

.. code-block:: python

   build_gamma(atom_symbols, coords, alpha, beta=0.0, eta_dict=HARDNESS_DICT)


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
     - **Bare Coulomb Representation**: Chooses how :math:`v(\mathbf{r}_1, \mathbf{r}_2)` is evaluated (semi-empirical atom-centered point charges vs. analytical Gaussian AO density-pair integrals).
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
   * - **Short-Range Limit** (:math:`r \to 0`)
     - Resta/DIM set :math:`W_{AA}=\gamma_{AA}` (unscreened, MNOK-damped); sBSE may screen on-site terms. The on-site conventions are not uniform across kernels.
     - Preserves full on-site Coulomb repulsion and atomic exchange splitting.
   * - **Bulk Limit** (:math:`R_{\mathrm{QD}} \to \infty`)
     - :math:`\Delta W \to 0` only if the QD screening model converges to the same bulk reference; the anchor model reaches :math:`\Delta_{\mathrm{bulk}}` by construction.
     - The long-range tail tends to :math:`-e^2/(\epsilon_\infty r)`. This is necessary but not sufficient for the Wannier–Mott limit, which also needs converged active spaces, band-edge masses and the full coupled BSE (not tested; see :doc:`/validation/bulk_exciton_limit`).


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

