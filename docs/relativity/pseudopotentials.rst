Pseudopotentials
================

Part of :doc:`/relativity/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.soc_utils``
* Callable: ``qdex.soc_utils.compute_spinor_subspace``
* CLI: ``--soc_flag, --gth_file``
* YAML: ``physics.soc, physics.soc_window_ev``

.. code-block:: python

   compute_spinor_subspace(atom_symbols, coords_ang, shells, C_AO, eps_Ha, S_AO, active_indices, gth_file, nthreads=1, soc_cache=None, assume_orthonormal=False, SC_AO=None, device='numpy', verbose=True)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part2_soc/index.rst:76-100``

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


.. rubric:: From ``docs/part2_soc/index.rst:101-110``

Why QDEX Chooses GTH SOC Pseudopotentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``QDEX`` is tailored for post-processing CP2K calculations. In CP2K Quickstep, the core electrons are integrated out using norm-conserving GTH pseudopotentials. By evaluating SOC within the exact same **separable GTH projector framework**:
1. **Rigorous Consistency**: The relativistic core potentials match the underlying Kohn-Sham ground state without parameter re-fitting.
2. **Semi-Local Projector Form**: The SOC Hamiltonian is completely factorized into separable inner products with atom-centered Gaussian projector functions.
3. **Analytic C++ Evaluation**: Matrix elements of the Gaussian projectors are evaluated analytically via Libint2, eliminating numerical integration grids.

---
