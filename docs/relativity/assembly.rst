Assembly
========

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

.. rubric:: From ``docs/part2_soc/index.rst:187-204``

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
