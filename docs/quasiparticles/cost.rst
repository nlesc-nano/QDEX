Cost
====

Part of :doc:`/quasiparticles/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:154-164``

4. The Nanocrystal Scaling Bottleneck
-------------------------------------

Although full :math:`G_0W_0` calculations resolve the band gap problem for small molecules or bulk crystals, computing :math:`G_0W_0` explicitly for colloidal quantum dots is computationally prohibitive:

* **Computational Complexity**: Standard :math:`G_0W_0` scales as :math:`O(N^4)` with real-space basis sets and :math:`O(N^5)` with plane-wave expansions.
* **Memory Footprint**: Calculating the polarizability matrix :math:`\chi_0(\mathbf{r}, \mathbf{r}'; \omega)` requires summing over thousands of unoccupied conduction states, demanding terabytes of RAM for clusters with :math:`> 500` atoms.
* **Nanocrystal Realities**: Chemically realistic colloidal quantum dots comprise 1,000 to 10,000 atoms, including surface passivation ligands, rendering direct *ab initio* :math:`G_0W_0` impossible for high-throughput screening or non-adiabatic molecular dynamics trajectories.

---
