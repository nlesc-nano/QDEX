Model selection
===============

Part of :doc:`/quasiparticles/index`.

.. important::

   The labels ``sgw``, ``evgw`` and ``qsgw`` in QDEX name reduced models. They should not be interpreted as a claim of numerical equivalence to conventional GW/QSGW calculations. Validate each against a common reference set.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:165-188``

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
   * - **Microscopic Dielectric Shift** (:math:`\Delta W`)
     - ``sgw-dim``, ``sgw-resta``, ``sgw``
     - Microscopic screened Coulomb difference :math:`\Delta W = W^{\mathrm{QD}} - W^{\mathrm{bulk}} + W^{\mathrm{solv}}`. Local :math:`v_{xc}` differences are neglected as a model approximation; bulk reference parameters remain required.
   * - **Self-Consistent Extensions**
     - ``evgw-dim``, ``evgw-resta``, ``qsgw-dim``, ``qsgw-resta``
     - Iterative eigenvalue self-consistency (:math:`evGW`) and static AO-basis orbital relaxation (``qsgw-*``) with empirical damping factor :math:`Z_p`.

---


.. rubric:: From ``docs/part3_gw_scissor/index.rst:1001-1062``

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
