Overview
========

Part of :doc:`/workflows/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system.*, physics.*, namd.*``

.. code-block:: python

   main()

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part7_examples/index.rst:1-14``


This tutorial series walks you through the entire ``QDEX`` workflow using a representative :math:`\text{CsPbBr}_3` perovskite quantum dot. The tutorials mirror the exact progressive curriculum from Parts 1 through 6:

1. Ground-state electronic structure (PDOS, IPR, COOP, Fuzzy Bands, Cube orbitals).
2. Relativistic Spin-Orbit Coupling (SOC).
3. Scaled GW Quasiparticle corrections and absolute band edge prediction.
4. Optical excitations across the four frameworks (DFT, QP, Diagonal BSE, and sTDA).
5. Wavefunction analysis (Plasser-Dreuw metrics and Natural Transition Orbitals).
6. Non-adiabatic carrier cooling dynamics, emission rates, and PLQY.

---
