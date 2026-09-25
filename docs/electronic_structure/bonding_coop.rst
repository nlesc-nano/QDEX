Bonding coop
============

Part of :doc:`/electronic_structure/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.pdos_coop``
* Callable: ``qdex.pdos_coop.compute_pdos_and_coop``
* CLI: ``--charge_type, --run_fuzzy``
* YAML: ``physics.charge_type, fuzzy.run``

.. code-block:: python

   compute_pdos_and_coop(C, S, eps_eV, shells, pdos_atoms, coop_pairs, ewin, sigma=0.03, is_soc=False, prefix='sf', pops=None, population_bars=None, device='numpy')

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part1_ground_state/index.rst:140-144``

4. Crystal Orbital Overlap Population (COOP)
--------------------------------------------

The **Crystal Orbital Overlap Population (COOP)** measures the nature and strength of chemical bonding between selected pairs of atoms :math:`A` and :math:`B` (e.g. Pb–Br or Cs–Br bonds) as a function of energy.


.. rubric:: From ``docs/part1_ground_state/index.rst:145-159``

Mathematical Definition
~~~~~~~~~~~~~~~~~~~~~~~

The COOP contribution of molecular orbital :math:`m` for atom pair :math:`(A, B)` is given by:

.. math::

   \mathrm{COOP}_{AB}(m) = 2 \sum_{\mu \in A} \sum_{\nu \in B} \operatorname{Re}\left[ C_{\mu m}^* S_{\mu \nu} C_{\nu m} \right]

For two-component relativistic spinors :math:`|\psi_k\rangle = \begin{pmatrix} \mathbf{C}_k^\alpha \\ \mathbf{C}_k^\beta \end{pmatrix}`, the COOP sums over both spin channels:

.. math::

   \mathrm{COOP}_{AB}^{\mathrm{SOC}}(k) = 2 \sum_{\mu \in A} \sum_{\nu \in B} \operatorname{Re}\left[ C_{\mu k}^{\alpha *} S_{\mu \nu} C_{\nu k}^\alpha + C_{\mu k}^{\beta *} S_{\mu \nu} C_{\nu k}^\beta \right]


.. rubric:: From ``docs/part1_ground_state/index.rst:160-170``

Bonding vs. Antibonding Character
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* :math:`\mathrm{COOP}_{AB} > 0`: **Bonding interaction** (constructive orbital interference, electronic charge accumulates in the internuclear region).
* :math:`\mathrm{COOP}_{AB} < 0`: **Antibonding interaction** (destructive interference, nodal plane between the nuclei).
* :math:`\mathrm{COOP}_{AB} \approx 0`: **Non-bonding interaction**.

In lead halide perovskites, the valence band maximum consists of an antibonding mixture between Pb :math:`6s` and halide :math:`np` orbitals, explaining the renowned defect tolerance of perovskites (since vacancies remove antibonding states).

---
