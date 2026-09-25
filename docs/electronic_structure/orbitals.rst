Orbitals
========

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

.. rubric:: From ``docs/part1_ground_state/index.rst:1-9``


The starting point of any ``QDEX`` calculation is the ground-state electronic structure of the nanocluster or quantum dot (QD), obtained from a Density Functional Theory (DFT) calculation performed with **CP2K / Quickstep**.

This section details the theoretical foundation of molecular orbital expansion, population analysis, Projected Density of States (PDOS), Inverse Participation Ratio (IPR), Crystal Orbital Overlap Population (COOP), supercell unfolding (Fuzzy Bands), and 3D volumetric orbital visualization via Gaussian ``.cube`` files.

---


.. rubric:: From ``docs/part1_ground_state/index.rst:10-31``

1. Molecular Orbitals and the Atomic Orbital Basis
--------------------------------------------------

In CP2K's Gaussian and Plane Waves (GPW) formalism, the one-electron Kohn-Sham molecular orbitals (MOs) :math:`\phi_m(\mathbf{r})` are expanded in a linear combination of atom-centered Gaussian-type basis functions (AOs) :math:`\chi_\mu(\mathbf{r})`:

.. math::

   \phi_m(\mathbf{r}) = \sum_{\mu=1}^{N_{\mathrm{ao}}} C_{\mu m} \chi_\mu(\mathbf{r})

where:
* :math:`N_{\mathrm{ao}}` is the total number of atomic orbitals in the spherical representation.
* :math:`C_{\mu m}` is the molecular orbital coefficient matrix element representing the contribution of AO :math:`\mu` to MO :math:`m`.
* :math:`\chi_\mu(\mathbf{r}) = R_{nl}(|\mathbf{r} - \mathbf{R}_I|) Y_{lm}(\theta, \phi)` are atom-centered Gaussian basis functions characterized by atomic center :math:`\mathbf{R}_I`, principal shell :math:`n`, and spherical harmonics :math:`Y_{lm}`.

The overlap between non-orthogonal atomic orbitals is defined by the symmetric overlap matrix :math:`\mathbf{S}`:

.. math::

   S_{\mu \nu} = \langle \chi_\mu | \chi_\nu \rangle = \int \chi_\mu^*(\mathbf{r}) \chi_\nu(\mathbf{r}) \, d\mathbf{r}

In ``QDEX``, :math:`\mathbf{S}` is computed analytically using the high-performance C++ backend powered by **Libint2**.


.. rubric:: From ``docs/part1_ground_state/index.rst:32-50``

MO Orthonormality and Diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The set of Kohn-Sham molecular orbitals must satisfy strict orthonormality:

.. math::

   \langle \phi_m | \phi_n \rangle = \sum_{\mu \nu} C_{\mu m}^* S_{\mu \nu} C_{\nu n} = (\mathbf{C}^\dagger \mathbf{S} \mathbf{C})_{mn} = \delta_{mn}

Before performing any further analysis or excited-state calculations, ``QDEX`` automatically computes the matrix product :math:`\mathbf{S} \mathbf{C}` and checks the Frobenius norm and maximum element-wise deviation from identity:

.. math::

   \Delta_{\mathrm{orth}} = \max_{\mu, \nu} |(\mathbf{C}^\dagger \mathbf{S} \mathbf{C})_{\mu \nu} - \delta_{\mu \nu}|

If :math:`\Delta_{\mathrm{orth}} > 10^{-5}`, the run is halted to avoid propagating numerical inconsistencies arising from mismatched basis set orders or spherical-to-Cartesian convention discrepancies.

---
