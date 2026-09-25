Foundations
===========

Part of :doc:`/excitons/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.solver``
* Callable: ``qdex.solver.solve``
* CLI: ``--excitation-mode, --include-direct-eh, --include-exchange``
* YAML: ``physics.excitation_mode, physics.include_direct_eh, physics.include_exchange``

.. code-block:: python

   solve(self, nroots=10, full_diag=False, tol=1e-05, excitation_mode='bse')

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part4_excited_states/index.rst:1-9``


The description of neutral optical excitations in semiconductor nanostructures requires treating the two-particle correlated motion of an electron promoted to the conduction band and the hole left behind in the valence band.

``QDEX`` provides four distinct theoretical frameworks for computing excited states, ranging from non-interacting single-particle transitions to the fully coupled **Bethe-Salpeter Equation (BSE)** under the Tamm-Dancoff Approximation (TDA), coupled with a clear separation between **Two-Electron Integral Representations** and **Dielectric Screening Kernels**.

---


.. rubric:: From ``docs/part4_excited_states/index.rst:10-22``

1. The Two-Particle Excitation Problem
--------------------------------------

In a neutral optical excitation, an electron is removed from an occupied valence orbital :math:`i` and placed into an unoccupied conduction orbital :math:`a`, creating an electron-hole pair configuration :math:`|ia\rangle = a_a^\dagger a_i |0\rangle`.

The exact correlated excited-state wavefunction :math:`|\Psi_S\rangle` is expressed as a linear superposition of electron-hole configurations:

.. math::

   |\Psi_S\rangle = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, |ia\rangle

where :math:`X_{ia}^S` are the configuration interaction amplitudes and :math:`\Omega_S` is the corresponding optical excitation energy.


.. rubric:: From ``docs/part4_excited_states/index.rst:23-51``

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


.. rubric:: From ``docs/part4_excited_states/index.rst:52-66``

Singlet vs. Triplet Matrix Elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The elements of the resonant matrix :math:`A_{ia, jb}` depend on the spin multiplicity:

.. math::

   A_{ia, jb}^{\mathrm{singlet}} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} + 2 K_{ia, jb}^x - K_{ia, jb}^d

.. math::

   A_{ia, jb}^{\mathrm{triplet}} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} - K_{ia, jb}^d.

The bare exchange term :math:`K_{ia, jb}^x` is strictly absent in triplet states because electrons with parallel spins experience identical spatial exchange cancellation. The factor of :math:`2 K_{ia, jb}^x` in singlets is responsible for the singlet-triplet exchange splitting.


.. rubric:: From ``docs/part4_excited_states/index.rst:67-79``

Relativistic 2-Component Spinor BSE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When spin-orbit coupling is enabled (``soc: true``), the single-particle spatial orbitals are expanded into 2-component Kramers spinors :math:`\psi_p(\mathbf{r}) = \begin{pmatrix} \phi_{p\alpha}(\mathbf{r}) \\ \phi_{p\beta}(\mathbf{r}) \end{pmatrix}`. The BSE Hamiltonian is transformed into the spinor electron-hole basis :math:`|IA\rangle = a_A^\dagger a_I |0\rangle`:

.. math::

   A_{IA, JB}^{\mathrm{spinor}} = \left( \varepsilon_A^{\mathrm{QP}} - \varepsilon_I^{\mathrm{QP}} \right) \delta_{IJ} \delta_{AB} + K_{IA, JB}^x - K_{IA, JB}^d

capturing fine-structure splittings, bright-dark exciton order inversion, and Rashba effects without empirical parameters.

---
