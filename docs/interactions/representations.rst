Representations
===============

Part of :doc:`/interactions/index`.

.. important::

   The ``xs`` backend evaluates analytical Gaussian integrals of the restricted form ``(mu mu | nu nu)``. It does not store the full four-index electron-repulsion tensor; the MO interactions are density-pair approximations. The MNOK heteronuclear damping currently uses the mean of inverse hardnesses, not ``2/(eta_A+eta_B)``.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.integrals``
* Callable: ``qdex.integrals.compute_two_electron_ao``
* CLI: ``--2e-integrals, --kernel, --eps-out``
* YAML: ``physics.2e-integrals, physics.kernel, physics.eps_out``

.. code-block:: python

   compute_two_electron_ao(shells, nthreads=1)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part4_excited_states/index.rst:143-145``

3. Two-Electron Integral Representations (``2e-integrals``)
-----------------------------------------------------------


.. rubric:: From ``docs/part4_excited_states/index.rst:146-179``

Semi-Empirical Atom-Centered Representation (``2e-integrals: mnok``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To achieve high efficiency for nanocrystals containing up to 10,000 atoms, ``QDEX`` contracts transition densities into atom-centered point charges using Mulliken or Lowdin population analysis:

.. math::

   q_A^{ia} = \sum_{\mu \in A} \sum_{\nu=1}^{N_{\mathrm{ao}}} C_{\mu i} S_{\mu \nu} C_{\nu a}.

For diagonal electron and hole charge densities:

* Hole density on atom :math:`A`: :math:`q_A^{ii} = \sum_{\mu \in A} \sum_\nu C_{\mu i} S_{\mu \nu} C_{\nu i}`
* Electron density on atom :math:`B`: :math:`q_B^{aa} = \sum_{\mu \in B} \sum_\nu C_{\mu a} S_{\mu \nu} C_{\nu a}`

The four-center integrals are replaced by pairwise contractions over atomic sites:

.. math::

   K_{ia, jb}^x = \sum_{A=1}^{N_{\mathrm{atoms}}} \sum_{B=1}^{N_{\mathrm{atoms}}} q_A^{ia} \, \gamma_{AB}^{\mathrm{Ohno}} \, q_B^{jb}

.. math::

   K_{ia, jb}^d = \sum_{A=1}^{N_{\mathrm{atoms}}} \sum_{B=1}^{N_{\mathrm{atoms}}} q_A^{ij} \, W_{AB} \, q_B^{ab}

where :math:`\gamma_{AB}^{\mathrm{Ohno}}` is the Mataga-Nishimoto-Ohno-Klopman (MNOK) damped Coulomb potential:

.. math::

   \gamma_{AB}^{\mathrm{Ohno}} = \frac{1}{\sqrt{R_{AB}^2 + a_{AB}^2}}

and :math:`a_{AB} = 2 / (\eta_A + \eta_B)` is the Ohno-Klopman damping parameter derived from atomic chemical hardnesses :math:`\eta_A` and :math:`\eta_B`.

**Advantages**: The atom-pair kernel stores :math:`O(N_{\mathrm{atoms}}^2)` elements; actual memory depends on system size and solver intermediates.


.. rubric:: From ``docs/part4_excited_states/index.rst:180-206``

Analytical AO Density-Pair Representation (``2e-integrals: xs``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For smaller clusters (:math:`\le 500` atoms) or benchmark comparisons, ``QDEX`` evaluates the **analytical AO density-pair Coulomb integrals** analytically via the C++ Libint2 library:

.. math::

   (\mu \nu | \lambda \sigma) = \iint \chi_\mu(\mathbf{r}_1) \chi_\nu(\mathbf{r}_1) \, \frac{1}{|\mathbf{r}_1 - \mathbf{r}_2|} \, \chi_\lambda(\mathbf{r}_2) \chi_\sigma(\mathbf{r}_2) \, d\mathbf{r}_1 \, d\mathbf{r}_2.

The two-electron atomic repulsion matrix is computed over AO density pairs:

.. math::

   \Gamma_{\mu \nu} = (\mu \mu | \nu \nu).

The molecular orbital matrix elements are constructed by density-pair AO-to-MO contraction:

.. math::

   K_{ia, jb}^x = \sum_{\mu, \nu=1}^{N_{\mathrm{ao}}} \left( C_{\mu i} C_{\nu a} \right) \Gamma_{\mu \nu} \left( C_{\mu j} C_{\nu b} \right)

.. math::

   K_{ia, jb}^d = \sum_{\mu, \nu=1}^{N_{\mathrm{ao}}} \left( C_{\mu i} C_{\nu j} \right) W_{\mu \nu} \left( C_{\mu a} C_{\nu b} \right).

**Advantages**: Preserves non-spherical orbital angular momentum components (e.g. anisotropic :math:`p` and :math:`d` orbital bonding interactions), eliminating any reliance on spherical atomic charge partitioning.


.. rubric:: From ``docs/part4_excited_states/index.rst:207-224``

Implementation in QDEX (``2e-integrals``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-electron integral representations are managed across :mod:`qdex.hardness`, :mod:`qdex.integrals`, and :mod:`qdex.exciton_hamiltonian`:

1. **Semi-Empirical Atom-Centered Points (``2e-integrals: mnok``)**:
   - Implemented via :func:`qdex.hardness.build_gamma` and :class:`qdex.exciton_hamiltonian.ExcitonHamiltonian`.
   - Projects molecular orbitals onto atomic centers using Löwdin or Mulliken populations: :math:`\mathbf{q}_A^{ia} = \sum_{\mu \in A, \nu} C_{\mu i} S_{\mu \nu} C_{\nu a}` via optimized BLAS ``DGEMM`` routines.
   - Pairs transition charges with the damped Ohno matrix :math:`\mathbf{\Gamma}_{AB}^{\mathrm{Ohno}}` for exchange and :math:`\mathbf{W}_{AB}` for direct screening.
   - Fast, memory-lean (< 10 MB RAM), and scales easily to nanocrystals containing 1,000 to 10,000 atoms.

2. **Analytical 4-Center GTO Integrals (``2e-integrals: xs``)**:
   - Evaluated via the C++ extension :mod:`libint_cpp` and wrapped in :func:`qdex.integrals.compute_two_electron_ao`.
   - Calculates exact two-electron Gaussian repulsion integrals :math:`(\mu \mu | \nu \nu)` over contracted GTO basis shells.
   - Preserves complete angular orbital anisotropy without spherical approximations. Recommended for molecular benchmarks and small nanoclusters (:math:`\le 500` atoms).

---
