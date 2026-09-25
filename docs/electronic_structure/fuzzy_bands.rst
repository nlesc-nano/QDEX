Fuzzy bands
===========

Part of :doc:`/electronic_structure/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.pdos_coop``
* Callable: ``qdex.pdos_coop.compute_pdos_and_coop``
* CLI: ``--charge_type, --run_fuzzy``
* YAML: ``physics.charge_type, fuzzy.run``

.. code-block:: python

   compute_pdos_and_coop(C, S, eps_eV, shells, pdos_atoms, coop_pairs, ewin, sigma=0.03, is_soc=False, prefix='sf', pops=None, population_bars=None, device='numpy')


5. Fuzzy Band Structure (Supercell Unfolding)
---------------------------------------------

Because quantum dots and nanocrystals are finite non-periodic clusters, they have discrete energy levels rather than continuous dispersion relations :math:`E(\mathbf{k})`. However, researchers frequently need to compare nanocrystal states against the parent bulk band structure (e.g., to observe quantum confinement shifts at the :math:`R`- or :math:`\Gamma`-point).

``QDEX`` implements the **Fuzzy Band** plane-wave projection algorithm to unfold cluster molecular orbitals onto effective bulk crystal wavevectors :math:`\mathbf{k}`.


Plane-Wave Fourier Projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each cluster molecular orbital :math:`\phi_m(\mathbf{r})` is projected onto a plane wave :math:`|\mathbf{k}\rangle = \frac{1}{\sqrt{V}} e^{i \mathbf{k} \cdot \mathbf{r}}`:

.. math::

   F_m(\mathbf{k}) = \langle e^{i \mathbf{k} \cdot \mathbf{r}} | \phi_m \rangle = \int e^{-i \mathbf{k} \cdot \mathbf{r}} \phi_m(\mathbf{r}) \, d\mathbf{r} = \sum_{\mu=1}^{N_{\mathrm{ao}}} C_{\mu m}^* \int e^{-i \mathbf{k} \cdot \mathbf{r}} \chi_\mu(\mathbf{r}) \, d\mathbf{r}

The Fourier transform of the Gaussian basis functions :math:`F_\mu(\mathbf{k}) = \int e^{-i \mathbf{k} \cdot \mathbf{r}} \chi_\mu(\mathbf{r}) \, d\mathbf{r}` is evaluated in closed analytical form via ``libint_cpp.ao_ft_complex``.

The spectral weight (fuzzy intensity) of state :math:`m` at wavevector :math:`\mathbf{k}` is:

.. math::

   I_m(\mathbf{k}) = |F_m(\mathbf{k})|^2 = \left| \sum_\mu C_{\mu m}^* F_\mu(\mathbf{k}) \right|^2


Brillouin Zone Folding
~~~~~~~~~~~~~~~~~~~~~~

For a finite nanocluster, momentum conservation is relaxed, spreading the spectral weight across reciprocal space. To reconstruct an effective band structure within the first Brillouin zone, ``QDEX`` allows summing over reciprocal lattice vectors :math:`\mathbf{G}` of the reference bulk cell:

.. math::

   I_m^{\mathrm{folded}}(\mathbf{k}) = \sum_{\mathbf{G} \in \mathcal{S}_g} |F_m(\mathbf{k} + \mathbf{G})|^2

where :math:`\mathcal{S}_g` denotes reciprocal shells controlled by the keyword ``g_shell``:

* ``g_shell: 0``: 1 replica (:math:`\mathbf{G} = \mathbf{0}`).
* ``g_shell: 1``: 27 reciprocal lattice replicas (:math:`h, k, l \in \{-1, 0, +1\}`).
* ``g_shell: 2``: 125 reciprocal lattice replicas.


Automated High-Symmetry Paths & PCA Alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given a reference crystal structure (``.cif`` file), ``QDEX``:

1. Determines the space group and high-symmetry :math:`k`-path using ``pymatgen`` (e.g., :math:`\Gamma \to X \to M \to \Gamma \to R`).
2. Scales the reciprocal lattice to match the core bond distances of the relaxed quantum dot.
3. Performs Principal Component Analysis (PCA) on the inertia tensors of the CIF and cluster geometries to automatically align rotational coordinate axes.

