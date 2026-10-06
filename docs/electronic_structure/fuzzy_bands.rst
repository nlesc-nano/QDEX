Fuzzy bands
===========

Part of :doc:`/electronic_structure/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.pdos_coop``
* Callable: ``qdex.pdos_coop.compute_pdos_and_coop``
* CLI: ``--charge_type, --run_fuzzy``
* YAML: ``integrals.charges``, ``analysis.run_fuzzy``

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

Folding is needed whenever a band edge has little weight at :math:`\mathbf{G} = 0`.
The CdSe valence-band maximum at :math:`\Gamma` is Se 4p-like: by symmetry its
:math:`\mathbf{G} = 0` component vanishes, so without folding the top of the valence
band is invisible at :math:`\Gamma`. Convergence for Cd16Se13Cl6 and
Cd68Se55Cl26 (PBE/DZVP, zinc-blende k-path, relative to ``g_shell: 3`` = 343
replicas; the map is the energy-smeared spin-free fuzzy map):

.. list-table::
   :header-rows: 1

   * - setting
     - map difference
     - VBM weight
     - VBM weight at Γ
     - CBM weight
   * - no folding
     - 80-83 %
     - 5-11 %
     - 0
     - 72 %
   * - ``g_shell: 1``
     - 7 %
     - 79-86 %
     - 76-83 %
     - 94-97 %
   * - ``g_shell: 2``
     - 1 %
     - 98 %
     - 97-99 %
     - 99.5 %

``fold_to_bz: true`` with ``g_shell: 2`` is the default. The cost grows
linearly with the number of replicas and stays small next to the rest of a run
(two seconds more than ``g_shell: 1`` on Cd68). The plane-wave transforms are
computed one replica at a time, so memory does not grow with ``g_shell``.


Automated High-Symmetry Paths & Lattice Orientation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given a reference crystal structure (``.cif`` file), ``QDEX``:

1. Determines the space group and high-symmetry :math:`k`-path using ``pymatgen`` (e.g., :math:`\Gamma \to X \to W \to K \to \Gamma \to L`).
2. Fits the orientation of the dot's lattice: the nearest-neighbour bond directions of the interior atoms are matched to the bond star of the crystal (Kabsch fit, started from every pair of crystal bonds whose angle matches two bonds of the most central atom). The rotation is fixed up to the symmetry of the bond star, which leaves the fuzzy weights unchanged (crystal point group plus :math:`\mathbf{k} \to -\mathbf{k}`). The log reports the mean cosine of the fit; below 0.95 the dot has no well-ordered core of the CIF structure and the k-directions are uncertain.
3. Scales the reciprocal lattice to the median interior bond length of the dot.

The k-path and the :math:`\mathbf{G}` replicas are expressed in the dot's frame, so the
folding uses the reciprocal lattice of the dot itself.


Bulk Band Overlay
~~~~~~~~~~~~~~~~~

The dashboards draw the bulk band structure of the material (``qdex/data/bulk_bands``,
CP2K PBE with the same basis and pseudopotentials) over the fuzzy map. The bulk segments
are placed on the fuzzy path by their k-coordinates.

The bulk bands are put on the dot's energy axis with a semicore level (from ``<name>.json``:
the cation d band, Cs 5p in the perovskites, otherwise the anion s band). For CdSe the Cd 4d
level of the dot is the Mulliken-weighted 4d energy of its *interior, bulk-like* Cd atoms (four
Se neighbours, no ligand, inner half by radius), and the bulk Cd 4d bands are placed there.
Cd bonded to Cl has its 4d level about 0.4 eV deeper, so averaging over all Cd atoms would
put the bulk bands too low. With the anchor in place the remaining offsets are physical:
confinement pushes the dot states away from the band extrema (down at the valence band
maximum, up at a valence band minimum such as the bottom of the p band at L).

The SOC dashboard uses spin-orbit bulk bands (``<name>_soc.bs.gz``) when they exist,
so that the :math:`\Gamma_8/\Gamma_7` splitting and the split bands at :math:`L` are in
the overlay as well. They are computed with the same GTH-SOC operator as the dots, from a
CP2K k-point run that prints the real-space Kohn-Sham and overlap matrices:

.. code-block:: text

   &PRINT
     &KS_CSR_WRITE
       REAL_SPACE .TRUE.
       UPPER_TRIANGULAR .FALSE.
     &END KS_CSR_WRITE
     &S_CSR_WRITE
       REAL_SPACE .TRUE.
       UPPER_TRIANGULAR .FALSE.
     &END S_CSR_WRITE
     &TREXIO
     &END TREXIO
     &BAND_STRUCTURE
       ...                      # the k-path of the fuzzy bands
     &END BAND_STRUCTURE
   &END PRINT

.. code-block:: bash

   python -m qdex.bulk_soc RUN_DIR --basis BASIS_MOLOPT_UZH --gth GTH_SOC_POTENTIALS \
       --name CdSe_zb --cif CdSe_zb.cif -d qdex/data/bulk_bands --gzip

The bulk bands of all the materials of the QDSpaceWebApp are listed in :doc:`bulk_bands`.

:math:`H(\mathbf{k}) = \sum_\mathbf{R} e^{i\mathbf{k}\cdot\mathbf{R}} \langle\chi_\mu(0)|H|\chi_\nu(\mathbf{R})\rangle`
and :math:`S(\mathbf{k})` are built at every point of the path, the SOC operator from
Bloch sums of the projector overlaps, and the two-component problem is solved in the full
AO basis. The run checks the AO basis against CP2K's overlap and the spin-free bands
against CP2K's band structure. For zinc-blende CdSe (PBE, DZVP-MOLOPT-PBE-GTH): the
spin-free bands agree with CP2K to 0.06 meV and :math:`\Delta_{so}(\Gamma) = 0.36` eV.
