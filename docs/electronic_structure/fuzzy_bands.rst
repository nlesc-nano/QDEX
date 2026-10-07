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
2. Fits the orientation of the dot's lattice: the nearest-neighbour bond directions of the interior atoms (only between element pairs bonded in the crystal, so Cs–X contacts of a perovskite are left out) are matched to the bond star of the crystal (Kabsch fit, started from every pair of crystal bonds whose angle matches two bonds of the most central atom). The rotation is fixed up to the symmetry of the bond star, which leaves the fuzzy weights unchanged (crystal point group plus :math:`\mathbf{k} \to -\mathbf{k}`). The log reports the mean cosine of the fit; below 0.95 the dot has no well-ordered core of the CIF structure and the k-directions are uncertain.
3. Scales the reciprocal lattice to the dot's lattice: the median nearest same-element distance of the interior atoms of the most coordinated element (Pb–Pb, the pseudo-cubic lattice constant, in a perovskite; Cd–Cd in CdSe) against the same distance in the CIF. The bond length is not used: octahedral tilts shorten the Pb–Pb distance but not the Pb–X bond, so in an orthorhombic CsPbBr\ :sub:`3` dot the bond is 2.8 % longer than in the cubic CIF while the lattice is 0.6 % shorter. The log warns when the bonds are more than 8 % off the CIF, which points to the wrong CIF or a geometry of another material.

The k-path and the :math:`\mathbf{G}` replicas are expressed in the dot's frame, so the
folding uses the reciprocal lattice of the dot itself.


Displaying the Map
~~~~~~~~~~~~~~~~~~

The k-width of a state is physical: a dot of diameter :math:`D` has
:math:`\Delta k \approx 2\pi/D`, about a quarter of :math:`\Gamma`-X for a 2.5 nm CdSe dot, and
a typical state carries weight on most of the path. The energy width is not (``fuzzy_sigma``,
0.01-0.03 eV). Two things blur the map on top of that, and the dashboard removes both:

* The total weight :math:`\sum_\mathbf{k} I_m(\mathbf{k})` grows with the localisation of the
  state in real space; a cation d or surface state carries up to ten times the weight of a band
  edge state and dominates the colours. Each state is therefore normalised along the path,
  :math:`P_m(\mathbf{k}) = I_m(\mathbf{k}) / \langle I_m \rangle_\mathbf{k}` (1 = spread evenly),
  before the energy smearing.
* A logarithmic colour scale over four decades turns the low-weight tails into a haze. The map
  is drawn on a square-root scale of the normalised weight instead (``fuzzy_display_mode``
  ``state_norm``, the default of ``generate_interactive_plot``; ``raw`` is the earlier log scale
  of :math:`I_m`).

The buttons above the map switch between views:

``Map``
   the energy-smeared map of :math:`P_m(\mathbf{k})`.
``Map + states`` / ``States``
   one marker per maximum of each state's k-profile, at the state's own energy, with the
   marker area set by the weight (maxima above 1.5 times the mean of the state and 25 % of
   its largest peak, searched separately on each continuous piece of the path). Nothing is
   broadened, so this is the sharp version of the map: in the effective-mass picture an
   envelope state :math:`nl` peaks at :math:`|\mathbf{k}| \approx x_{nl}/R`, and the markers
   trace the bulk dispersion :math:`E(\mathbf{k})`.
``j-character`` (SOC)
   the map coloured by the spin-orbit energy of the spinors (below).

**j character of the spinors.** In the basis
:math:`(\phi_i\alpha, \phi_i\beta)` the spinor Hamiltonian is
:math:`\mathrm{diag}(\varepsilon) + V_{SOC}`, so the spin-orbit energy of spinor :math:`n`
follows from its eigenvalue and coefficients alone:

.. math::

   \langle V_{SOC} \rangle_n = E_n - \sum_i |U_{in}|^2 \varepsilon_i .

For anion p states it is :math:`+\Delta_{so}/3` for :math:`j = 3/2` (heavy and light holes)
and :math:`-2\Delta_{so}/3` for :math:`j = 1/2` (split-off); :math:`V_{SOC}` is traceless, so
it averages to zero over the SOC window. The SOC markers are coloured by it (red
:math:`j = 3/2`, blue :math:`j = 1/2`, white for the spin-free states outside the window), and
the ``j-character`` view colours the map by the weighted mean of
:math:`\langle V_{SOC}\rangle` in each pixel, with the brightness of the map. The value
per state is stored as ``soc/fuzzy/soc_energy_ev`` in ``qdex_electronic.h5``.

In the Cd\ :sub:`156`\ Se\ :sub:`111`\ Cl\ :sub:`90` dot (2.5 nm) the top 0.5 eV of the valence
band is :math:`j = 3/2` (:math:`\langle V_{SOC}\rangle` = +0.04 to +0.075 eV, against
:math:`\Delta_{so}/3` = +0.12 eV for bulk CdSe), but no state comes near the split-off value
:math:`-2\Delta_{so}/3 \approx -0.24` eV: the :math:`j = 1/2` character is spread over the
dense valence band from about 1.5 eV below the HOMO down (-0.005 eV per state on average). In an effective-mass
estimate confinement puts the lowest split-off envelope state about 1 eV below the valence
band maximum, among heavy- and light-hole states about 10 meV apart. In the bulk, k
conservation keeps the split-off band apart from the heavy- and light-hole states at the same
energy; in a dot the surface breaks it and the spin-orbit coupling mixes them. The mixing
weakens as the dot grows and k becomes a better quantum number.


Bulk Band Overlay
~~~~~~~~~~~~~~~~~

The dashboards draw the bulk band structure of the material (``qdex/data/bulk_bands``,
CP2K PBE with the same basis and pseudopotentials) over the fuzzy map. The bulk segments
are placed on the fuzzy path by their k-coordinates. A control bar above the plot shows or
hides the overlay and sets its colour, line width and opacity (for the unfolded perovskite
overlay the width scales the markers, whose opacity still follows the unfolding weight);
clicking the overlay's legend entry also toggles it. In the configuration, ``fuzzy.bulk_alignment``
(``core_level``, the default, ``midgap`` or ``vbm``) picks the energy alignment, ``fuzzy.bulk_unfolded:
off`` replaces the unfolded orthorhombic perovskite bands by the cubic ones and ``fuzzy.bulk_overlay:
false`` leaves the overlay out.

The anchor needs the semicore band in the MO file: for CsPbX\ :sub:`3` dots print all occupied
MOs (``MO_INDEX_RANGE 1 <HOMO + n_virtual>``); a window of the upper valence band (the Br 4p and
Pb 6s states, about 7 eV deep) stops above the semicore bands, and the bulk bands are then aligned
at mid-gap.

**Perovskites are anchored on the halide s band, not on Cs 5p.** The Cs 5p level follows the cage
around the A site: it lies 0.6–0.9 eV closer to the VBM in the orthorhombic bulk than in the cubic one,
and in a 5 × 5 × 5 CsPbX\ :sub:`3` cube it drifts by 0.6 eV from the core (−8.5 eV) to the surface
(−9.1 eV). The halide s level is flat across the dot to 0.1 eV. On the PBE-relaxed Cl/Br/I cubes
anchored on Cs 5p, the bulk VBM ended up 0.23 eV above the CsPbCl\ :sub:`3` HOMO (a confined hole
above the bulk band edge); on the halide s it is 0.12 eV above it, and for Br and I the holes are
confined by 0.12 and 0.05 eV, the electrons by about 0.27 eV (spin-free). With ``bulk_anchor: auto``
(the default) the manifolds of the A-site cation, which has no bond in the crystal bond star, are
tried last; ``bulk_anchor: Br-s`` (or any label of the bulk data) forces a manifold. Bulk-like atoms
are those with the crystal's number of *bonded* neighbours (two Pb for a halide), only crystal
elements within the bond cutoff, in the inner half of the dot.

The bulk bands are put on the dot's energy axis with a semicore level (from ``<name>.json``:
the cation d band, otherwise the anion s band; the perovskite A-site Cs 5p only as a last resort). For CdSe the Cd 4d
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

**Tilted perovskites: unfolded bands.** Lead halide perovskite nanocrystals relax with tilted
octahedra, as in the orthorhombic (Pnma) bulk phase. The fuzzy map keeps the cubic path (the
pseudo-cubic Brillouin zone, as in the literature): Pnma is a four-fold supercell of the cubic cell
in which the cubic R and M points both fold onto Γ, and a map on the orthorhombic path piles every
band edge there. The overlay is the orthorhombic band structure *unfolded* onto the cubic path
(``qdex.bulk_unfold``; Popescu and Zunger, *Phys. Rev. B* 85, 085201 (2012)): each supercell state
:math:`\psi_{K n}` gets the weight

.. math::

   P_{Kn}(\mathbf{k}) = \frac{\sum_{\mathbf{g} \in \text{cubic}} |\langle \mathbf{k}+\mathbf{g}|\psi_{Kn}\rangle|^2}
                             {\sum_{\mathbf{G} \in \text{supercell}} |\langle \mathbf{k}+\mathbf{G}|\psi_{Kn}\rangle|^2},

with the plane-wave amplitudes from the analytic Fourier transforms of the basis (the fuzzy-band
machinery) and the cubic lattice the average pseudo-cubic one of the supercell (nearest Pb–Pb
vectors). The dashboard draws one point per (k, band), with an opacity set by the weight: the
R-derived band edges are bright at R, their folded images faint elsewhere. Validation: untilted
cubic CsPbBr\ :sub:`3` written in the same supercell unfolds onto its cubic bands within 3.8 meV, with
weights 0 or 1 (shared only within degenerate sets). The orthorhombic structures are those of
:doc:`bulk_bands` (CsPbBr\ :sub:`3`, γ-CsPbI\ :sub:`3` measured; CsPbCl\ :sub:`3` with PBE-relaxed tilts),
scaled to the PBE pseudo-cubic volume like the dots; the files are
``qdex/data/bulk_bands/<cubic CIF stem>_unfolded.npz`` and are used whenever they exist
(``get_aligned_bulk_bands(..., unfolded="off")`` restores the cubic bands). They are anchored on
the Cs 5p level, like the cubic bands.

:math:`H(\mathbf{k}) = \sum_\mathbf{R} e^{i\mathbf{k}\cdot\mathbf{R}} \langle\chi_\mu(0)|H|\chi_\nu(\mathbf{R})\rangle`
and :math:`S(\mathbf{k})` are built at every point of the path, the SOC operator from
Bloch sums of the projector overlaps, and the two-component problem is solved in the full
AO basis. The run checks the AO basis against CP2K's overlap and the spin-free bands
against CP2K's band structure. For zinc-blende CdSe (PBE, DZVP-MOLOPT-PBE-GTH): the
spin-free bands agree with CP2K to 0.06 meV and :math:`\Delta_{so}(\Gamma) = 0.36` eV.
