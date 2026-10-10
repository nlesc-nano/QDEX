Configuration Reference (YAML)
==============================

A QDEX input file is organised by what each block controls:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Section
     - Controls
   * - ``system``
     - input files, material, threads, device
   * - ``environment``
     - the dielectric environment (``eps_out``), seen by the QP correction and the BSE kernel
   * - ``quasiparticles``
     - the QP correction of the orbital energies (:doc:`/quasiparticles/index`)
   * - ``integrals``
     - the representation of the two-electron integrals, MNOK or ZDO xs
       (:doc:`/integrals/index`)
   * - ``excitations``
     - the excited-state framework, the BSE kernel and the active space (:doc:`/excitons/index`)
   * - ``soc``
     - spin–orbit coupling
   * - ``analysis``, ``output``
     - populations, PDOS/COOP, fuzzy bands; spectra, cubes, CSV and NTO files
   * - ``periodic``, ``auger``, ``namd``
     - periodic images, Auger rates, non-adiabatic dynamics

Every key can also be given on the command line; the flag is shown in brackets below, and a flag
overrides the file. Files in the old layout (one ``physics:`` block with the command-line key names,
plus ``solver:`` and ``fuzzy:``) are still read, with a notice.

Example
-------

.. code-block:: yaml

   system:
     mo_file: "MOs_cleaned_20ang.txt"
     xyz: "geom.xyz"
     basis_txt: "BASIS_MOLOPT_UZH"
     basis_name: "DZVP-MOLOPT-PBE-GTH"
     material: "CDSE"
     gth_file: "GTH_SOC_POTENTIALS.txt"
     nthreads: 12
     skip_orthonormality_check: true

   environment:
     eps_out: 2.24                  # toluene; 1.0 = vacuum

   quasiparticles:
     model: "sgw-resta"             # QP correction (see below)

   integrals:
     representation: "mnok"         # mnok or xs
     charges: "mulliken"            # mulliken or lowdin

   excitations:
     mode: "bse"
     nhomos: 25
     nlumos: 25
     nroots: 40

   soc:
     enabled: true

   output:
     plot: true
     write_csv: true

system
------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key [flag]
     - Meaning
   * - ``mo_file``, ``mo_file_beta``
     - CP2K MO file (text, ``.gz``, ``.mbse`` or TREXIO HDF5 ``.h5``); beta MOs for unrestricted runs.
       For an unrestricted ``.h5`` give the same file for both; the α and β channels are read
       separately
   * - ``xyz``
     - geometry; optional with an ``.h5`` MO file, whose ``nucleus`` group is read instead
   * - ``basis_txt``, ``basis_name``
     - basis-set file and basis name; ``basis_name: per-atom`` reads a per-atom basis file written by
       ``qdex.xtb.molden`` (g-xTB orbitals, :doc:`/dynamics/gxtb`)
   * - ``material`` [``--material``]
     - ``MATERIAL_DB`` entry (bulk gaps, ε∞, effective masses), e.g. ``CDSE``, ``CSPBBR3``. CdSe and CdS
       have one entry per polytype: ``CDSE``/``CDS`` (= ``CDSE_ZB``/``CDS_ZB``) zinc blende,
       ``CDSE_WZ``/``CDS_WZ`` wurtzite (different bulk limits; :doc:`/electronic_structure/bulk_bands`)
   * - ``gth_file``
     - GTH SOC pseudopotential file
   * - ``nthreads``, ``device``
     - CPU threads; ``auto``, ``cpu``, ``cuda`` or ``mps``
   * - ``lowdin_device`` [``--lowdin-device``]
     - where the overlap matrix S is diagonalized for the Löwdin quantities: ``auto`` (the compute device),
       ``cpu`` or ``cuda``/``cuda:N``. On the CPU, LAPACK divide and conquer; past 32 767 AOs the 64-bit-integer
       ``dsyevd_64`` of MKL (oneMKL ≥ 2023), else scipy's slower syevr. On a GPU, cuSOLVER ``cusolverDnXsyevd``
       (cuSOLVER ≥ 12 beyond ~32k AOs; library from ``$QDEX_CUSOLVER`` or next to torch), when ~40 n² bytes
       are free on the GPU, otherwise the CPU. 42.6k AOs: syevr > 2000 s, MKL ILP64 900 s (192 cores), one
       H200 28 s. ``$QDEX_LOWDIN_CACHE`` (a directory) keeps the result for later runs on the same orbitals
   * - ``skip_orthonormality_check``
     - skip the Cᵀ S C test of the MO file (saves one n_ao³ product)
   * - ``inorganic_elements``
     - elements seen by SAXS for the reported size (default: all but H, C, N, O, P, B, Si, F;
       :doc:`/reference/cluster_size`)
   * - ``cache_mos``, ``log_file``, ``orthonormality_tol``
     - binary MO cache, log file name, tolerance of the orthonormality test

environment
-----------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key [flag]
     - Meaning
   * - ``eps_out`` [``--eps-out``]
     - optical dielectric constant outside the dot (default 2.0; vacuum 1, toluene 2.24). It enters the
       sphere reaction field of ΔW, hence both the QP correction and the BSE kernel.

quasiparticles
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 16 54

   * - Key [flag]
     - Default
     - Meaning
   * - ``model`` [``--qp_gap``]
     - ``brus``
     - ``sgw-resta``, ``sgw-dim``, ``evgw-resta``, ``evgw-dim``, ``qsgw-resta``, ``qsgw-dim``; ``none``
       (DFT energies as they are; also ``pbe``, ``dft``), ``bulk`` (PBE energies + bulk GW correction,
       PBE orbitals only; with a bulk kernel this is the sBSE), ``brus``, or a gap
       in eV (:doc:`/quasiparticles/models`)
   * - ``reference`` [``--qp-reference``]
     - ``pbe``
     - orbitals the ``bulk`` shift corrects: ``pbe``, or ``gxtb`` for g-xTB frames in the NAMD
       precompute (:doc:`/dynamics/gxtb`)
   * - ``bulk_vertex`` [``--bulk-vertex``]
     - ``none``
     - vertex correction of the bulk QSGW opening: ``none`` (pure QSGW), ``full`` (bulk factor at every
       size) or ``scaled`` (times the Penn fraction of bulk screening; :doc:`/quasiparticles/gw`)
   * - ``bulk_vertex_factor`` [``--bulk-vertex-factor``]
     - 0.8
     - bulk factor a: Δ_Σ → a Δ_Σ in the bulk; ``material``: the factor that puts the bulk limit on the
       experimental gap of the material (:doc:`/quasiparticles/gw`)
   * - ``bulk_residual`` [``--bulk-residual``]
     - ``none``
     - ``experimental``: a constant shift that puts the bulk limit on the room-temperature gap; with
       ``bulk_vertex: scaled`` and a numeric factor, the split model (:doc:`/quasiparticles/gw`)
   * - ``bulk_edge_split`` [``--bulk-edge-split``]
     - ``symmetric``
     - HOMO share f_b of the bulk QP correction; places IP and EA, leaves the gap and the BSE unchanged:
       ``symmetric`` (1/2), ``cluster`` (calibrated on evGW of small clusters), ``bulk`` (bulk GW and
       experiment) or a number (:doc:`/workflows/qp_edges`)
   * - ``ip_ea`` [``--qp-ip-ea``]
     - ``resta``
     - IP/EA of ``bulk``: ``resta`` adds the ``sgw-resta`` self-energy of the HOMO and LUMO (report only),
       ``rigid`` the bulk shift only (:doc:`/quasiparticles/models`)
   * - ``bulk_geometry`` [``--bulk-geometry``]
     - ``strain``
     - geometry correction of the bulk shift for PBE-relaxed dots: ``strain`` (bulk PBE gap change between
       a_exp and the lattice measured on the dot's interior bonds), ``full`` (dot at the bulk PBE lattice)
       or ``none`` (:doc:`/electronic_structure/bulk_bands`)
   * - ``z`` [``--qp-z``]
     - ``derived``
     - quasiparticle weight: plasmon pole of the model's ε, or a number
   * - ``selfenergy`` [``--qp-selfenergy``]
     - ``cohsex``
     - one-shot ΔCOHSEX or ``classical`` ½ qᵀΔWq
   * - ``levels`` [``--qp-levels``]
     - ``orbital``
     - correct every orbital, or one ``rigid`` scissor
   * - ``window`` [``--qp-window``], ``window_size`` [``--qp-window-size``]
     - ``active``
     - orbitals evaluated explicitly by ΔCOHSEX: the BSE active space (``active``) or ``all``
   * - ``cohsex_occ_window`` [``--cohsex-occ-window``]
     - all occupied
     - screened-exchange density of ΔCOHSEX from the occupied MOs within this many eV of the HOMO only; 5 eV
       changes the QP gap and E₁ by ≤ 5 meV (CdSe, 1.5–3.5 nm) and lets the h5 hold only those occupied MOs
   * - ``solvent_term`` [``--qp-solvent-term``]
     - ``sphere``
     - environment part of ΔW: dielectric-sphere reaction field, or the older ``born`` form
   * - ``radius`` [``--qp-radius``]
     - ``saxs``
     - radius of the dielectric sphere and of the Brus confinement: SAXS-equivalent (``saxs``) or core
       hull + 1.25 Å (``hull``); the two-anchor ``gw`` always uses the hull
   * - ``energy_reference``
     - ``vacuum``
     - QP energies relative to the vacuum or to the Fermi level (``fermi``)

integrals
---------

.. list-table::
   :header-rows: 1
   :widths: 30 16 54

   * - Key [flag]
     - Default
     - Meaning
   * - ``representation`` [``--2e-integrals``]
     - ``mnok``
     - ``mnok`` (atom-condensed, MNOK γ) or ``xs`` (ZDO, exact (μμ|νν)). Used for ΔW in the QP correction
       and for K\ :sup:`x`, K\ :sup:`d` in the BSE.
   * - ``charges`` [``--charge_type``]
     - ``mulliken``
     - population partition of the MNOK densities: ``mulliken`` or ``lowdin`` (xs always uses Löwdin)
   * - ``mnok_exponent`` [``--mnok-exponent``]
     - 2
     - exponent β of the MNOK interaction (r^β + a^β)^(−1/β), all interactions
   * - ``mnok_exponent_exchange`` [``--mnok-exponent-exchange``]
     - = mnok_exponent
     - β of the exchange interaction only
   * - ``mnok_onsite`` [``--mnok-onsite``]
     - ``ip_ea``
     - on-site value IP − EA = 2η (``ip_ea``) or η (``eta``, earlier convention);
       :doc:`/integrals/representation`
   * - ``beta`` [``--beta``]
     - 0
     - MNOK parameter of the bare γ

excitations
-----------

.. list-table::
   :header-rows: 1
   :widths: 30 16 54

   * - Key [flag]
     - Default
     - Meaning
   * - ``mode`` [``--excitation-mode``]
     - ``bse``
     - ``independent_dft``, ``independent_qp``, ``diagonal_bse``, ``bse``, ``stda``, ``diagonal_stda``; ``sbse`` and ``diagonal_sbse``
       are the same solvers, named for use with ``quasiparticles.model: bulk`` and a bulk kernel
       (:doc:`/excitons/index`)
   * - ``functional`` [``--stda-functional``], ``ax`` [``--stda-ax``]
     - —
     - sTDA only: functional of the MO file (sets a_x; a range-separated one such as ``wb97m-v`` or
       ``gxtb`` sets a_x, α and β), or a_x as a number or ``dielectric`` (1/ε∞) (:doc:`/excitons/stda`)
   * - ``stda_alpha`` [``--stda-alpha``], ``stda_beta`` [``--stda-beta``]
     - —
     - sTDA only: explicit exponents of γ\ :sup:`K` and γ\ :sup:`J`; override the preset or the
       global-hybrid formulas
   * - ``kernel`` [``--kernel``]
     - model default
     - W of K\ :sup:`d`. Set by the QP model (``qp`` for Resta and DIM). Models without W (``none``, ``bulk``,
       ``brus``, ``pbe``, a gap) take ``resta``, ``dim``, ``rpa``, ``sbse``, ``xs-*`` or ``bse``
       (:doc:`/excitons/kernel`)
   * - ``nhomos``, ``nlumos``
     - all
     - active occupied and virtual orbitals
   * - ``e_thresh``, ``f_thresh``
     - —
     - transition-energy cutoff (eV); minimum oscillator strength printed
   * - ``nroots``, ``full_diag``, ``tol``
     - 10, false, 10⁻⁵
     - roots, dense instead of Davidson, Davidson tolerance
   * - ``selection`` [``--selection``], ``selection_energy``, ``selection_pt``
     - ``none``, 7.0, 10⁻⁴
     - ``perturbative``: Grimme's selection of transitions from the active space for the coupled
       solvers; E_thr in eV, t in hartree (:doc:`/excitons/bse`)
   * - ``triplet``
     - false
     - triplet BSE (no K\ :sup:`x`)
   * - ``include_exchange``, ``include_direct_eh``
     - true
     - switch K\ :sup:`x` or K\ :sup:`d` off
   * - ``kernel_scaling`` [``--alpha``]
     - 1.0
     - scale factor of the uniform ``bse`` kernel
   * - ``allow_inconsistent_kernel``
     - false
     - allow a kernel different from the QP model's W (legacy results only)
   * - ``energy_shift`` [``--soc``]
     - 0
     - empirical shift subtracted from the excitation energies

soc
---

* ``enabled`` [``--soc_flag``]: two-component spinor BSE.
* ``window`` [``--soc_window``]: energy window (eV) around the Fermi level for full SOC mixing in the fuzzy bands (default: largest ``|ewin|`` + 1 eV).
* ``bse_window`` [``--soc_bse_window``]: SOC window the BSE spinors are projected from (default: ``window`` when set, else largest ``|ewin|`` + 2 eV).
* The GTH SOC file is ``system.gth_file``.

analysis and output
-------------------

* ``analysis``: ``run_fuzzy``, ``cif``, ``pdos_atoms``, ``coop_pairs``, ``pdos_sigma``, ``fuzzy_sigma``,
  ``ewin``, ``fold_to_bz``, ``g_shell``, ``population_print_range``, ``dashboard_energy_mode``,
  ``bulk_overlay``, ``bulk_alignment``, ``bulk_anchor``, ``bulk_unfolded``
  (:doc:`/electronic_structure/index`).
* ``output``: ``plot``, ``show``, ``broadening``, ``sigma``, ``write_csv``, ``csv_roots``, ``save_xia``,
  ``time``, ``cube``, ``cube_spacing``, ``cube_nhomos``, ``cube_nlumos``, ``cube_nhomos_deloc``,
  ``cube_nlumos_deloc``, ``nbse``, ``bse_states``,
  ``nto``, ``nto_states``, ``nto_top``, ``nto_csv`` (:doc:`/exciton_analysis/index`), ``verbosity``.

Database output
^^^^^^^^^^^^^^^

For runs that feed a database (many dots, results read by programs rather than
people), ``output`` has:

* ``h5: true`` writes two HDF5 files (energies in eV, lengths in Å, transition
  dipoles in e·bohr; the QDEX commit and all settings as file attributes):

  - ``qdex_electronic.h5``: ``structure``; ``qp`` (DFT and QP gaps, HOMO/LUMO
    vs vacuum, spin-free and SOC); ``sf/mo`` (every MO of the orbital file:
    DFT and QP energies, occupation, element / element(l) / surface fractions,
    IPR, COOP per pair, PDOS on a grid); ``soc/bse_spinor`` (the exciton active
    space spinors, same projections, absolute and QP energies) and
    ``soc/spinor`` (the fuzzy-window spinors); ``sf/fuzzy`` and ``soc/fuzzy``
    (raw weights :math:`|\langle\phi_n|k\rangle|^2` along the k-path, to be
    smeared with the stored ``sigma_ev``).
  - ``qdex_excitations.h5``: ``sf`` and ``soc``, every exciton up to
    ``excitations_emax`` (default 4.5 eV, about 275 nm, measured on the exciton
    energy) and at least ``excitations_min_states`` (default 20): energy,
    oscillator strength, transition dipole, hole and electron orbital (MO or
    spinor index), D / Kx / −Kd, d_eh, d_CT, σ_h, σ_e, CT character, type,
    singlet fraction (SOC); and the absorption spectrum of every computed state on
    a fixed grid (``spectrum_grid``, default 0.5-6.0 eV in 5 meV steps) for each of
    ``spectrum_sigmas`` (default 0.03 and 0.10 eV). In the diagonal modes
    every transition is kept whatever ``nroots`` is, at no extra cost; the
    eigenvectors are not stored (one hole-electron pair per state).

* ``html: false`` skips the HTML dashboards and the CSV / NPZ files they read.
* ``mo_cubes: true`` writes the spin-free ``cube_nhomos`` + ``cube_nlumos`` MOs
  around the gap (HOMO-1, HOMO, LUMO, LUMO+1 by default) on a coarse
  ``mo_cube_spacing`` grid (default 0.8 Å), without the spinor and exciton
  cubes of ``cube``.
* ``cube_nhomos_deloc`` / ``cube_nlumos_deloc`` (default 0) also write the first
  N band states from the delocalized HOMO down and from the delocalized LUMO up
  (the trap detector of the fuzzy bands, so ``analysis.run_fuzzy`` and
  ``analysis.cif`` are needed), as ``spatial_MO_HOMO-20_dHOMO.cube``,
  ``spatial_MO_HOMO-21_dHOMO-1.cube``, ... An MO that is already among the
  nominal cubes (``cube`` / ``mo_cubes`` with ``cube_nhomos`` and ``cube_nlumos``)
  is not written again: without traps the delocalized and nominal edges coincide
  and no extra cube is made. The spin-free dashboard shows them in a second row
  of its 3D panel.
* ``verbosity: quiet`` keeps the console to warnings; the log file keeps everything.

The HTML pages of a database run can be made later from its two HDF5 files,
without recomputing anything::

   qdex dashboards <run dir> [-o <out dir>] [--no-cubes] [--keep-data]

This writes ``fuzzy_dashboard_sf/soc.html`` (with the MO cubes of the run
directory in a 3D panel, several MB; ``--no-cubes`` leaves it out),
``exciton_analysis_sf/soc.html`` (the stored excitons, at most 100) and the
spectrum plots. The data drawn are those of a direct run with ``html: true`` to
float32 rounding.

.. code-block:: yaml

   output:
     h5: true
     html: false
     plot: false
     write_csv: false
     mo_cubes: true
     verbosity: quiet

**Output verbosity.** ``output.verbosity`` (``--verbosity``) sets what reaches the console:
``full`` (default, everything), ``normal`` (without iteration traces, timings and diagnostics) or
``quiet`` (warnings and errors only). The log file (``system.log_file``, default ``minibse.log``)
always receives the full output, so a quiet run, for example a NAMD precompute, still leaves a complete
record. The messages go through Python's ``logging`` module (logger ``qdex``); scripts that import QDEX
can change the level with ``qdex.logging_setup.set_verbosity``.

auger
~~~~~
* **run** (*bool*): Enable Auger recombination calculations.
* **sigma** (*float*): Gaussian energy conservation broadening width in eV (default: ``0.05``).
* **channel** (*str*): Recombination channel: ``"all"``, ``"eeh"``, or ``"hhe"``.
* **n_initial_states** (*int*): Number of frontier band-edge states to consider as initial carriers (default: ``1``).
* **lineshape** (*str*): Energy conservation model: ``"gaussian"`` or ``"fcwd"`` (Marcus multi-phonon line shape).

namd
~~~~
* **trajectory_dir** (*str*): Path to folder containing trajectory MD frames.
* **dt_fs** (*float*): Nuclear time step in femtoseconds between frames.
* **temperature_k** (*float*): Lattice temperature in Kelvin for detailed balance.
* **method** / **engine** (*str*): ``"master_equation"`` (deterministic tensorized PME), ``"cpa_fssh"`` (CPA-FSSH with continuous EDC), or ``"dish"`` (Decoherence-Induced Surface Hopping).
* **tau_dec_fs** (*str or float*): ``"edc"`` / ``"cumulant"`` (*ab initio* from energy gap fluctuations) or fixed float in fs.
* **decoherence** (*str*): Decoherence scheme for FSSH (``"edc"`` for Granucci-Persico continuous energy-based damping).
* **recombination.include_ground_state** (*bool*): Couple excited states to the ground state.
* **recombination.radiative** (*bool*): Use *ab initio* Einstein spontaneous emission formula for :math:`k_{\mathrm{rad}}`.
* **recombination.tau_nr_ns** (*float, optional*): Defect trap non-radiative lifetime in nanoseconds for PLQY computation.
