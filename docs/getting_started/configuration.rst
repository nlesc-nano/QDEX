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
       (:doc:`/quasiparticles/representation`)
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
     - CP2K MO file (text, ``.gz`` or ``.mbse``); beta MOs for unrestricted runs
   * - ``xyz``
     - geometry
   * - ``basis_txt``, ``basis_name``
     - basis-set file and basis name
   * - ``material`` [``--material``]
     - ``MATERIAL_DB`` entry (bulk gaps, ε∞, anchor, masses), e.g. ``CDSE``, ``CSPBBR3``
   * - ``gth_file``
     - GTH SOC pseudopotential file
   * - ``nthreads``, ``device``
     - CPU threads; ``auto``, ``cpu``, ``cuda`` or ``mps``
   * - ``skip_orthonormality_check``
     - skip the Cᵀ S C test of the MO file (saves one n_ao³ product)
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
       (KS + bulk GW correction, no ΔW; with a bulk kernel this is the sBSE), ``brus``, ``pbe``, or a gap
       in eV (:doc:`/quasiparticles/models`)
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
   * - ``solvent_term`` [``--qp-solvent-term``]
     - ``sphere``
     - environment part of ΔW: dielectric-sphere reaction field, or the older ``born`` form
   * - ``anchor_residual`` [``--qp-anchor-residual``]
     - ``on``
     - add the calibrated evGW\@PBE0 residual (:doc:`/quasiparticles/anchor`)
   * - ``residual_scaling`` [``--qp-residual-scaling``], ``residual_power``
     - ``econf``
     - size scaling of the residual: E_conf(R)/E_conf(R₀), or ``power`` (R₀/R)^p
   * - ``anchor_calibrate`` [``--qp-anchor-calibrate``], ``anchor_table``
     - off
     - run on the anchor cluster to store the model's residual
   * - ``edge_split`` [``--qp-edge-split``]
     - ``anchor``
     - HOMO/LUMO split for absolute IP/EA
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
     - ``independent_dft``, ``independent_qp``, ``diagonal_bse``, ``bse``; ``sbse`` and ``diagonal_sbse``
       are the same solvers, named for use with ``quasiparticles.model: none`` and a bulk kernel
       (:doc:`/excitons/frameworks`)
   * - ``kernel`` [``--kernel``]
     - model default
     - W of K\ :sup:`d`. Set by the QP model (``qp`` for Resta and DIM). Models without W (``none``,
       ``brus``, ``pbe``, a gap) take ``resta``, ``dim``, ``rpa``, ``sbse``, ``xs-*`` or ``bse``
       (:doc:`/excitons/screened_kernel`)
   * - ``nhomos``, ``nlumos``
     - all
     - active occupied and virtual orbitals
   * - ``e_thresh``, ``f_thresh``
     - —
     - transition-energy cutoff (eV); minimum oscillator strength printed
   * - ``nroots``, ``full_diag``, ``tol``
     - 10, false, 10⁻⁵
     - roots, dense instead of Davidson, Davidson tolerance
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
* ``window`` [``--soc_window``]: energy window (eV) around the Fermi level for full SOC mixing.
* The GTH SOC file is ``system.gth_file``.

analysis and output
-------------------

* ``analysis``: ``run_fuzzy``, ``cif``, ``pdos_atoms``, ``coop_pairs``, ``pdos_sigma``, ``fuzzy_sigma``,
  ``ewin``, ``fold_to_bz``, ``g_shell``, ``population_print_range``, ``dashboard_energy_mode``
  (:doc:`/electronic_structure/index`).
* ``output``: ``plot``, ``show``, ``broadening``, ``sigma``, ``write_csv``, ``csv_roots``, ``save_xia``,
  ``time``, ``cube``, ``cube_spacing``, ``cube_nhomos``, ``cube_nlumos``, ``nbse``, ``bse_states``,
  ``nto``, ``nto_states``, ``nto_top``, ``nto_csv`` (:doc:`/exciton_analysis/index`).

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
