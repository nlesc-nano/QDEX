Configuration Reference (YAML)
==============================

`QDEX` uses structured YAML files to control calculations. Below is the complete reference of all supported configuration sections and keywords.

Complete Example Configuration
------------------------------

.. code-block:: yaml

   system:
     mo_file: "CsPbBr3_MOs.mbse"
     xyz: "CsPbBr3_QD.xyz"
     basis_txt: "BASIS_MOLOPT_UZH"
     basis_name: "DZVP-MOLOPT-PBE-GTH"
     material: "CSPBBR3"
     nthreads: 12
     device: "auto"                    # "auto", "cpu", "cuda", or "mps"
     gth_file: "GTH_SOC_POTENTIALS.txt" # Required if soc: true
     cif: "CsPbBr3_bulk.cif"           # Required for fuzzy band unfolding

   physics:
     excitation_mode: "diagonal_bse"   # "bse", "diagonal_bse", "independent_qp", "independent_dft"
     qp_gap: "sgw-anchor"              # "sgw-anchor" (or "gw"), "sgw-dim", "evgw-dim", "qsgw-dim", "brus", "pbe"
     2e-integrals: "mnok"              # "mnok" atom charges or "xs" analytical AO density-pair integrals
     kernel: "resta"                   # "resta", "dim", "rpa", "sbse", or "bse"
     dynamic_z: false                  # Empirical state-dependent Z_p damping
     update_orbitals: false            # Static AO-basis orbital relaxation model
     include_direct_eh: true           # Attractive screened electron-hole term Kd
     include_exchange: true            # Repulsive bare electron-hole term Kx
     eps_out: 2.25                     # Solvent / external dielectric constant
     nhomos: 50                        # Active occupied MOs (or all if omitted)
     nlumos: 50                        # Active virtual MOs (or all if omitted)
     soc: true                         # Enable 2-component spinor Spin-Orbit Coupling
     soc_window_ev: 8.0                # Energy window around Fermi level for SOC active space
     triplet: false                    # Perform triplet BSE (omitting 2*K_x)
     charge_type: "mulliken"           # "mulliken" or "lowdin"

   bse:
     nroots: 10                        # Number of exciton roots to compute
     full_diag: false                  # Use Davidson iterative solver (false) or dense LAPACK (true)
     tol: 1.0e-5                       # Convergence tolerance for Davidson residual norm

   fuzzy:
     run: true                         # Run supercell unfolding and PDOS/COOP
     pdos_atoms: ["Pb", "Br", "Cs"]    # Element symbols for projected DOS
     coop_pairs: ["Pb-Br", "Cs-Br"]    # Atom pairs for Crystal Orbital Overlap Population
     ewin: [-4.0, 4.0]                 # Energy window relative to Fermi level (in eV)
     pdos_sigma: 0.08                  # Gaussian broadening for PDOS (in eV)
     fuzzy_sigma: 0.03                 # Gaussian broadening along energy axis for fuzzy bands (in eV)
     fold_to_bz: true                  # Fold into 1st Brillouin zone by summing reciprocal replicas
     g_shell: 1                        # Shell of reciprocal replicas (0: 1, 1: 27, 2: 125 replicas)
     dashboard_energy_mode: "both"     # "dft", "qp", or "both"
     qp_energy_reference: "vacuum"     # "vacuum" (absolute IP/EA) or "fermi"

   cube:
     export: true                      # Generate volumetric Gaussian .cube files
     spacing_ang: 0.4                  # Grid spacing in Angstroms
     nhomos: 2                         # Number of HOMO orbitals to export
     nlumos: 2                         # Number of LUMO orbitals to export
     bse_states: [1, 2]                # Specific BSE exciton states to export (hole, electron, diff)

   analysis:
     nto: true                         # Natural Transition Orbital analysis
     nto_states: [1, 2, 3]             # Specific states for NTO analysis (1-indexed)
     nto_top: 3                        # Dominant NTO pairs to print and tabulate
     nto_csv: true                     # Write detailed NTO descriptors to CSV
     plot: true                        # Generate publication plots and Plotly HTML dashboards

   auger:
     run: true                         # Compute non-radiative Auger recombination rates
     sigma: 0.05                       # Energy conservation broadening in eV
     channel: "all"                    # "all", "eeh", or "hhe"
     n_initial_states: 1               # Number of frontier band-edge carriers
     lineshape: "gaussian"             # "gaussian" or "fcwd"

   namd:
     trajectory_dir: "./trajectory"    # Directory containing frame_0001, frame_0002...
     dt_fs: 1.0                        # Nuclear MD time step in femtoseconds
     temperature_k: 300.0              # Lattice temperature for detailed balance
     engine: "master_equation"         # "master_equation" (PME) or "surface_hopping" (FSSH)
     tau_dec_fs: "cumulant"            # "cumulant" (ab initio) or fixed float in fs
     
     recombination:
       include_ground_state: true
       radiative: true                 # Uses Einstein spontaneous emission formula
       tau_nr_ns: 25.0                 # Non-radiative defect trap lifetime in nanoseconds
     
     storage:
       precompute_dir: "namd_precomputed"
       output_dir: "namd_results"

Detailed Keyword Reference
--------------------------

system
~~~~~~
* **mo_file** (*str*): Path to CP2K binary molecular orbitals (``.mbse``) or formatted text file.
* **xyz** (*str*): Path to Cartesian coordinates (``.xyz``) of the system.
* **basis_txt** (*str*): Path to CP2K Gaussian basis set file (e.g. ``BASIS_MOLOPT``).
* **basis_name** (*str*): Name of basis set to extract (e.g. ``DZVP-MOLOPT-PBE-GTH``).
* **material** (*str*): Target material key in ``MATERIAL_DB`` (e.g. ``CSPBBR3``, ``CDSE``, ``INAS``).
* **cif** (*str, optional*): Path to reference bulk crystallographic unit cell (for high-symmetry :math:`k`-path generation in Fuzzy Bands).
* **nthreads** (*int*): Number of OpenMP / BLAS CPU threads.
* **device** (*str*): Computation device: ``"auto"``, ``"cpu"``, ``"cuda"``, or ``"mps"``.
* **gth_file** (*str, optional*): Path to custom GTH SOC pseudopotential parameter file.

physics
~~~~~~~
* **excitation_mode** (*str*):
  - ``"bse"``: Full Tamm-Dancoff Bethe-Salpeter Equation with configuration interaction.
  - ``"diagonal_bse"``: Diagonal bare exchange :math:`K_x` and screened direct attraction :math:`K_d` without off-diagonal coupling. Optimal for dense NAMD.
  - ``"independent_qp"``: Non-interacting single-particle transitions with scaled GW scissor gap.
  - ``"independent_dft"``: Non-interacting single-particle transitions with bare DFT gap.
* **qp_gap** (*str or float*):
  - ``"sgw-anchor"`` (or ``"gw"``): Scaled GW model with two anchors (vacuum cluster and bulk limit) and dielectric polarization.
  - ``"sgw-dim"``: Microscopic atomistic polarizable dipole screening difference model.
  - ``"sgw-resta"``: Resta electronic Thomas-Fermi screening model with Penn size-dependent dielectric scaling.
  - ``"evgw-dim"`` / ``"evgw-resta"``: Effective-gap self-consistent screening models.
  - ``"qsgw-dim"`` / ``"qsgw-resta"``: Static AO-basis COHSEX-like orbital-relaxation models; these are not conventional QSGW.
  - ``"brus"``: Brus effective mass confinement model.
  - ``"pbe"``: Uncorrected DFT eigenvalues.
  - *float*: Explicit user-defined target band gap in eV.
* **2e-integrals** (*str*): Two-electron integral representation: ``"mnok"`` (semi-empirical atom-centered transition charges) or ``"xs"`` (analytical Libint2 integrals restricted to AO density pairs :math:`(\mu\mu|\nu\nu)`). Also accepted as ``two_electron_integrals`` or legacy ``kernel_type``.
* **kernel** (*str*): Screened interaction :math:`W` of the BSE direct term.

  - ``"qp"``: the :math:`W` built by the QP model. It is the default, and the only allowed choice, for the Delta-W QP models (``sgw-*``, ``evgw-*``, ``qsgw-*``, ``sgw``), so that GW and BSE share one :math:`W`. Any other kernel with these models is rejected; see :doc:`/validation/model_comparison`.
  - ``"resta"``, ``"dim"``, ``"xs-resta"``, ``"sbse"``, ``"bse"`` (unscreened MNOK, legacy default): independent kernels for the gap-only QP models (``pbe``, ``brus``, ``gw``, a numeric gap).
* **qp_z** (*str or float*): Quasiparticle weight for the Delta-W models. ``"derived"`` (default for all Delta-W models) uses one plasmon pole built from the model's own dielectric constant and the material's valence plasmon; a number (e.g. ``1.0`` or ``0.8``) is used as a fixed value. The same Z scales Delta-W in the BSE kernel. Rejected for gap-only models. CLI: ``--qp-z``.

* **qp_polarization** (*str*): finite-size term of the ``gw`` model: ``"sphere"`` (default; classical surface polarization of a dielectric sphere) or ``"legacy"`` (11.52 eV Å (1/ε_out − 1/ε_∞)/(R + ℓ)). CLI: ``--qp-polarization``.
* **qp_edge_split** (*str*): HOMO/LUMO split for absolute IP/EA: ``"anchor"`` (default; per-edge two-anchor curves) or ``"model"`` (the Delta-W model's own split). CLI: ``--qp-edge-split``.
* **allow_inconsistent_kernel** (*bool*): Allow a BSE kernel different from the QP model's :math:`W`. Only for reproducing results obtained before this check existed. CLI: ``--allow-inconsistent-kernel``.
* **dynamic_z** (*bool*): Same as ``qp_z: derived`` (empirical state-dependent :math:`Z_p`; no frequency-dependent self-energy or :math:`f`-sum rule is evaluated).
* **update_orbitals** (*bool*): Relax molecular orbitals in the AO basis using the static COHSEX-like ``qsgw-*`` model.
* **include_direct_eh** (*bool*): Include screened attractive :math:`K^d` in ``bse`` or ``diagonal_bse`` (default ``true``). CLI: ``--include-direct-eh`` / ``--no-direct-eh``.
* **include_exchange** (*bool*): Include bare repulsive :math:`K^x` independently of :math:`K^d` (default ``true``). CLI: ``--include-exchange`` / ``--no-exchange``.
* **soc** (*bool*): Enable fully relativistic 2-component spinor Hamiltonian.
* **soc_window_ev** (*float*): Energy window in eV around the Fermi level for selecting active MOs in SOC.
* **eps_out** (*float*): Surrounding solvent or matrix dielectric constant (default 2.0).
* **nhomos** / **nlumos** (*int*): Number of occupied and virtual frontier molecular orbitals to include in the active space.
* **triplet** (*bool*): Perform triplet BSE calculation (omits repulsive exchange :math:`2K^x`).
* **charge_type** (*str*): Method for computing transition charges: ``"mulliken"`` or ``"lowdin"``.

bse
~~~
* **nroots** (*int*): Number of lowest excited states to compute.
* **full_diag** (*bool*): If ``true``, uses dense LAPACK diagonalization; if ``false``, uses the memory-efficient Davidson iterative solver.
* **tol** (*float*): Convergence tolerance for the Davidson solver residual norm (default :math:`10^{-5}`).

fuzzy
~~~~~
* **run** (*bool*): Enable supercell unfolding (Fuzzy Bands) and PDOS/COOP analysis.
* **pdos_atoms** (*list of str*): Atomic species for PDOS projection (e.g. ``["Pb", "Br"]``).
* **coop_pairs** (*list of str*): Atom pairs for COOP bonding analysis (e.g. ``["Pb-Br"]``).
* **ewin** (*list of float*): Energy window :math:`[E_{\min}, E_{\max}]` in eV relative to Fermi level.
* **pdos_sigma** (*float*): Gaussian broadening standard deviation for continuous PDOS in eV.
* **fuzzy_sigma** (*float*): Gaussian broadening along energy axis for fuzzy bands in eV.
* **fold_to_bz** (*bool*): Fold spectral weights into the first Brillouin zone using reciprocal replicas.
* **g_shell** (*int*): Reciprocal lattice vector shell for folding (0: 1 replica, 1: 27 replicas, 2: 125 replicas).
* **dashboard_energy_mode** (*str*): ``"dft"``, ``"qp"``, or ``"both"``.
* **qp_energy_reference** (*str*): ``"vacuum"`` (absolute IP/EA) or ``"fermi"`` (:math:`E_F = 0`).

cube
~~~~
* **export** (*bool*): Generate 3D volumetric Gaussian ``.cube`` files.
* **spacing_ang** (*float*): 3D grid spacing in Angstroms (default 0.5 Å).
* **nhomos** / **nlumos** (*int*): Number of frontier spatial/spinor MOs to export.
* **bse_states** (*list of int*): Specific 1-indexed exciton roots to export (hole, electron, and difference densities).

analysis
~~~~~~~~
* **nto** (*bool*): Perform Natural Transition Orbital (NTO) analysis after BSE.
* **nto_states** (*list of int*): Specific exciton states to analyze.
* **nto_top** (*int*): Number of dominant NTO pairs to print per state.
* **nto_csv** (*bool*): Export NTO compactness metrics and weights to ``nto_results.csv``.
* **plot** (*bool*): Generate publication figures and interactive Plotly HTML dashboards.

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
