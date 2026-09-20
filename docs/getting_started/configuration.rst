Configuration Reference (YAML)
==============================

`miniBSE` uses structured YAML files to control calculations. Below is the complete reference of all supported configuration options.

Complete Example Configuration
------------------------------

.. code-block:: yaml

   system:
     basis_txt: "BASIS_MOLOPT_UZH"
     basis_name: "DZVP-MOLOPT-PBE-GTH"
     material: "CSPBBR3"
     nthreads: 12
     device: "auto"                    # "auto", "cpu", "cuda", or "mps"
     gth_file: "GTH_SOC_POTENTIALS.txt" # Required if soc: true

   physics:
     excitation_mode: "diagonal_bse"   # "bse", "diagonal_bse", "independent_qp", "independent_dft"
     qp_gap: "gw"                      # "gw", "brus", or numeric value in eV
     exchange: true                    # Include bare electron-hole exchange
     kernel: "resta"                   # "resta" or "gamma"
     eps_out: 2.4                      # Solvent / external dielectric constant
     nhomos: 1275                      # Active occupied MOs (or all if omitted)
     nlumos: 522                       # Active virtual MOs (or all if omitted)
     soc: true                         # Enable 2-component spinor Spin-Orbit Coupling

   namd:
     trajectory:
       dir: "./trajectory"             # Folder containing frame_000001, frame_000002...
       frame_pattern: "frame_*"
       xyz_file: "frame.xyz"
       mo_file: "MOs.mbse"
       dt_nuc_fs: 2.0                  # Nuclear MD time step in femtoseconds
       start_frame: 1
       end_frame: 500

     storage:
       precompute_dir: "namd_precomputed"
       active_energy_window_ev: [2.0, 7.5] # Retain exciton pairs in [E_min, E_max]
       compress: true

     tracking:
       phase_correction: true          # Fix random U(1) wavefunction signs
       hungarian_tracking: true        # Track state crossings via Hungarian matching
       completeness_threshold: 0.99

     integration:
       integrator: "unitary_matrix_exp" # "unitary_matrix_exp", "strang", "rk4"
       n_substeps: 100                 # Electronic substeps per nuclear time step

     dynamics:
       method: "master_equation"       # "master_equation" (deterministic) or "cpa_fssh"
       n_trajectories: 1000            # Number of trajectories (for cpa_fssh)
       temperature_k: 300.0            # Lattice temperature for detailed balance
       detailed_balance: true
       decoherence: "cumulant"         # "cumulant", "edc", or "none"
       tau_dec_fs: "cumulant"          # "cumulant" (ab initio lowest state) or float in fs
       tau_nr_ns: 50.0                 # Optional non-radiative trap lifetime in ns

     initial_excitation:
       mode: "ratio_eg"                # "ratio_eg" or "energy_ev"
       ratio: 2.0                      # Initial excitation energy = ratio * Eg
       pulse_fwhm_ev: 0.08
       filter_dark_states: true        # Exclude states with zero oscillator strength

     output:
       cooling_curve_csv: "carrier_cooling.csv"
       populations_npz: "populations.npz"
       plot_cooling: true
       plot_file: "carrier_cooling.png"

Section Breakdown
-----------------

system
~~~~~~
* **basis_txt** (*str*): Path to CP2K MOLOPT Gaussian basis set file.
* **basis_name** (*str*): Exact basis set name (e.g. `DZVP-MOLOPT-PBE-GTH`).
* **material** (*str*): Material key in `MATERIAL_DB` (e.g. `CSPBBR3`, `CDSE`, `INAS`, `PBS`).
* **nthreads** (*int*): OpenMP / BLAS CPU threads.
* **device** (*str*): Computation device (`cpu`, `cuda`, `mps`, or `auto`).
* **gth_file** (*str*): Path to GTH pseudopotential file containing spin-orbit parameters.

physics
~~~~~~~
* **excitation_mode** (*str*):
  - `bse`: Full resonant Bethe-Salpeter Equation with configuration mixing.
  - `diagonal_bse`: Includes diagonal bare exchange $K_x$ and screened direct attraction $K_d$ without off-diagonal mixing. Recommended for dense NAMD.
  - `independent_qp`: Single-particle quasiparticle transitions ($\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}$).
  - `independent_dft`: Bare DFT orbital energy differences.
* **qp_gap** (*str or float*):
  - `gw`: Scaled GW model with finite vacuum anchor and dielectric polarization.
  - `brus`: Quantum confinement model.
  - *float*: User-specified target quasiparticle gap in eV.
* **soc** (*bool*): Enables 2-component relativistic spinor Hamiltonian.
* **kernel** (*str*): Dielectric screening model (`resta` or `gamma`).
* **eps_out** (*float*): Surrounding solvent / matrix dielectric constant.

namd
~~~~
* **trajectory.dt_nuc_fs** (*float*): MD time step between frames.
* **dynamics.method** (*str*):
  - `master_equation`: Deterministic Pauli Master Equation with BLAS tensor decomposition. Instant, noise-free, handles millions of states.
  - `cpa_fssh`: Classical Path Approximation Fewest Switches Surface Hopping with stochastic trajectory sampling.
* **dynamics.decoherence** (*str*):
  - `cumulant`: Evaluates second-order cumulant expansion of energy gap fluctuations from lowest excited state.
  - `edc`: Energy Decoherence Correction with kinetic energy scaling.
* **dynamics.tau_nr_ns** (*float, optional*): Defect trap non-radiative lifetime in nanoseconds for PLQY prediction.
