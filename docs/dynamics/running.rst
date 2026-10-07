Running NAMD: workflow, commands and keys
=========================================

Part of :doc:`/dynamics/index`.

Workflow
--------

One YAML file serves all stages; ``namd.storage.precompute_dir`` connects them.

.. code-block:: bash

   qdex --config config.yaml --namd-precompute     # excitons, overlaps, tracking  (sequential over frames)
   qdex --config config.yaml --namd-decoherence    # pair dephasing times          (~1 s)
   qdex --config config.yaml --namd-nac            # logm couplings                (parallel over steps)
   qdex --config config.yaml --namd-run            # FSSH / DISH / PME; repeat for each method, pump, T

The precompute is sequential (each frame is aligned to the previous one); ``--namd-nac`` and the dynamics
are not. On a 192-core node (AMD Genoa, OpenBLAS) for a 4 nm CsPbX\ :sub:`3` cube with 1300 × 800
orbitals:

.. list-table::
   :header-rows: 1

   * - threads
     - 8
     - 16
     - 24
     - 32
     - 48
     - 64
     - 96
   * - spin-free, s/frame
     - 6.0
     - 3.4
     - 2.55
     - 2.0
     - 1.55
     - 1.8
     - 1.6
   * - SOC (2600 × 1600 spinors), s/frame
     - 35.5
     - 27.3
     - 24.5
     - 23.3
     - 22.3
     -
     -

Spin-free saturates at 32–48 threads (7.7 GB); the SOC frame is dominated by the diagonalisation of the
spinor Hamiltonian (4200 × 4200), which OpenBLAS does not parallelise beyond ~16 threads. Several
precomputes therefore share a node well: four chains of 48 threads (two materials, spin-free and SOC) take
about 40 min (spin-free) and 6 h (SOC) for 1000 frames. Disk: about 20 MB per step spin-free and 150 MB
per step with SOC.

Initial conditions and ensembles
--------------------------------

* ``initial_excitation``: ``mode: ratio_eg`` pumps at ``ratio`` × the lowest exciton energy (mean over
  the trajectory), with a Gaussian line of ``pulse_fwhm_ev``; ``filter_dark_states`` weights the pairs by
  their oscillator strength. Excess energy above the gap = (ratio − 1) E\ :sub:`g`. The active window must
  hold the carriers with about 1 eV to spare (:doc:`states_couplings`).
* ``initial_conditions: multiple``: ``n_origins`` runs of ``window_fs`` each, started from frames spread
  over the trajectory (spacing at least 200 fs and twice the gap correlation time); every origin is pumped
  at the same excess energy above its own lowest exciton. ``n_trajectories`` is the total, divided over
  the origins.
* ``seed`` fixes the random numbers (initial sampling, hops, decoherence events).

Command-line flags
------------------

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Flag
     - Default
     - Description
   * - ``--namd-precompute``
     -
     - excitons of every frame, overlaps, tracking (:doc:`states_couplings`)
   * - ``--namd-decoherence [dir]``
     -
     - pair dephasing times into ``decoherence_times.npz`` (:doc:`decoherence`); needed by DISH, the FSSH
       damping and the PME line widths
   * - ``--namd-nac [dir]``
     -
     - couplings of every step into ``nac_<k>_to_<k+1>.npz`` (``namd.nac.scheme``)
   * - ``--namd-run``
     -
     - dynamics from the precomputed data
   * - ``--namd-method <m>``
     - YAML
     - ``cpa_fssh``, ``dish``, ``pme`` (``master_equation``), ``cpa_fssh_gdc``
   * - ``--namd-initial-conditions <mode>``, ``--namd-multi-init``
     - ``single``
     - ``single`` (from frame 0) or ``multiple`` origins
   * - ``--namd-origins <int>``, ``--namd-window-fs <float>``
     - auto
     - number of origins and length of each run (auto: from the gap correlation time and a pilot cooling
       time)
   * - ``--namd-soc``
     -
     - spinors with SOC in the precompute (or ``soc.enabled``)
   * - ``--namd-compact [dir]``
     -
     - compress a precompute directory
   * - ``--namd-ta``, ``--namd-ta-sigma``, ``--namd-ta-plot``
     -
     - transient absorption from the dynamics (:doc:`/spectroscopy/index`)
   * - ``--namd-ecsh-auger``, ``--namd-biexciton``, ``--namd-trajectory-loops``
     -
     - Auger two-body hops, biexciton initial state, looping the trajectory (the wrap from the last frame
       to the first has no coupling)

YAML example
------------

.. code-block:: yaml

   system:
     material: CSPBBR3
     basis_txt: BASIS_MOLOPT_UZH
     basis_name: DZVP-MOLOPT-PBE-GTH
     nthreads: 48
   quasiparticles:
     model: bulk
     bulk_vertex: scaled
     bulk_vertex_factor: material
     bulk_residual: experimental
   excitations:
     mode: diagonal_sbse
     kernel: resta
     include_exchange: false           # K^x: a few meV, 2.7x the cost of a frame
     nhomos: 1300                      # window: ~1 eV beyond the largest carrier excess energy
     nlumos: 800
   soc:
     enabled: false

   namd:
     trajectory:
       dir: ./traj                     # frame_* subdirectories with frame.xyz and MOs.mbse
       frame_pattern: frame_*
       xyz_file: frame.xyz
       mo_file: MOs.mbse               # .mbse, text, or TREXIO HDF5 (geometry read from it)
       dt_nuc_fs: 2.0
       start_frame: 1
       end_frame: 1000
     storage:
       precompute_dir: ./precomputed
     tracking:
       phase_correction: true
       hungarian_tracking: true        # relabel trivial crossings (|S_ii| < 0.5)
       degeneracy_tol_ev: 1.0e-5       # Kramers pairs (SOC): parallel transport inside the pair
     decoherence:
       method: cumulant                # or gaussian (hbar/sigma)
       max_lag_fs: 1000
     nac:
       scheme: logm                    # or hst
       workers: 6
       threads_per_worker: 8
     integration:
       integrator: strang
       n_substeps: 20                  # electronic sub-steps per nuclear step (default 2)
     dynamics:
       method: dish                    # master_equation, dish, cpa_fssh
       temperature_k: 300.0
       detailed_balance: true
       decoherence: cumulant
       nac_scheme: logm
       pme_tau: pairs                  # PME line widths = DISH pair times
       initial_conditions: multiple
       n_origins: 6
       window_fs: 500.0
       n_trajectories: 3000            # total, divided over the origins
       seed: 2026
     initial_excitation:
       mode: ratio_eg
       ratio: 2.0                      # pump at 2 Eg: excess energy 1 Eg
       pulse_fwhm_ev: 0.08
       filter_dark_states: true
     output:
       cooling_curve_csv: cooling.csv
       populations_npz: populations.npz
       plot_cooling: true
       plot_file: cooling.png
