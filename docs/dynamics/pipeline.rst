Pipeline
========

Part of :doc:`/dynamics/index`.

.. figure:: /_static/figures/namd_workflow.svg
   :width: 100%
   :alt: namd workflow

   Non-adiabatic dynamics pipeline and its observables.


.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.precompute``
* Callable: ``qdex.namd.precompute.precompute_namd_data``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   precompute_namd_data(config)


Following the photoexcitation of a semiconductor nanocrystal or quantum dot by an ultrashort laser pulse, high-energy ("hot") electrons and holes rapidly dissipate their excess energy through electron-phonon scattering and non-adiabatic transitions. Carriers cascade down the dense ladder of excited states, cooling toward the band edges before recombining to the ground state.

``QDEX`` features an advanced, high-throughput **Non-Adiabatic Molecular Dynamics (NAMD)** engine designed to simulate carrier relaxation, phonon bottleneck phenomena, surface defect trapping/de-trapping, and photoluminescence recombination along *ab initio* molecular dynamics (AIMD) trajectories. ``--namd-run`` is a cooling calculation. Auger recombination, energy-conserving two-body hops, and a biexciton initial state are off unless ``--auger``, ``--namd-ecsh-auger``, or ``--namd-biexciton`` is set.


1. Overview & The NAMD Pipeline Architecture
--------------------------------------------

The NAMD simulation workflow is decoupled into three modular stages:

.. code-block:: text

   AIMD Trajectory Frames (frame_0001 ... frame_N)
                 │
                 ▼
   [Stage 1: Trajectory Precomputation]
     ├─ Evaluate Quasiparticle & Diagonal BSE Exciton States per frame
     ├─ Compute Cross-Frame Overlaps S(t, t+Δt) analytically via Libint2
     ├─ Extract Non-Adiabatic Couplings (NAC) d_IJ(t) via finite differences
     ├─ Eliminate Random Phase Jumps e^{iθ} via Geometric Phase Alignment
     ├─ Repair trivial crossings only (|S_ii| < 0.5); avoided crossings stay adiabatic
     └─ Compress & Cache Precomputed Data into step_*.npz
                 │
                 ▼
   [Stage 2: Dynamical Propagation]
     ├─ Select Engine: Pauli Master Equation (PME) or CPA-FSSH
     ├─ Compute Ab Initio Cumulant Decoherence from Gap Fluctuations
     ├─ Propagate Electronic Populations & Wavepacket Amplitudes
     └─ Couple to Ground State (Einstein Emission k_rad & Multi-Phonon k_nr)
                 │
                 ▼
   [Stage 3: Deep Analysis & Visualization]
     ├─ Extract Hot Carrier Cooling Lifetimes (τ_elec vs. τ_hole)
     ├─ Compute Phonon Spectral Density J(ω) & Identify Active Modes
     ├─ Correlate NAC vs. Energy Gap (Testing the Energy Gap Law)
     ├─ Track Trap Hopping / De-Hopping Kinetics
     └─ Compute Photoluminescence Quantum Yield (PLQY) & 6-Panel Figure

