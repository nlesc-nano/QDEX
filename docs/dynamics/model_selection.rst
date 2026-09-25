Model selection
===============

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.precompute``
* Callable: ``qdex.namd.precompute.precompute_namd_data``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   precompute_namd_data(config)


3. Method Selection: PME vs. CPA-FSSH-EDC vs. DISH
--------------------------------------------------

A central methodological decision in non-adiabatic dynamics is choosing among a **deterministic Master Equation (PME)**, **stochastic Fewest Switches Surface Hopping with continuous Energy-Based Decoherence (CPA-FSSH-EDC)**, and **Decoherence-Induced Surface Hopping (DISH)**. All three frameworks are implemented in ``QDEX``, and each possesses distinct physical domains of applicability:

.. list-table::
   :widths: 18 26 28 28
   :header-rows: 1

   * - Criterion / Regime
     - Pauli Master Equation (PME)
     - CPA-FSSH with EDC
     - Decoherence-Induced Surface Hopping (DISH)
   * - **Density of States (DOS)**
     - Dense, quasi-continuous manifolds (:math:`> 10^3 - 10^6` exciton states in large QDs).
     - Discrete or sparse frontier level manifolds (:math:`1S_e, 1P_e` states, small molecules).
     - Dense or intermediate manifolds where single-trajectory hopping events and dwell times are needed.
   * - **Quantum Coherence**
     - Fast dephasing regime: :math:`\tau_{\mathrm{dec}} \ll \tau_{\mathrm{transfer}}`. Coherences decay before population builds up.
     - Coherent regime: quantum interference, state superpositions, and phase memory persist.
     - Fast dephasing with stochastic decoherence-driven wavepacket bifurcation and collapse.
   * - **Phonon Bottleneck**
     - May overestimate relaxation if multi-phonon wavepacket dynamics are approximated by broad rates.
     - **Essential**: captures discrete quantum transitions and coherent vibrational wavepacket motion.
     - **Accurate**: captures rare stochastic transitions and discrete dwell times without Quantum Zeno freezing.
   * - **Surface Traps**
     - Provides average, memoryless Markovian trapping rates; cannot capture stochastic residence times.
     - **Essential**: tracks explicit hopping and de-hopping fluctuations, barrier crossings, and bifurcation.
     - **Essential**: tracks explicit trajectory bifurcation, dwell times, and thermally activated de-trapping.
   * - **Computational Cost**
     - **Ultra-fast**: Tensor decomposition propagates :math:`10^6` states in :math:`< 0.2\text{ s}` per step.
     - **Heavy**: Requires sub-stepping and averaging over :math:`10^3 - 10^4` stochastic trajectories.
     - **Fast Stochastic**: Unitary sub-stepping without continuous damping, with vectorized Poisson branching.


The Phonon Bottleneck Case (:math:`1P_e \to 1S_e`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In quantum-confined semiconductor nanocrystals, quantum confinement shifts atomic orbitals into discrete, shell-like states labeled by atomic-like symmetries (:math:`1S_e, 1P_e, 1D_e`). While higher-lying states form a dense, continuous manifold where cooling is ultrafast, the energy separation between the lowest unoccupied conduction state (:math:`1S_e`) and the next state (:math:`1P_e`) can be several tenths of an electronvolt:

.. math::

   \Delta E(1P_e - 1S_e) \gg \hbar \omega_{\mathrm{LO}}

Because this energy gap greatly exceeds the energy of a single longitudinal optical (LO) phonon (:math:`\hbar \omega_{\mathrm{LO}} \approx 15 - 35\text{ meV}`), relaxation cannot occur through single-phonon scattering. This phenomenon is known as the **phonon bottleneck**.

* **Why FSSH or DISH is Recommended for Bottlenecks**: The :math:`1P_e \to 1S_e` transition is mediated by rare multi-phonon wavepacket coincidences, non-adiabatic surface crossings, or Auger-type electron-hole energy exchange. CPA-FSSH and DISH explicitly evolve the time-dependent Schrödinger equation, capturing quantum interference between discrete electronic states and vibrational wavepackets, and resolving whether the bottleneck persists or is bypassed.


Surface Trap States: Hopping & De-Hopping Kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Colloidal quantum dots frequently contain under-coordinated surface atoms, halide vacancies, or localized ligand termination defects. These defects introduce discrete electronic levels located inside the fundamental band gap (surface traps).

When a hot carrier cools to the band edge:

1. **Carrier Trapping (Hopping into Trap)**: The carrier transitions from a delocalized core state (:math:`1S`) into a spatially localized defect level. This is accompanied by strong local lattice distortion (large polaron or Jahn-Teller rearrangement).
2. **Carrier De-Trapping (Hopping out of Trap)**: Thermal fluctuations from the nuclear bath can impart sufficient energy to kick the carrier back from the defect into the delocalized band states (thermally activated de-trapping).

* **Why Stochastic Methods (FSSH / DISH) are Recommended for Traps**: Trapping and de-trapping are stochastic, trajectory-dependent barrier-crossing events. A deterministic rate equation (PME) treats trapping as an irreversible, memoryless Markovian decay that washes out individual carrier dwell times and trapping/detrapping equilibrium fluctuations. Surface hopping tracks individual stochastic trajectories: some trajectories get trapped permanently, while others hop into the trap, reside there for several picoseconds, and subsequently de-hop back into the band. Capturing this physics accurately requires either **CPA-FSSH-EDC** or **DISH** along **extended AIMD trajectories** (typically :math:`> 10 - 50\text{ ps}`).


Summary Decision Rule
~~~~~~~~~~~~~~~~~~~~~

* Choose **``method: "master_equation"``** when screening carrier cooling lifetimes across dense manifolds in medium-to-large quantum dots (:math:`> 500` atoms) where the density of states is high and fast dephasing dominates.
* Choose **``method: "cpa_fssh"``** (with continuous EDC) when investigating discrete frontier level transitions (:math:`1P \to 1S`), quantum coherence, or strong non-adiabatic coupling spikes near avoided crossings.
* Choose **``method: "dish"``** when simulating dense nanocrystals where Tully's derivative flux causes overcoherence or stiffness, recovering PME cooling rates while preserving single-trajectory stochastic statistics, dwell times, and trap residence kinetics.


