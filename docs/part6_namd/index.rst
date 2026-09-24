Part 6: Carrier Cooling Dynamics & Photoluminescence (NAMD)
===========================================================

Following the photoexcitation of a semiconductor nanocrystal or quantum dot by an ultrashort laser pulse, high-energy ("hot") electrons and holes rapidly dissipate their excess energy through electron-phonon scattering and non-adiabatic transitions. Carriers cascade down the dense ladder of excited states, cooling toward the band edges before recombining to the ground state.

``QDEX`` features an advanced, high-throughput **Non-Adiabatic Molecular Dynamics (NAMD)** engine designed to simulate carrier relaxation, phonon bottleneck phenomena, surface defect trapping/de-trapping, and photoluminescence recombination along *ab initio* molecular dynamics (AIMD) trajectories. ``--namd-run`` is a cooling calculation. Auger recombination, energy-conserving two-body hops, and a biexciton initial state are off unless ``--auger``, ``--namd-ecsh-auger``, or ``--namd-biexciton`` is set.

---

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

---

2. Multi-Timescale Integration: Separating Nuclear and Electronic Time Steps
----------------------------------------------------------------------------

A fundamental challenge in simulating non-adiabatic carrier dynamics is the dramatic **timescale mismatch** between nuclear vibrations and electronic phase oscillations.

The Timescale Mismatch
~~~~~~~~~~~~~~~~~~~~~~

* **Nuclear Motion (:math:`\sim 1 - 2\text{ fs}`)**: 
  Atomic nuclei are thousands of times heavier than electrons (:math:`M_{\mathrm{Pb}} / m_e \approx 3.8 \times 10^5`). Nuclear motion is governed by acoustic and optical phonon frequencies (:math:`\omega_{\mathrm{ph}} \approx 50 - 300\text{ cm}^{-1}`), corresponding to vibrational periods of :math:`T_{\mathrm{vib}} \approx 100 - 600\text{ fs}`. A nuclear time step of :math:`\Delta t_{\mathrm{nuc}} \approx 1.0 - 2.0\text{ fs}` is therefore fully sufficient to integrate classical Newton's equations of motion with energy conservation.

* **Electronic Wavefunction Oscillations (:math:`\sim 0.005 - 0.05\text{ fs}`)**:
  In contrast, the electronic wavepacket oscillates at the Bohr transition frequencies:

  .. math::

     \omega_{IJ} = \frac{|E_I - E_J|}{\hbar}

  For electronic energy differences of :math:`\Delta E \approx 1.0 - 4.0\text{ eV}`, the characteristic quantum phase oscillation period is:

  .. math::

     \tau_{\mathrm{elec}} = \frac{2\pi \hbar}{\Delta E} \approx 1.0 - 4.1\text{ fs}

Attempting to propagate the electronic Schrödinger equation using the coarse nuclear step :math:`\Delta t_{\mathrm{nuc}} \sim 1\text{ fs}` violates the Nyquist-Shannon sampling theorem, causing severe numerical instability, catastrophic loss of norm conservation (:math:`\sum_I |c_I|^2 \neq 1`), and unphysical population blowup.

The Classical Path Approximation (CPA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To solve this timescale separation, ``QDEX`` operates within the **Classical Path Approximation (CPA)**. 

In semiconductor nanoclusters, the electronic transition involves one or two electrons out of thousands of valence electrons. To first order, the nuclear trajectory :math:`\mathbf{R}(t)` is driven primarily by the ground-state lattice potential, and the back-reaction of single-carrier relaxation on the heavy nuclear motion is negligible compared to thermal kinetic fluctuations at 300 K.

This provides an immense computational advantage:
1. **Decoupled Workflow**: The heavy *ab initio* DFT molecular dynamics simulation is performed **only once** to generate the classical trajectory :math:`\mathbf{R}(t)`.
2. **Post-Processing Reusability**: All non-adiabatic electronic calculations (FSSH with thousands of stochastic trajectories, or PME at multiple temperatures and decoherence models) are executed in post-processing without ever re-evaluating expensive DFT self-consistent field cycles or nuclear forces.

Electronic Sub-Stepping in CPA-FSSH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Within each nuclear time interval :math:`[t_k, t_{k+1}]` of duration :math:`\Delta t_{\mathrm{nuc}}`, ``QDEX`` divides the interval into :math:`N_{\mathrm{sub}}` fine electronic sub-steps (typically :math:`N_{\mathrm{sub}} = 100 - 500`):

.. math::

   \delta t_{\mathrm{elec}} = \frac{\Delta t_{\mathrm{nuc}}}{N_{\mathrm{sub}}} \approx 0.002 - 0.02\text{ fs}

Along the sub-steps :math:`\tau_m = (m + 0.5) \delta t_{\mathrm{elec}}`, the adiabatic energies are linearly interpolated:

.. math::

   E_I(\tau_m) = E_I(t_k) + \frac{\tau_m}{\Delta t_{\mathrm{nuc}}} \left( E_I(t_{k+1}) - E_I(t_k) \right)

The effective non-adiabatic Hamiltonian driving electronic evolution is constructed in ``qdex.namd.integrator``:

.. math::

   \mathbf{H}_{\mathrm{eff}}(\tau_m) = \operatorname{diag}\left( \mathbf{E}(\tau_m) \right) - i \hbar \, \mathbf{d}(t_k)

Because the non-adiabatic coupling matrix :math:`\mathbf{d}` is anti-Hermitian (:math:`d_{IJ} = -d_{JI}^*`), the product :math:`-i\hbar \mathbf{d}` is **strictly Hermitian**, ensuring that :math:`\mathbf{H}_{\mathrm{eff}}` is Hermitian.

Ensemble propagator
"""""""""""""""""""

The single-wavefunction routine ``step_unitary_matrix_exp`` diagonalizes :math:`\mathbf{H}_{\mathrm{eff}} = \mathbf{V} \boldsymbol{\Lambda} \mathbf{V}^\dagger` and applies the exact unitary exponential

.. math::

   \mathbf{c}(\tau + \delta t_{\mathrm{elec}}) = \mathbf{V} \, \exp\left( -i \boldsymbol{\Lambda} \frac{\delta t_{\mathrm{elec}}}{\hbar} \right) \mathbf{V}^\dagger \, \mathbf{c}(\tau)

which conserves :math:`\sum_I |c_I|^2` to the accuracy of the diagonalization. The ensemble path used by surface hopping does not call that routine. It uses second-order Strang splitting,

.. math::

   \mathbf{U}(\delta t) = e^{-i E \delta t / 2\hbar} \, e^{-\mathbf{d}\, \delta t} \, e^{-i E \delta t / 2\hbar}

with :math:`E` the diagonal-BSE pair energy, binding :math:`K_d` included. Because :math:`\mathbf{d}` is anti-Hermitian, :math:`e^{-\mathbf{d}\,\delta t}` is unitary. Set ``integrator: strang``. The name ``unitary_matrix_exp`` is rejected on the ensemble path.

Tully Hopping Flux Accumulation
"""""""""""""""""""""""""""""""

Tully's fewest switches hopping probabilities are accumulated incrementally across the electronic sub-steps:

.. math::

   g_{I \to J} = \sum_{m=1}^{N_{\mathrm{sub}}} \max\left( 0, \, \frac{2 \delta t_{\mathrm{elec}} \, \operatorname{Re}\left( c_I^*(\tau_m) c_J(\tau_m) d_{IJ} \right)}{|c_I(\tau_m)|^2} \right) \times B_{IJ}(T)

The sum is the hop probability for that nuclear step. If the electron and hole channels together exceed 1, they are scaled so the total is 1 and the event is counted. Each nuclear step of a surface-hopping trajectory starts from the active orbital. The cumulant time :math:`\tau_{\mathrm{dec}}` is computed and printed. It is not applied again inside a step that already begins in a pure active state. Upward one-body hops are multiplied by :math:`B_{IJ}` after the flux is accumulated. That factor is the classical-path stand-in for a rejected velocity rescaling (Parandekar and Tully, J. Chem. Phys. 2005; Jain, Alguire, and Subotnik, J. Chem. Phys. 2016). Two-body Auger hops, when requested, do not carry it.

Electronic Sub-Stepping in the Pauli Master Equation (PME)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although the Pauli Master Equation propagates real-valued populations :math:`P_I(t)` rather than oscillating complex amplitudes :math:`c_I(t)`, sub-stepping is equally vital for numerical stability.

In dense manifolds where non-adiabatic couplings are large, single-step forward Euler integration of :math:`\frac{d\mathbf{P}}{dt} = \mathbf{R} \mathbf{P}` with a coarse step :math:`\Delta t_{\mathrm{nuc}} \sim 1\text{ fs}` can lead to **stiffness instabilities**, causing populations to oscillate or become negative (:math:`P_I < 0`).

In ``qdex.namd.master_equation``, two robust solutions are provided:

1. **Exact Matrix Exponential**:
   For modest state spaces, the population vector is propagated analytically over the nuclear step:

   .. math::

      \mathbf{P}(t + \Delta t_{\mathrm{nuc}}) = \exp\left( \mathbf{R} \, \Delta t_{\mathrm{nuc}} \right) \mathbf{P}(t)

2. **Tensorized row-stochastic step (Diagonal BSE)**:
   For the pair manifold actually propagated, each sub-step builds a row-stochastic matrix :math:`T = I + \delta t K` from the electron channel and from the hole channel. A row whose leaving probability would exceed 1 is renormalized onto its outgoing transitions. Applying :math:`\mathbf{P} \leftarrow \mathbf{P} T_e` and then :math:`\mathbf{P} \leftarrow T_h^\mathsf{T} \mathbf{P}` keeps every entry non-negative and conserves probability. Recombination :math:`\exp(-k_{\mathrm{loss}}\Delta t)` is applied once after the sub-steps, without putting that lost population back.

---

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


---

4. Theoretical Foundations of the Dynamical Engines
---------------------------------------------------

1. Derivation of the Pauli Master Equation (PME)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The exact quantum dynamics of a coupled electron-nuclear system is governed by the Liouville-von Neumann equation for the total density operator :math:`\hat{\rho}(t)`:

.. math::

   i\hbar \frac{\partial \hat{\rho}(t)}{\partial t} = \left[ \hat{H}(t), \hat{\rho}(t) \right]

Applying the **Nakajima-Zwanzig projection operator technique**, the density matrix is partitioned into diagonal populations :math:`\mathcal{P} \hat{\rho} = \sum_I \rho_{II} |I\rangle \langle I|` and off-diagonal electronic coherences :math:`\mathcal{Q} \hat{\rho} = \sum_{I \neq J} \rho_{IJ} |I\rangle \langle J|`.

In a condensed-phase environment containing many nuclear degrees of freedom, thermal fluctuations of the nuclear bath induce rapid random phase fluctuations that cause the off-diagonal coherences :math:`\rho_{IJ}(t)` to decay exponentially with characteristic decoherence time :math:`\tau_{\mathrm{dec}}`:

.. math::

   \rho_{IJ}(t) \approx \rho_{IJ}(0) \, \exp\left( -\frac{t}{\tau_{\mathrm{dec}}} \right) \exp\left( -i \frac{\Delta E_{IJ}}{\hbar} t \right)

Under the **Markovian approximation** (where the bath correlation time is much shorter than the population relaxation time), integrating out the rapidly decaying coherences yields the closed **Pauli Master Equation** for populations :math:`P_I(t) \equiv \rho_{II}(t)`:

.. math::

   \frac{d P_I(t)}{dt} = \sum_{J \neq I} \left[ k_{J \to I}(t) P_J(t) - k_{I \to J}(t) P_I(t) \right]

Fermi's Golden Rule with Lorentzian Broadening
""""""""""""""""""""""""""""""""""""""""""""""

The instantaneous rate constant :math:`k_{I \to J}` is obtained from second-order time-dependent perturbation theory (Fermi's Golden Rule). The standard energy-conserving Dirac delta function :math:`\delta(E_I - E_J)` is broadened by the finite electronic decoherence time :math:`\tau_{\mathrm{dec}}`, yielding a Lorentzian line shape:

.. math::

   k_{I \to J}(t) = 2 |d_{IJ}(t)|^2 \, \left[ \frac{\tau_{\mathrm{dec}}}{1 + \left( \frac{(E_J(t) - E_I(t)) \tau_{\mathrm{dec}}}{\hbar} \right)^2} \right] \times B_{IJ}(T)

where :math:`d_{IJ}(t) = \langle \psi_I | \frac{\partial}{\partial t} | \psi_J \rangle` is the non-adiabatic coupling, and :math:`B_{IJ}(T)` enforces thermodynamic **detailed balance** at lattice temperature :math:`T`:

.. math::

   B_{IJ}(T) = \begin{cases}
     1 & \text{for downward transitions } (E_J \le E_I) \\
     \exp\left( -\frac{E_J - E_I}{k_B T} \right) & \text{for upward thermal activation } (E_J > E_I)
   \end{cases}

Vectorized Tensor Decomposition for Diagonal BSE
""""""""""""""""""""""""""""""""""""""""""""""""

In a two-particle excitonic manifold with :math:`N_{\mathrm{occ}}` occupied orbitals and :math:`N_{\mathrm{virt}}` virtual orbitals, the total number of electron-hole pairs is :math:`N_{\mathrm{pairs}} = N_{\mathrm{occ}} \times N_{\mathrm{virt}}`. Constructing and multiplying an :math:`(N_{\mathrm{pairs}} \times N_{\mathrm{pairs}})` rate matrix scales as :math:`O(N_{\mathrm{pairs}}^2) = O(N_{\mathrm{occ}}^2 N_{\mathrm{virt}}^2)`. For :math:`N_{\mathrm{occ}} = N_{\mathrm{virt}} = 500`, this corresponds to an intractable :math:`250,000 \times 250,000` dense matrix (:math:`500\text{ GB}` of RAM).

Under the **Diagonal BSE** representation, the exciton state :math:`|ia\rangle` factorizes into an independent occupied hole state :math:`i` and an independent virtual electron state :math:`a`. An exciton relaxes either via an electron transition (:math:`a \to b`) with rate :math:`K_e(a \to b)` or a hole transition (:math:`i \to j`) with rate :math:`K_h(i \to j)`:

.. math::

   \frac{d P_{ia}(t)}{dt} = \sum_{b \neq a} \left[ K_e(b \to a) P_{ib} - K_e(a \to b) P_{ia} \right] + \sum_{j \neq i} \left[ K_h(j \to i) P_{ja} - K_h(i \to j) P_{ia} \right]

In matrix notation, this decomposes into an exact **BLAS Level-3 tensor product**:

.. math::

   \frac{\partial \mathbf{P}}{\partial t} = \left( \mathbf{P} \, \mathbf{K}_e - \mathbf{P} \operatorname{diag}(\mathbf{L}_e) \right) + \left( \mathbf{K}_h^\mathsf{T} \, \mathbf{P} - \operatorname{diag}(\mathbf{L}_h) \, \mathbf{P} \right)

where :math:`\mathbf{L}_e = \sum_b K_e(a \to b)` and :math:`\mathbf{L}_h = \sum_j K_h(i \to j)` are the total state loss vectors.

This breakthrough reduces the computational scaling from :math:`O(N_{\mathrm{pairs}}^2)` to :math:`O(N_{\mathrm{occ}}^2 + N_{\mathrm{virt}}^2)`. A million exciton configurations are propagated in **less than 0.2 seconds per nuclear time step**.

2. Classical Path Approximation Surface Hopping with Energy-Based Decoherence (CPA-FSSH-EDC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In CPA-FSSH, the classical nuclei follow a precomputed ground-state molecular dynamics trajectory :math:`\mathbf{R}(t)`. The quantum electronic wavefunction :math:`|\Psi(t)\rangle` evolves according to the Time-Dependent Schrödinger Equation (TDSE):

.. math::

   i\hbar \frac{\partial |\Psi(t)\rangle}{\partial t} = \hat{H}_{\mathrm{elec}}(\mathbf{r}; \mathbf{R}(t)) |\Psi(t)\rangle

Expanding the electronic state in the instantaneous adiabatic Kohn-Sham or diagonal-BSE pair basis :math:`|\psi_I(\mathbf{R}(t))\rangle`:

.. math::

   |\Psi(t)\rangle = \sum_I c_I(t) |\psi_I(\mathbf{R}(t))\rangle

Substituting this expansion into the TDSE and projecting onto :math:`\langle \psi_I|` yields the coupled equations of motion for the complex amplitudes :math:`c_I(t)`:

.. math::

   i\hbar \frac{d c_I(t)}{dt} = E_I(t) c_I(t) - i\hbar \sum_J d_{IJ}(t) c_J(t)

where :math:`E_I(t) = \langle \psi_I | \hat{H}_{\mathrm{elec}} | \psi_I \rangle` is the instantaneous adiabatic energy and :math:`d_{IJ}(t) = \langle \psi_I | \frac{\partial}{\partial t} | \psi_J \rangle` is the non-adiabatic coupling matrix element.

In matrix notation, this forms the effective non-adiabatic Schrödinger equation:

.. math::

   \frac{d \mathbf{c}(t)}{dt} = -\frac{i}{\hbar} \mathbf{H}_{\mathrm{eff}}(t) \mathbf{c}(t), \qquad \mathbf{H}_{\mathrm{eff}}(t) = \operatorname{diag}\left(\mathbf{E}(t)\right) - i\hbar \mathbf{d}(t)

Because :math:`\mathbf{d}` is anti-Hermitian (:math:`d_{IJ} = -d_{JI}^*`), the off-diagonal coupling term :math:`-i\hbar \mathbf{d}` is **strictly Hermitian**, ensuring that :math:`\mathbf{H}_{\mathrm{eff}}` is Hermitian and the total electronic probability is conserved:

.. math::

   \sum_I |c_I(t)|^2 = 1

Unitary Time Evolution via Second-Order Strang Splitting
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To propagate the electronic coefficients across the nuclear time step :math:`[t_k, t_{k+1}]`, ``QDEX`` divides the interval into fine sub-steps :math:`\delta t = \Delta t_{\mathrm{nuc}} / N_{\mathrm{sub}}`. Within each sub-step, the propagator is evaluated using second-order Strang splitting:

.. math::

   \mathbf{U}(\delta t) = \exp\left( -i \mathbf{E}(\tau_m) \frac{\delta t}{2\hbar} \right) \exp\left( -\mathbf{d}(t_k) \, \delta t \right) \exp\left( -i \mathbf{E}(\tau_m) \frac{\delta t}{2\hbar} \right)

where :math:`\mathbf{E}(\tau_m)` is linearly interpolated between :math:`\mathbf{E}(t_k)` and :math:`\mathbf{E}(t_{k+1})`. Because :math:`\mathbf{d}` is anti-Hermitian, :math:`\exp(-\mathbf{d}\,\delta t)` is an exact unitary rotation, preserving norm conservation to machine precision.

Tully's Fewest Switches Hopping Probability
"""""""""""""""""""""""""""""""""""""""""""

In Fewest Switches Surface Hopping (Tully, *J. Chem. Phys.* 93, 1061, 1990), an ensemble of :math:`N_{\mathrm{traj}}` classical trajectories is propagated. Each trajectory :math:`tr` resides on a specific "active surface" :math:`K(t)`.

At each time step, the probability for a trajectory on active surface :math:`K` to switch to an inactive surface :math:`J \neq K` is determined by the outward probability flux accumulated across the electronic sub-steps:

.. math::

   g_{K \to J} = \sum_{m=1}^{N_{\mathrm{sub}}} \max\left( 0, \, \frac{2 \delta t \, \operatorname{Re}\left( c_K^*(\tau_m) c_J(\tau_m) d_{KJ}(t_k) \right)}{|c_K(\tau_m)|^2} \right) \times B_{KJ}(T)

where :math:`B_{KJ}(T)` enforces detailed balance at lattice temperature :math:`T`:

.. math::

   B_{KJ}(T) = \min\left( 1, \, \exp\left( -\frac{\max(E_J - E_K, 0)}{k_B T} \right) \right)

In standard molecular dynamics, detailed balance is maintained through momentum rescaling along the non-adiabatic coupling vector :math:`\mathbf{d}_{KJ}` (rejecting upward hops if nuclear kinetic energy is insufficient). In the Classical Path Approximation, nuclear velocities follow a fixed ground-state trajectory; multiplying upward hops by :math:`B_{KJ}(T)` provides the rigorous classical-path counterpart to velocity rescaling (Parandekar & Tully, *J. Chem. Phys.* 122, 094102, 2005; Jain, Alguire, & Subotnik, *J. Chem. Phys.* 144, 214110, 2016).

Hopping Decision & Wavepacket Collapse
""""""""""""""""""""""""""""""""""""""

A uniform random number :math:`\xi \in [0, 1)` is generated:
* If :math:`\sum_{L=1}^{J-1} g_{K \to L} < \xi \le \sum_{L=1}^J g_{K \to L}`, trajectory :math:`tr` hops to state :math:`J`.
* The active surface is reassigned: :math:`K \leftarrow J`.
* The electronic wavepacket undergoes **projective collapse** onto the new active state:

  .. math::

     c_J \leftarrow 1.0, \qquad c_{L \neq J} \leftarrow 0.0

This projective collapse models the immediate decoherence of the trajectory as it separates along the newly occupied potential surface.

The Overcoherence Problem & Continuous Energy-Based Decoherence (EDC)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In standard FSSH, when a trajectory does **not** hop, electronic amplitudes :math:`c_J` remain populated indefinitely because classical nuclei in the CPA follow a shared trajectory. In dense quantum-dot manifolds with hundreds of states, this causes artificial state superpositions to linger, producing unphysical hopping loops (**overcoherence**).

Conversely, if the wavepacket is reset to a pure state at every nuclear step (:math:`c_K = 1, c_{J \neq K} = 0`), frequent projective measurements collapse the wavepacket before transitions can develop, completely freezing carrier relaxation—a pathological numerical artifact known as the **Quantum Zeno effect**.

To resolve both problems, ``QDEX`` implements the **Energy-Based Decoherence (EDC)** scheme (Granucci & Persico, *J. Chem. Phys.* 126, 134114, 2007). 

In EDC, electronic wavepacket amplitudes persist continuously across nuclear steps. When no hop occurs, inactive states :math:`J \neq K` are continuously damped over the nuclear step :math:`\Delta t`:

.. math::

   c_J(t + \Delta t) \leftarrow c_J(t + \Delta t) \, \exp\left( -\frac{\Delta t}{\tau_{KJ}} \right), \qquad \forall J \neq K

where :math:`\tau_{KJ}` is the state-pair pure-dephasing time. To strictly preserve total probability conservation (:math:`\sum_L |c_L|^2 = 1`) without altering the quantum phase of the active state, the active amplitude :math:`c_K` is renormalized:

.. math::

   c_K(t + \Delta t) \leftarrow \sqrt{ 1 - \sum_{J \neq K} |c_J(t + \Delta t)|^2 } \; \frac{c_K(t + \Delta t)}{|c_K(t + \Delta t)|}

This continuous damping smoothly suppresses off-diagonal coherences on physical dephasing timescales (:math:`10 - 40\text{ fs}`) while allowing short-time quantum interference to drive hopping.

State-Pair Dephasing Times: Ab Initio Covariance vs. Granucci-Persico
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``QDEX`` supports two rigorous evaluations of :math:`\tau_{KJ}`:

1. **Ab Initio Pure-Dephasing Covariance Matrix**:
   Evaluated from the trajectory energy fluctuations:

   .. math::

      \tau_{KJ} = \min\left( \tau_{\max}, \, \frac{\hbar \sqrt{2}}{\sigma_{KJ}} \right), \qquad \sigma_{KJ}^2 = \operatorname{Var}\left( E_K(t) - E_J(t) \right)

   Diagonal elements are capped at :math:`\tau_{\max} = 500\text{ fs}`. Precomputed across all orbital pairs and cached in ``decoherence_times.npz``.

2. **Granucci-Persico Energy-Gap Formula**:
   Evaluated instantaneously from the energy gap and nuclear kinetic energy:

   .. math::

      \tau_{KJ} = \frac{\hbar}{|E_K - E_J|} \left( 1 + \frac{C}{E_{\mathrm{kin}}} \right)

   where :math:`C = 0.1\text{ Hartree} = 2.721\text{ eV}`, and :math:`E_{\mathrm{kin}} = \frac{3}{2} N_{\mathrm{atoms}} k_B T` is the classical thermal kinetic energy of the nuclear lattice.

---

3. Decoherence-Induced Surface Hopping (DISH)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While CPA-FSSH with EDC uses Tully's derivative coupling flux to drive transitions and applies decoherence as an extrinsic damping correction, **Decoherence-Induced Surface Hopping (DISH)** (Jaeger, Fischer, & Prezhdo, *J. Chem. Phys.* 137, 22A545, 2012; Akimov & Prezhdo, *J. Chem. Phys.* 138, 124102, 2013) introduces a fundamentally different physical paradigm.

Physical Foundations of DISH
""""""""""""""""""""""""""""

In condensed matter systems (colloidal quantum dots, perovskite nanocrystals, organic semiconductors), an electronic excitation couples to thousands of nuclear vibrational degrees of freedom. Thermal phonon fluctuations destroy electronic phase coherence within :math:`5 - 25\text{ fs}`.

In this fast-dephasing regime, electronic transitions are **not driven by instantaneous derivative coupling spikes**, but rather by **environment-induced decoherence (wavepacket branching into the bath)**. DISH operationalizes this insight by formulating surface hopping directly in terms of quantum measurement theory and stochastic wavepacket collapse.

Piecewise Unitary Electronic Propagation
""""""""""""""""""""""""""""""""""""""""

Between stochastic collapse events, the electronic wavepacket :math:`\mathbf{c}(t)` evolves **strictly unitarily** according to the Time-Dependent Schrödinger Equation:

.. math::

   \mathbf{c}(t + \Delta t) = \mathbf{U}(t, t + \Delta t) \, \mathbf{c}(t)

Crucially, **no artificial continuous exponential damping** is applied to :math:`\mathbf{c}(t)` during unitary propagation. Quantum superpositions and phase interference evolve naturally.

Step 1: Poisson Stochastic Dephasing
""""""""""""""""""""""""""""""""""""

At each nuclear step :math:`\Delta t`, for every inactive state :math:`J \neq K` (where :math:`K` is the active surface of trajectory :math:`tr`), the occurrence of a decoherence event is governed by a Poisson arrival process:

.. math::

   P_{\mathrm{dec}, J} = 1 - \exp\left( -\frac{\Delta t}{\tau_{KJ}} \right)

where :math:`\tau_{KJ}` is the state-pair pure-dephasing time (loaded from ``decoherence_times.npz``).

For each state :math:`J \neq K`, a uniform random number :math:`R_1 \in [0, 1)` is sampled:
* If :math:`R_1 \ge P_{\mathrm{dec}, J}`: State :math:`J` remains coherent with active state :math:`K`. No collapse attempt is made for state :math:`J`.
* If :math:`R_1 < P_{\mathrm{dec}, J}`: A dephasing event has occurred. The nuclear wavepacket associated with state :math:`J` has spatially separated from the wavepacket on surface :math:`K`. The trajectory must now undergo stochastic branching!

Step 2: Stochastic Branching (Collapse vs. Quenching)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When state :math:`J` dephases, the trajectory reaches a quantum bifurcation point. In accordance with the Born rule, the probability that the system collapses into state :math:`J` is given by its instantaneous electronic population :math:`|c_J|^2`, scaled by detailed balance:



.. math::

   P_{\mathrm{hop}, J} = |c_J|^2 \times \min\left( 1, \, \exp\left( -\frac{\max(E_J - E_K, 0)}{k_B T} \right) \right)

A second independent uniform random number :math:`R_2 \in [0, 1)` is sampled:

* **Case A: Hop Accepted (:math:`R_2 < P_{\mathrm{hop}, J}`)**:
  The trajectory successfully transitions to state :math:`J`. The active surface switches :math:`K \leftarrow J`, and the wavepacket undergoes complete projective collapse onto state :math:`J`:

  .. math::

     c_J \leftarrow 1.0, \qquad c_{L \neq J} \leftarrow 0.0

  If multiple states simultaneously qualify for a hop in a single time step, one state is selected with probability proportional to :math:`P_{\mathrm{hop}, J}`.

* **Case B: Hop Rejected (:math:`R_2 \ge P_{\mathrm{hop}, J}`)**:
  The trajectory remains on active surface :math:`K`. Because a dephasing event did occur, coherence between state :math:`J` and active state :math:`K` has been irreversibly lost to the nuclear bath. Consequently, state :math:`J` is **quenched**:

  .. math::

     c_J \leftarrow 0.0

  The remaining surviving amplitudes are renormalized to conserve total probability:

  .. math::

     \mathbf{c} \leftarrow \frac{\mathbf{c}}{\sqrt{\sum_L |c_L|^2}}

This stochastic quenching removes off-diagonal population without continuous damping, completely preventing overcoherence and avoiding the Quantum Zeno effect.

Analytical Equivalence of DISH to the Pauli Master Equation
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A profound theoretical property of DISH is that in the condensed-phase limit where dephasing is fast compared to electronic transitions (:math:`\tau_{KJ} \ll \tau_{\mathrm{transfer}}`), DISH **analytically converges to the Pauli Master Equation with Lorentzian line broadening**!

*Proof Sketch*:
Consider a two-level system initialized in active state :math:`K` (:math:`c_K(0) = 1`, :math:`c_J(0) = 0`). Over a short nuclear interval :math:`\Delta t`, first-order perturbation theory on the TDSE gives:

.. math::

   c_J(\Delta t) \approx - \frac{d_{KJ} \, \Delta t}{1 + i \frac{\Delta E_{KJ} \Delta t}{2\hbar}}

The population amplitude built up in state :math:`J` during interval :math:`\Delta t` is:

.. math::

   |c_J(\Delta t)|^2 \approx \frac{|d_{KJ}|^2 \Delta t^2}{1 + \left( \frac{\Delta E_{KJ} \Delta t}{2\hbar} \right)^2}

In DISH, the transition rate :math:`k_{K \to J}^{\mathrm{DISH}}` is the product of the dephasing frequency :math:`\Gamma_{\mathrm{dec}} = 1/\tau_{KJ}` and the branching probability :math:`P_{\mathrm{hop}, J} \approx |c_J|^2`:

.. math::

   k_{K \to J}^{\mathrm{DISH}} = \frac{P_{\mathrm{dec}, J} \cdot P_{\mathrm{hop}, J}}{\Delta t} \approx \frac{1}{\tau_{KJ}} \left[ \frac{|d_{KJ}|^2 \tau_{KJ}^2}{1 + \left( \frac{\Delta E_{KJ} \tau_{KJ}}{\hbar} \right)^2} \right] = |d_{KJ}|^2 \left[ \frac{\tau_{KJ}}{1 + \left( \frac{\Delta E_{KJ} \tau_{KJ}}{\hbar} \right)^2} \right]

Apart from a standard factor of 2 arising from the full two-sided integration of the bath autocorrelation function, this expression is **mathematically identical to Fermi's Golden Rule with Lorentzian broadening** used in the Pauli Master Equation:

.. math::

   k_{K \to J}^{\mathrm{PME}} = 2 |d_{KJ}|^2 \left[ \frac{\tau_{\mathrm{dec}}}{1 + \left( \frac{\Delta E_{KJ} \tau_{\mathrm{dec}}}{\hbar} \right)^2} \right] \times B_{KJ}(T)

**Conclusion**: DISH provides a rigorous theoretical unification of wavepacket quantum dynamics and statistical master equations. It recovers PME cooling rates in dense manifolds while retaining single-trajectory stochastic statistics, individual dwell times, and branching kinetics.

---

4. Comparative Synthesis: PME vs. CPA-FSSH-EDC vs. DISH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following comprehensive comparison synthesizes the mathematical foundations, computational characteristics, and physical applicability of the three dynamical engines in ``QDEX``:

.. list-table::
   :widths: 20 26 27 27
   :header-rows: 1

   * - Theoretical Dimension
     - Pauli Master Equation (PME)
     - CPA-FSSH with EDC
     - Decoherence-Induced Surface Hopping (DISH)
   * - **Primary Dynamical Variable**
     - Real-valued population vector :math:`\mathbf{P}(t) \in \mathbb{R}^{N}`.
     - Complex amplitude vector :math:`\mathbf{c}(t) \in \mathbb{C}^{N}` plus active surface index :math:`K(t)`.
     - Complex amplitude vector :math:`\mathbf{c}(t) \in \mathbb{C}^{N}` plus active surface index :math:`K(t)`.
   * - **Electronic Equation of Motion**
     - Markovian rate equation: :math:`\frac{d\mathbf{P}}{dt} = \mathbf{K} \mathbf{P}`.
     - Continuous TDSE: :math:`i\hbar \dot{\mathbf{c}} = \mathbf{H}_{\mathrm{eff}} \mathbf{c}` with continuous EDC damping.
     - Piecewise unitary TDSE: :math:`\mathbf{c}(t+\Delta t) = \mathbf{U} \mathbf{c}(t)` (no continuous damping).
   * - **Hopping / Transition Mechanism**
     - Continuous probability flux between populations via FGR rate matrix.
     - Tully's Fewest Switches non-adiabatic coupling flux: :math:`g_{K \to J} \propto \operatorname{Re}(c_K^* c_J d_{KJ})`.
     - Environment-driven stochastic branching: Poisson dephasing followed by Born-rule collapse :math:`P_{\mathrm{hop}} \propto |c_J|^2`.
   * - **Decoherence Treatment**
     - Implicit: Lorentzian broadening of energy conservation delta function.
     - Continuous exponential damping of inactive states: :math:`c_J \leftarrow c_J e^{-\Delta t / \tau_{KJ}}`.
     - Discrete stochastic quenching of inactive states: :math:`c_J \leftarrow 0` upon rejected dephasing.
   * - **Detailed Balance Enforcement**
     - Transition rates scaled by Boltzmann factor :math:`B_{IJ}(T) = \exp(-\Delta E / k_B T)`.
     - Hopping fluxes scaled by Boltzmann factor :math:`B_{IJ}(T)` (classical-path velocity rescaling stand-in).
     - Hopping probabilities scaled by Boltzmann factor :math:`B_{IJ}(T)`.
   * - **Quantum Zeno Vulnerability**
     - **Immune**: Propagates macroscopic populations; no wavepacket resetting.
     - **Immune**: Wavepacket amplitudes persist across steps; EDC damps smoothly without abrupt collapse.
     - **Immune**: Decoherence events are stochastic and Poisson-distributed; active state remains untouched.
   * - **Equivalence in Dense Limit**
     - Canonical baseline (Lorentzian FGR).
     - Diverges if dephasing times are miscalibrated or couplings are stiff.
     - **Analytically proves equivalence** to Lorentzian FGR in the fast-dephasing limit.
   * - **Computational Complexity**
     - :math:`O(N_{\mathrm{occ}}^2 + N_{\mathrm{virt}}^2)` via BLAS-3 tensor contraction (< 0.2 s per step for :math:`10^6` pairs).
     - :math:`O(N_{\mathrm{traj}} \cdot N_{\mathrm{sub}} \cdot N_{\mathrm{dyn}}^2)` (heavy; requires fine sub-stepping).
     - :math:`O(N_{\mathrm{traj}} \cdot N_{\mathrm{sub}} \cdot N_{\mathrm{dyn}}^2)` (vectorized batching across trajectories).
   * - **Recommended Regime**
     - Ultrafast screening across huge state manifolds (:math:`> 1,000` states), high DOS.
     - Discrete frontier states (:math:`1P \to 1S`), persistent quantum coherence, small systems.
     - Nanocrystals with dense bands, avoided crossings, defect trapping, and single-carrier dwell times.


---

5. Trajectory Precomputation & Wavefunction Tracking
-----------------------------------------------------

Numerical Non-Adiabatic Couplings (NAC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Along the classical nuclear trajectory :math:`\mathbf{R}(t)`, the non-adiabatic coupling is evaluated numerically via finite differences:

.. math::

   d_{IJ}(t + \frac{\Delta t}{2}) = \langle \psi_I(t) | \frac{\partial}{\partial t} | \psi_J(t) \rangle \approx \frac{S_{IJ}(t, t+\Delta t) - S_{JI}^*(t, t+\Delta t)}{2 \Delta t}

where :math:`S_{IJ}(t, t+\Delta t) = \langle \psi_I(t) | \psi_J(t+\Delta t) \rangle` is the cross-frame state overlap. In ``QDEX``, the underlying atomic orbital cross-overlaps :math:`S_{\mu \nu}(t, t+\Delta t) = \int \chi_\mu(\mathbf{r}; \mathbf{R}(t)) \chi_\nu(\mathbf{r}; \mathbf{R}(t+\Delta t)) d\mathbf{r}` are evaluated analytically via Libint2.

Eliminating Gauge Phase Discontinuities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard electronic eigensolvers determine eigenvectors up to an arbitrary global phase factor :math:`e^{i \theta_I(t)}`. If left uncorrected, random sign flips between successive MD frames cause :math:`S_{II}(t, t+\Delta t) \approx -1`, producing spurious non-adiabatic couplings that are orders of magnitude too large.

``QDEX`` eliminates gauge discontinuities by applying a phase rotation:

.. math::

   |\psi_I(t+\Delta t)\rangle \leftarrow |\psi_I(t+\Delta t)\rangle \, e^{-i \theta_I}

where :math:`\theta_I = \operatorname{arg}(S_{II}(t, t+\Delta t))`. This guarantees that the diagonal overlap is strictly real and positive:

.. math::

   \operatorname{Re}(S_{II}(t, t+\Delta t)) \ge 0, \quad \operatorname{Im}(S_{II}(t, t+\Delta t)) = 0

Hungarian Matching for Trivial Crossings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A trivial crossing is a label swap between two orbitals that do not interact. Inside one nuclear step their diagonal overlap collapses and the character sits on an off-diagonal element. Following the old energy label then puts a numerical spike into :math:`d_{IJ}`.

``QDEX`` repairs only that case. The Hungarian assignment uses the cost matrix

.. math::

   C_{IJ} = 1 - |S_{IJ}(t, t+\Delta t)|^2

and a state is allowed to leave its own column only when :math:`|S_{II}| < 0.5`. An avoided crossing that still overlaps its own adiabatic label stays in the adiabatic basis that surface hopping propagates. The assignment is not a global diabatization.

---

6. Electronic Decoherence: Origin, Computation, and Rationale
-------------------------------------------------------------

Physical Origin of Electronic Decoherence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In true quantum dynamics, an electronic superposition state :math:`|\Psi\rangle = c_1 |\psi_1\rangle + c_2 |\psi_2\rangle` is accompanied by nuclear wavepackets moving on the respective potential energy surfaces:

.. math::

   |\Psi_{\mathrm{total}}(t)\rangle = c_1 |\psi_1\rangle |\chi_1(t)\rangle + c_2 |\psi_2\rangle |\chi_2(t)\rangle

Because surfaces 1 and 2 exert different forces (:math:`-\boldsymbol{\nabla} E_1 \neq -\boldsymbol{\nabla} E_2`), the nuclear wavepackets :math:`|\chi_1(t)\rangle` and :math:`|\chi_2(t)\rangle` accelerate differently and rapidly separate in nuclear configuration space.

The electronic coherence is proportional to the nuclear wavepacket overlap:

.. math::

   \rho_{12}(t) \propto c_1 c_2^* \, \langle \chi_2(t) | \chi_1(t) \rangle

As soon as the wavepackets separate spatially (:math:`\langle \chi_2 | \chi_1 \rangle \to 0`), electronic coherence is destroyed. This process is called **electronic decoherence** or **quantum dephasing**.

Why Dephasing is Ultrafast in Nanocrystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A quantum dot possesses thousands of vibrational normal modes. Each phonon mode slightly modulates the electronic transition energy. Because these fluctuations are largely uncorrelated, their destructive phase interference leads to **ultrafast dephasing within 5 to 25 femtoseconds**. Once dephased, the system behaves as a statistical mixture of classical probabilities, justifying the Master Equation.

Automated Ab Initio Cumulant Decoherence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than guessing an empirical decoherence time (e.g. :math:`\tau_{\mathrm{dec}} = 10\text{ fs}`), ``QDEX`` computes :math:`\tau_{\mathrm{dec}}` *ab initio* from the second-order cumulant expansion of energy gap fluctuations.

Let :math:`\delta E_1(t) = E_1(t) - \langle E_1 \rangle` be the instantaneous fluctuation of the lowest excited state gap along the trajectory. The normalized gap autocorrelation function is:

.. math::

   C(t) = \frac{\langle \delta E_1(0) \delta E_1(t) \rangle}{\sigma_E^2}

where :math:`\sigma_E^2 = \langle \delta E_1^2 \rangle` is the variance.

The bath line shape function :math:`g(t)` is evaluated by double integration:

.. math::

   g(t) = \frac{\sigma_E^2}{\hbar^2} \int_0^t dt_1 \int_0^{t_1} dt_2 \, C(t_2)

The electronic dephasing decay function is:

.. math::

   D(t) = \exp\left( -g(t) \right)

``QDEX`` solves for the characteristic decoherence time :math:`\tau_{\mathrm{dec}}` satisfying:

.. math::

   D(\tau_{\mathrm{dec}}) = \frac{1}{e}

Why We Compute It This Way
""""""""""""""""""""""""""

1. **Parameter-Free**: Eliminates arbitrary empirical fitting parameters from NAMD simulations.
2. **Temperature & Lattice Sensitive**: Soft, anharmonic lattices (such as lead halide perovskites) exhibit large thermal gap fluctuations (:math:`\sigma_E \approx 50 - 100\text{ meV}`), correctly yielding short dephasing times (:math:`\tau_{\mathrm{dec}} \approx 7 - 12\text{ fs}`), whereas rigid covalent quantum dots (like InAs or Si) yield longer dephasing times (:math:`\tau_{\mathrm{dec}} \approx 20 - 40\text{ fs}`).

State-Pair Pure-Dephasing Matrices (tau_ij)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In multi-state quantum dynamics, electronic dephasing times vary substantially across orbital pairs depending on energy gaps, orbital spatial localization, and lattice phonon coupling. Rather than assuming an arbitrary uniform dephasing time, ``QDEX`` precomputes the complete state-pair pure-dephasing matrices :math:`\boldsymbol{\tau}_{\mathrm{occ}} \in \mathbb{R}^{N_{\mathrm{occ}} \times N_{\mathrm{occ}}}` and :math:`\boldsymbol{\tau}_{\mathrm{virt}} \in \mathbb{R}^{N_{\mathrm{virt}} \times N_{\mathrm{virt}}}` directly from the AIMD trajectory energy fluctuations.

Microscopic Derivation from Short-Time Cumulant Expansion
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For any pair of adiabatic states :math:`(i, j)`, the instantaneous energy difference along the trajectory is:

.. math::

   \Delta \varepsilon_{ij}(t) = \varepsilon_i(t) - \varepsilon_j(t)

The thermal fluctuation around its equilibrium trajectory mean is:

.. math::

   \delta \Delta \varepsilon_{ij}(t) = \Delta \varepsilon_{ij}(t) - \langle \Delta \varepsilon_{ij} \rangle

with variance:

.. math::

   \sigma_{ij}^2 = \operatorname{Var}\left( \Delta \varepsilon_{ij} \right) = \langle \left( \delta \Delta \varepsilon_{ij}(t) \right)^2 \rangle

In second-order cumulant expansion of the reduced density operator, the pure-dephasing decay function :math:`D_{ij}(t)` is driven by the bath line-shape function :math:`g_{ij}(t)`:

.. math::

   D_{ij}(t) = \exp\left( -g_{ij}(t) \right), \qquad g_{ij}(t) = \frac{1}{\hbar^2} \int_0^t dt_1 \int_0^{t_1} dt_2 \, \langle \delta \Delta \varepsilon_{ij}(0) \delta \Delta \varepsilon_{ij}(t_2) \rangle

On timescales shorter than the characteristic nuclear phonon correlation time (:math:`t \ll \tau_{\mathrm{bath}} \sim 50\text{ fs}`), the energy fluctuation autocorrelation is essentially static: :math:`\langle \delta \Delta \varepsilon_{ij}(0) \delta \Delta \varepsilon_{ij}(t_2) \rangle \approx \sigma_{ij}^2`.

The double time integral simplifies analytically to:

.. math::

   g_{ij}(t) \approx \frac{\sigma_{ij}^2}{\hbar^2} \int_0^t dt_1 \, t_1 = \frac{\sigma_{ij}^2 \, t^2}{2 \hbar^2}

yielding a Gaussian dephasing profile:

.. math::

   D_{ij}(t) \approx \exp\left( -\frac{\sigma_{ij}^2 \, t^2}{2 \hbar^2} \right) = \exp\left( - \left(\frac{t}{\tau_{ij}}\right)^2 \right)

Defining the characteristic dephasing time :math:`\tau_{ij}` at the standard :math:`1/e` decay threshold (:math:`D_{ij}(\tau_{ij}) = 1/e`):

.. math::

   \frac{\sigma_{ij}^2 \, \tau_{ij}^2}{2 \hbar^2} = 1 \implies \tau_{ij} = \frac{\hbar \sqrt{2}}{\sigma_{ij}}

To prevent numerical singularities along the diagonal (:math:`i = j` where :math:`\sigma_{ii} = 0`), the diagonal elements are set to a maximum physical cutoff :math:`\tau_{\max} = 500\text{ fs}`:

.. math::

   \tau_{ij} = \min\left( \tau_{\max}, \, \frac{\hbar \sqrt{2}}{\sigma_{ij}} \right)

The resulting pairwise dephasing matrices are precomputed via ``qdex --namd-decoherence`` and cached in ``decoherence_times.npz``, where both **CPA-FSSH-EDC** and **DISH** dynamically consume them at run time.


Surface Hopping Schemes: FSSH-EDC vs. DISH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To prevent the Quantum Zeno effect (where resetting wavepackets at every nuclear step freezes carrier cooling) while curing FSSH overcoherence, ``QDEX`` provides two physically rigorous stochastic algorithms:

1. **Continuous CPA-FSSH with Energy-Based Decoherence Correction (EDC)**
   (Granucci & Persico, *J. Chem. Phys.* 126, 134114, 2007)
   Electronic amplitudes :math:`\mathbf{C}(t)` persist continuously across nuclear steps. At each nuclear step:
   * Electronic TDSE is integrated over the nuclear interval to accumulate Tully's Fewest Switches flux :math:`g_{I \to J}`.
   * If a hop occurs: the active state switches to :math:`J` and the wavepacket collapses onto state :math:`J` (:math:`c_J = 1`, :math:`c_{K \neq J} = 0`).
   * If no hop occurs: off-diagonal amplitudes are continuously damped by the state-pair dephasing factor:
     
     .. math::
     
        c_J \leftarrow c_J \exp\left( -\frac{\Delta t}{\tau_{IJ}} \right), \quad \forall J \neq I
        
     and the active amplitude :math:`c_I` is renormalized to conserve total probability :math:`\sum_K |c_K|^2 = 1`.

2. **Decoherence-Induced Surface Hopping (DISH)**
   (Jaeger, Fischer, Prezhdo, *J. Chem. Phys.* 137, 22A545, 2012; Akimov & Prezhdo, *J. Chem. Phys.* 138, 124102, 2013)
   In DISH, electronic propagation is strictly unitary between stochastic collapse events (no continuous artificial damping). At each nuclear step:
   * For each state :math:`J \neq I`, a Poisson dephasing event occurs with probability:
     
     .. math::
     
        P_{\mathrm{dec}, J} = 1 - \exp\left( -\frac{\Delta t}{\tau_{IJ}} \right)
        
   * If a dephasing event occurs for state :math:`J`, the system attempts a collapse (hop) to state :math:`J` with probability:
     
     .. math::
     
        P_{\mathrm{hop}, J} = |c_J|^2 \times \min\left( 1, \, e^{-\beta \max(E_J - E_I, 0)} \right)
        
   * If accepted, the trajectory hops to state :math:`J` and collapses into a pure state (:math:`c_J = 1`).
   * If rejected, the coherence with state :math:`J` is quenched (:math:`c_J \to 0`) and the remaining states are renormalized.

In dense semiconductor nanocrystals, DISH analytically recovers the Pauli Master Equation cooling rates while maintaining individual trajectory statistics.

---

7. Phonon Spectral Density J(ω): Mapping Electron-Phonon Coupling
-----------------------------------------------------------------

Mathematical Definition
~~~~~~~~~~~~~~~~~~~~~~~

The **Phonon Spectral Density** :math:`J(\omega)` is the Fourier transform of the energy gap autocorrelation function:

.. math::

   J(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} C(t) \, e^{i \omega t} \, dt

``QDEX`` stores the non-negative real part of a Hann-windowed Fourier transform of :math:`C(t)`, in :math:`\mathrm{cm}^{-1}`. That is the cosine transform of a real, even autocorrelation. It is not :math:`|\mathrm{FFT}|^2`. On a record of length :math:`T` the Rayleigh spacing is :math:`1/T` (about :math:`33\,\mathrm{cm}^{-1}` for a 1 ps trajectory at 2 fs), and that spacing is printed next to the peak.

Physical Meaning
~~~~~~~~~~~~~~~~

:math:`J(\omega)` provides the **frequency-resolved spectrum of nuclear vibrations that couple to the electronic transitions**. The area under a peak at frequency :math:`\omega` is proportional to the electron-phonon coupling strength (Huang-Rhys parameter :math:`S_\alpha`) for that vibrational mode:

.. math::

   J(\omega) = \pi \sum_\alpha \omega_\alpha^2 \, S_\alpha \, \delta(\omega - \omega_\alpha)

Identifying Active Phonon Modes During Cooling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By inspecting the peaks in :math:`J(\omega)`, researchers can directly identify which phonon modes facilitate carrier cooling and dephasing:

* **Low-Frequency Acoustic Modes (:math:`< 50\text{ cm}^{-1}`)**: Acoustic phonons provide a continuous low-energy bath that mediates intra-band thermalization within dense manifolds.
* **Optical Modes (e.g. Pb–Br / Pb–I Stretching, :math:`60 - 150\text{ cm}^{-1}`)**: Polar optical phonons create strong macroscopic electric fields (Fröhlich interaction), driving fast non-adiabatic transitions across intermediate energy gaps.
* **Organic Cation / Ligand Modes (:math:`200 - 300\text{ cm}^{-1}`)**: In hybrid perovskites (:math:`\text{MAPbI}_3`), rotational and librational motions of methylammonium cations produce high-frequency peaks in :math:`J(\omega)` that help bridge larger energy spacings.

---

8. Radiative & Non-Radiative Recombination Mechanisms
-----------------------------------------------------

Once photoexcited hot carriers have cooled to the band edge (forming the lowest :math:`1S` exciton state), they recombine to the electronic ground state :math:`|S_0\rangle` on longer timescales (nanoseconds). In ``QDEX``, carrier recombination is coupled directly to the population dynamics:

.. math::

   \frac{d P_I(t)}{dt} = \sum_{J \neq I} \left[ k_{J \to I} P_J - k_{I \to J} P_I \right] - \left( k_{I \to 0}^{\mathrm{rad}} + k_{I \to 0}^{\mathrm{nr}} \right) P_I

where :math:`k_{I \to 0}^{\mathrm{rad}}` and :math:`k_{I \to 0}^{\mathrm{nr}}` represent the radiative and non-radiative recombination rates from excited state :math:`|I\rangle` to the ground state.

1. Microscopic Origin: From DFT and BSE to Recombination
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To avoid empirical approximations, all transition energies and optical matrix elements originate directly from the ground-state DFT and excited-state frameworks:

1. **Ground-State DFT (CP2K)**:
   Self-consistent field calculations yield the Kohn-Sham molecular orbitals :math:`\psi_i(\mathbf{r}) = \sum_\mu C_{\mu i} \chi_\mu(\mathbf{r})` and orbital eigenvalues :math:`\varepsilon_i`. The atomic-orbital transition dipole matrices are computed analytically via Libint:

   .. math::

      \boldsymbol{\mu}_{\mu\nu}^{\mathrm{AO}} = \langle \chi_\mu | e \mathbf{r} | \chi_\nu \rangle

2. **Spin-Orbit Coupling (SOC)**:
   When SOC is active, spatial orbitals are projected into two-component spinors:

   .. math::

      \Psi_n^{\mathrm{SP}}(\mathbf{r}) = U_{n\alpha} \psi_\alpha(\mathbf{r}) + U_{n\beta} \psi_\beta(\mathbf{r})

   The transition dipole matrix is rotated into the spinor basis:

   .. math::

      \boldsymbol{\mu}_{ia}^{\mathrm{SP}} = \mathbf{U}_{\mathrm{occ},\alpha}^\dagger \mathbf{M} \mathbf{U}_{\mathrm{virt},\alpha} + \mathbf{U}_{\mathrm{occ},\beta}^\dagger \mathbf{M} \mathbf{U}_{\mathrm{virt},\beta}

   where :math:`\mathbf{M} = \mathbf{C}_{\mathrm{act}}^T \boldsymbol{\mu}^{\mathrm{AO}} \mathbf{C}_{\mathrm{act}}`.

3. **Quasiparticle Corrections & Excitation Energies**:
   Applying the GW scissor shift :math:`\Delta_{\mathrm{scissor}}` yields quasiparticle eigenvalues :math:`\tilde{\varepsilon}_i = \varepsilon_i - \Delta_{\mathrm{scissor}} f_i^{\mathrm{HOMO}}` and :math:`\tilde{\varepsilon}_a = \varepsilon_a + \Delta_{\mathrm{scissor}} f_a^{\mathrm{LUMO}}`.

   * **In Diagonal BSE (``diagonal_bse``)**:
     The transition energy of an electron-hole pair :math:`|i \to a\rangle` incorporates the direct screened Coulomb attraction:

     .. math::

        E_{ia} = \tilde{\varepsilon}_a - \tilde{\varepsilon}_i - K_d(ia, ia)

     where :math:`K_d(ia, ia) = \sum_{AB} q_{ii}^A W_{AB} q_{aa}^B` is the Ohno-Klopman or Resta-MNOK screened kernel.
   * **In Full BSE / sTDA (``bse``, ``stda``)**:
     Exciton states are coherent superpositions of single-particle transitions:

     .. math::

        |\Psi_I\rangle = \sum_{ia} c_{ia}^I |i \to a\rangle

     with total transition dipole :math:`\boldsymbol{\mu}_I = \sum_{ia} c_{ia}^I \boldsymbol{\mu}_{ia}^{\mathrm{SP}}`.

4. **Oscillator Strengths**:
   The dimensionless oscillator strength for transition :math:`|I\rangle \to |S_0\rangle` is evaluated as:

   .. math::

      f = \begin{cases}
      \dfrac{4}{3}\, E_{\mathrm{Ha}} \, |\boldsymbol{\mu}|^2 & \text{closed-shell spatial orbital} \\
      \dfrac{2}{3}\, E_{\mathrm{Ha}} \, |\boldsymbol{\mu}|^2 & \text{spinor}
      \end{cases}

   with the excitation energy in Hartree and the dipole in :math:`ea_0`. This is the atomic-unit reduction of :math:`f = (2/3)(m_e/\hbar^2) E |\mu|^2` per spin-orbital, doubled for a singlet that uses a spatial orbital. Older precomputes that stored :math:`(2/3) E_{\mathrm{eV}} |\mu|^2` are rescaled when the dynamics read them.

2. Einstein Radiative Rate: Single-Frame vs. NAMD Trajectory Averaging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spontaneous photon emission into the vacuum radiation field inside a dielectric medium of refractive index :math:`n_{\mathrm{r}}` is given by the Einstein A coefficient:

.. math::

   k_{I \to 0}^{\mathrm{rad}} = \left[ \frac{2 e^2}{4\pi \epsilon_0 m_e c^3 \hbar^2} \right] n_{\mathrm{r}} \, E_I^2 \, f_I = C_{\mathrm{rad}} \, n_{\mathrm{r}} \, E_I^2 \, f_I

where :math:`C_{\mathrm{rad}} = 4.3391988 \times 10^7\text{ s}^{-1}\text{ eV}^{-2}` (:math:`4.3391988 \times 10^{-8}\text{ fs}^{-1}\text{ eV}^{-2}`) and :math:`n_{\mathrm{r}}` is loaded from ``REFRACTIVE_INDEX_DICT`` in ``qdex.hardness`` (e.g. :math:`n_{\mathrm{r}} = 2.19` for :math:`\text{CsPbBr}_3`, :math:`n_{\mathrm{r}} = 3.30` for :math:`\text{GaAs}`).

Single-Frame vs. Trajectory Ensemble Averaging
""""""""""""""""""""""""""""""""""""""""""""""

A critical question is whether :math:`f_I` and :math:`E_I` should be taken from a single static snapshot (frame 0) or averaged along the NAMD trajectory:

1. **Static / Single-Frame Rate (:math:`k_{\mathrm{rad}}(t=0)`)**:
   Evaluates :math:`E_I(0)` and :math:`f_I(0)` at the relaxed ground-state equilibrium geometry.
2. **Thermalized Band-Edge Rate at Frame 0**:
   Because the fine-structure splitting between band-edge exciton states (:math:`1-10\text{ meV}`) is much smaller than thermal energy (:math:`k_B T \approx 25.8\text{ meV}` at :math:`300\text{ K}`), carriers rapidly reach thermal equilibrium among low-lying states before radiating:

   .. math::

      k_{\mathrm{rad}}^{\mathrm{therm}}(t=0) = \frac{\sum_I k_{\mathrm{rad}, I}(0) \, \exp\left( -\frac{E_I(0) - E_0(0)}{k_B T} \right)}{\sum_I \exp\left( -\frac{E_I(0) - E_0(0)}{k_B T} \right)}

3. **What the dynamics actually use**:
   Dipoles are computed on frame 0. The photoluminescence yield uses the frame-0 thermal average above, with the lowest-exciton lifetime taken from :math:`\arg\min E_I`, not from pair index 0. A trajectory average :math:`\langle k_{\mathrm{rad}}\rangle_{\mathrm{MD}}` would require the dipole on every frame. That average is not computed. A dark band-edge exciton at the first geometry therefore stays dark for the reported yield. Herzberg–Teller intensity borrowing along the trajectory is a real physical effect and is not yet in the rate.

3. Non-Radiative Decay Across Large Gaps: Englman-Jortner Energy Gap Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Direct non-radiative recombination across a wide semiconductor band gap (:math:`E_g > 1.0\text{ eV}`) requires dissipating enormous electronic energy into the nuclear lattice. Because this involves dozens of vibrational quanta (:math:`p = E_g / \hbar \omega_{\mathrm{LO}} \sim 30 - 80`), multi-phonon perturbation theory (Englman & Jortner, 1970) yields the **Energy Gap Law**:

.. math::

   k_{\mathrm{nr}} = A_{\mathrm{nr}} \exp\left( -\gamma \frac{E_g}{\hbar \omega_{\mathrm{LO}}} \right)

where:

* :math:`\hbar \omega_{\mathrm{LO}}` is the dominant accepting optical phonon energy.
* :math:`\gamma` is the multi-phonon coupling parameter:

  .. math::

     \gamma = \ln\left( \frac{E_g}{S \, \hbar \omega_{\mathrm{LO}}} \right) - 1 = \ln\left( \frac{E_g}{\lambda} \right) - 1

  where :math:`S` is the dimensionless Huang-Rhys factor and :math:`\lambda = S \hbar \omega_{\mathrm{LO}}` is the nuclear reorganization energy.
* :math:`A_{\mathrm{nr}}` is the electronic prefactor:

  .. math::

     A_{\mathrm{nr}} = \frac{C_{\mathrm{el}}^2}{\hbar} \sqrt{\frac{2\pi}{\hbar \omega_{\mathrm{LO}} E_g}} \sim 10^{12} - 10^{13}\text{ s}^{-1}

  where :math:`C_{\mathrm{el}} \approx V_{\mathrm{el}}` is the non-adiabatic coupling matrix element between excited and ground states.

4. Non-Empirical Extraction of Optical Phonon Energy from NAMD Spectral Density
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than relying on empirical phonon frequencies, ``QDEX`` extracts :math:`\hbar \omega_{\mathrm{LO}}` directly from the **Phonon Spectral Density** :math:`J(\omega)` of the NAMD trajectory:

1. The instantaneous energy gap fluctuation of the lowest transition along the MD trajectory is tracked:

   .. math::

      \delta E_g(t) = E_g(t) - \langle E_g \rangle

2. The normalized gap time-autocorrelation function is evaluated:

   .. math::

      C(t) = \frac{\langle \delta E_g(0) \delta E_g(t) \rangle}{\sigma_E^2}

3. A Hann-windowed Fourier transform stores :math:`\max(\mathrm{Re}\,\mathrm{FFT}[C], 0)` as :math:`J(\omega)` in wavenumbers (:math:`\text{cm}^{-1}`). The spacing of an independent frequency bin is :math:`1/T`:

   .. math::

      J(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} C(t) W(t) \, e^{i \omega t} \, dt

4. The dominant optical phonon mode is identified from the primary peak of :math:`J(\omega)` (excluding low-frequency acoustic noise :math:`< 30\text{ cm}^{-1}`):

   .. math::

      \tilde{\nu}_{\mathrm{LO}} = \operatorname{argmax}_{\tilde{\nu} \ge 30\text{ cm}^{-1}} J(\tilde{\nu})

5. Converting from wavenumber to energy gives the optical phonon quantum:

   .. math::

      \hbar \omega_{\mathrm{LO}} = h c \, \tilde{\nu}_{\mathrm{LO}} = (1.23984 \times 10^{-4}\text{ eV}\cdot\text{cm}) \times \tilde{\nu}_{\mathrm{LO}}

   *For example, in lead halide perovskites (:math:`\text{CsPbBr}_3`), the dominant peak at :math:`\tilde{\nu} \approx 150\text{ cm}^{-1}` yields :math:`\hbar \omega_{\mathrm{LO}} = 18.6\text{ meV}`. In CdSe nanocrystals (:math:`\tilde{\nu} \approx 210\text{ cm}^{-1}`), it yields :math:`\hbar \omega_{\mathrm{LO}} = 26.0\text{ meV}`.*

5. Derivation of Huang-Rhys Factor S and Reorganization Energy λ from Trajectory Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A central parameter in multi-phonon transitions is the **Huang-Rhys factor** :math:`S`, which quantifies the average number of phonons emitted during electronic transition. In ``QDEX``, :math:`S` and the nuclear reorganization energy :math:`\lambda` are determined non-empirically via the **Fluctuation-Dissipation Theorem / Marcus linear response theory**:

1. **Thermal Gap Variance**:
   Along the ab initio trajectory at temperature :math:`T`, the classical variance of the energy gap is computed directly:

   .. math::

      \sigma_E^2 = \langle (E_g(t) - \langle E_g \rangle)^2 \rangle = \frac{1}{N_{\mathrm{frames}}} \sum_{k=1}^{N_{\mathrm{frames}}} \delta E_g(t_k)^2

2. **Nuclear Reorganization Energy (:math:`\lambda`)**:
   In linear response theory for a harmonic bath in the classical limit (:math:`k_B T \gg \hbar \omega / 2`), the energy gap variance is directly proportional to the reorganization energy:

   .. math::

      \sigma_E^2 = 2 \lambda \, k_B T \implies \lambda = \frac{\sigma_E^2}{2 \, k_B T}

   *(Quantum mechanically, this corresponds to :math:`\sigma_E^2 = \int_0^\infty \frac{2}{\pi} \hbar \omega J(\omega) \coth\left(\frac{\hbar \omega}{2 k_B T}\right) d\omega`).*

3. **Huang-Rhys Factor (:math:`S`)**:
   Because the total reorganization energy partitioned into the dominant optical phonon mode of frequency :math:`\hbar \omega_{\mathrm{LO}}` is :math:`\lambda = S \, \hbar \omega_{\mathrm{LO}}`, we solve directly for :math:`S`:

   .. math::

      S = \frac{\lambda}{\hbar \omega_{\mathrm{LO}}} = \frac{\sigma_E^2}{2 \, k_B T \, \hbar \omega_{\mathrm{LO}}}

**Physical Significance**: Both :math:`\lambda` and :math:`S` are extracted directly from the NAMD trajectory without any adjustable parameters. Soft, polar perovskite lattices with large thermal gap fluctuations (:math:`\sigma_E \approx 60\text{ meV}`) yield :math:`\lambda \approx 70\text{ meV}` and :math:`S \approx 3.8`, reflecting significant electron-phonon coupling, whereas rigid covalent nanocrystals exhibit :math:`S \approx 0.5 - 1.5`.

6. Intermediate & Narrow Gap Decay: Franck-Condon Weighted Density of States (FCWD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In narrow-gap quantum dots (:math:`E_g < 1.0\text{ eV}`, e.g. infrared PbS, PbSe, HgTe, or InAs) or near crossing points with shallow defect states, the number of required accepting phonons is small (:math:`p < 15`). In this regime, the discrete multi-phonon expansion transitions into a continuous **Franck-Condon Weighted Density of States (FCWD)** evaluated from Fermi's Golden Rule:

.. math::

   k_{\mathrm{nr}}^{\mathrm{FCWD}} = \frac{2\pi}{\hbar} |V_{\mathrm{el}}|^2 \, \mathrm{FCWD}(E_g)

In the Marcus-Levich framework, the nuclear Franck-Condon factor takes a Gaussian line-shape:

.. math::

   \mathrm{FCWD}(E_g) = \frac{1}{\sqrt{2\pi \sigma^2}} \exp\left( -\frac{(E_g - \lambda)^2}{2\sigma^2} \right)

Data-Driven Extraction of FCWD Parameters
""""""""""""""""""""""""""""""""""""""""""

Every variable in this rate expression is evaluated directly from the NAMD trajectory:

1. **Gaussian Broadening (:math:`\sigma`)**:
   :math:`\sigma` is **not an arbitrary broadening parameter**! It is the exact standard deviation of the energy gap fluctuations from the MD trajectory:

   .. math::

      \sigma = \sigma_E = \sqrt{\langle (E_g(t) - \langle E_g \rangle)^2 \rangle} = \sqrt{2 \lambda k_B T}

2. **Nuclear Reorganization Energy (:math:`\lambda`)**:
   Computed from the trajectory variance: :math:`\lambda = \frac{\sigma_E^2}{2 k_B T}`.
3. **Electronic prefactor**:
   The wide-gap loss uses the Englman–Jortner law with the trajectory Huang–Rhys factor :math:`S`, the phonon energy taken from the peak of :math:`J(\omega)`, and a configured prefactor :math:`A_{\mathrm{nr}}` (default :math:`10^{13}\,\mathrm{s}^{-1}`). The root-mean-square intraband coupling :math:`\hbar\langle|d|\rangle` connects excited states to each other. It is not the exciton-to-ground matrix element, and it is not inserted as :math:`V_{\mathrm{el}}`. The Marcus–Levich Gaussian is appropriate only when the gap is a few :math:`\sigma`. It is not the loss used for a :math:`\mathrm{CsPbBr}_3` gap of about 3 eV.

7. Defect Trap-Assisted Recombination (Shockley-Read-Hall)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In real quantum dots with unpassivated surfaces or vacancies, non-radiative recombination is overwhelmingly accelerated by **Shockley-Read-Hall (SRH) mid-gap traps**. Instead of bridging a single large gap of :math:`2.0\text{ eV}`, carriers drop into an intermediate trap state (:math:`\Delta E \approx 0.5 - 1.0\text{ eV}`), where multi-phonon tunneling is orders of magnitude faster.

In ``QDEX``, trap-assisted recombination can be configured via the YAML input:

.. code-block:: yaml

   namd:
     recombination:
       include_ground_state: true
       tau_nr_ns: 25.0  # Trap non-radiative lifetime in nanoseconds

8. Photoluminescence Quantum Yield (PLQY)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total Photoluminescence Quantum Yield (PLQY) represents the branching ratio between radiative photon emission and non-radiative dissipation:

.. math::

   \mathrm{PLQY} = \frac{\langle k_{\mathrm{rad}} \rangle}{\langle k_{\mathrm{rad}} \rangle + k_{\mathrm{nr}}} \times 100\%

* In pristine, defect-free nanocrystals where :math:`E_g \gg \hbar \omega_{\mathrm{LO}}`, the intrinsic non-radiative rate is negligible (:math:`k_{\mathrm{nr}} \ll 10^3\text{ s}^{-1}`), leading to near-unity PLQY (:math:`\sim 99\%`).
* In the presence of surface traps (:math:`\tau_{\mathrm{nr}} \sim 10 - 50\text{ ns}`), non-radiative decay competes directly with radiative emission (:math:`\tau_{\mathrm{rad}} \sim 5 - 20\text{ ns}`), yielding realistic PLQYs between :math:`20\%` and :math:`70\%`.

---

9. In-Depth Analysis of NAMD Simulations
----------------------------------------

``QDEX`` includes a dedicated analysis module (``qdex.namd.analysis``) that automatically processes precomputed and dynamic trajectory data.

1. Carrier Cooling Curves, Lifetimes, and Band Edge Arrival Times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Carrier relaxation is monitored by tracking the time-dependent excess energy of electrons (:math:`\Delta E_e`), holes (:math:`\Delta E_h`), and total exciton (:math:`\Delta E_{\mathrm{exc}}`) above their respective band edges:

.. math::

   \Delta E_e(t) = \sum_{a \in \mathrm{virt}} P_a(t) \left( \varepsilon_a(t) - \varepsilon_{\mathrm{LUMO}}(t) \right)

.. math::

   \Delta E_h(t) = \sum_{i \in \mathrm{occ}} P_i(t) \left( \varepsilon_{\mathrm{HOMO}}(t) - \varepsilon_i(t) \right)

Cooling Rates and Exponential Lifetimes
"""""""""""""""""""""""""""""""""""""""

The cooling lifetimes :math:`\tau_{\mathrm{exc}}`, :math:`\tau_e`, and :math:`\tau_h` are extracted by linear regression on the logarithmic excess energy decay:

.. math::

   \Delta E(t) \approx \Delta E(0) \, \exp\left( -\frac{t}{\tau_{\mathrm{cooling}}} \right)

The corresponding carrier cooling rate is:

.. math::

   k_{\mathrm{cool}} = \frac{1}{\tau_{\mathrm{cooling}}} \quad [\text{ps}^{-1}]

Time to Reach the Band Edge: Analytical Estimates vs. Actual Trajectory
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In experiments and device modeling, researchers require not only the initial relaxation slope :math:`\tau`, but also the **total time required for hot carriers to fully reach the band edge**:

1. **Analytical Estimates from Exponential Fit**:
   * **95% Excess Energy Dissipated**:

     .. math::

        t_{\mathrm{est}}^{95\%} = -\ln(0.05) \, \tau \approx 3.0 \, \tau

   * **99% Excess Energy Dissipated (Complete Thermalization)**:

     .. math::

        t_{\mathrm{est}}^{99\%} = -\ln(0.01) \, \tau \approx 4.6 \, \tau

2. **Actual Numerical Arrival Times from Simulation Trajectory**:
   * :math:`t_{\mathrm{act}}^{95\%}`: The first simulation timestamp :math:`t` where :math:`\Delta E(t) \le 0.05 \, \Delta E(0)`.
   * :math:`t_{\mathrm{act}}^{99\%}`: The first simulation timestamp :math:`t` where :math:`\Delta E(t) \le 0.01 \, \Delta E(0)`.
   * :math:`t_{\mathrm{act}}^{\mathrm{therm}}`: The first timestamp where excess energy drops below the thermal energy of the lattice bath:

     .. math::

        \Delta E(t) \le k_B T \approx 25.8\text{ meV at } 300\text{ K}

**Physical Insight**: While an ideal exponential decay satisfies :math:`t_{\mathrm{act}} \approx t_{\mathrm{est}}`, realistic atomistic trajectories often exhibit non-exponential behaviors—such as an initial **phonon bottleneck** across discrete sub-bands or delayed cascades through intermediate surface states. Comparing :math:`t_{\mathrm{est}}` with :math:`t_{\mathrm{act}}` immediately diagnostics whether carrier cooling proceeds smoothly or is delayed by bottlenecks.

2. State-Resolved Population Kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Transient populations :math:`P_I(t)` are exported to ``carrier_cooling_populations.csv``, displaying the sequential decay of initial hot excitons into intermediate states and finally into the emitting :math:`1S` state.

3. NAC vs. Energy Gap Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To verify whether non-adiabatic transitions obey the energy-gap law, ``QDEX`` samples pairs of states across trajectory frames and plots non-adiabatic coupling magnitudes :math:`|d_{IJ}|` against energy differences :math:`|E_J - E_I|`. This distinguishes smooth exponential decay from resonant vibronic enhancements.

4. 6-Panel Publication Figures & Dashboards
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Executing the analysis workflow generates a comprehensive 6-panel summary figure (``namd_analysis_6panel.png``) and an interactive Plotly HTML dashboard:
* **Panel A**: Carrier cooling curves (:math:`\Delta E_e(t)` vs. :math:`\Delta E_h(t)` with fitted lifetimes).
* **Panel B**: Time-dependent populations of frontier exciton states.
* **Panel C**: Band-gap thermal fluctuation trajectory :math:`E_g(t)`.
* **Panel D**: Energy gap autocorrelation function :math:`C(t)` and cumulant dephasing decay :math:`D(t)`.
* **Panel E**: Phonon Spectral Density :math:`J(\omega)` in :math:`\text{cm}^{-1}`.
* **Panel F**: Non-adiabatic coupling distribution :math:`|d_{IJ}|` vs. :math:`\Delta E_{IJ}`.

---

10. Ultrafast Pump-Probe Transient Absorption (TA) Spectroscopy
---------------------------------------------------------------

In ultrafast optical experiments on semiconductor nanocrystals, **pump-probe transient absorption (TA)** spectroscopy tracks the non-equilibrium evolution of photoexcited carriers. A high-energy pump pulse creates an initial non-thermal carrier distribution, while a broadband, time-delayed probe pulse monitors the differential absorbance:

.. math::

   \Delta A(E, t) = A_{\mathrm{pump-on}}(E, t) - A_{\mathrm{pump-off}}(E)

as a function of probe photon energy :math:`E` and pump–probe delay time :math:`t`.

Physical Mechanisms in Nanocrystal Transient Absorption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transient absorption signal :math:`\Delta A(E, t)` comprises three distinct physical mechanisms:

1. **Ground-State Bleach (GSB) & State-Filling (Pauli Blocking)** (:math:`\Delta A < 0`):
   When electrons occupy virtual conduction states :math:`a` and holes occupy valence states :math:`i`, optical transitions between these occupied states are blocked by the Pauli exclusion principle. As carriers cascade down the non-adiabatic ladder toward the band edge, absorption into the lowest exciton state (:math:`1S`) is progressively quenched, forming the prominent negative **1S Bleach**.

2. **Stimulated Emission (SE)** (:math:`\Delta A < 0`):
   Probe photons passing through the sample stimulate coherent radiative transitions from populated excited levels back to the ground state. Because stimulated emission adds photons to the transmitted probe beam, it appears as a negative differential absorption feature at the same transition frequencies as the ground-state bleach.

3. **Excited-State Absorption (ESA)** (:math:`\Delta A > 0`):
   A populated excited carrier can absorb a second probe photon to transition into higher-lying continuum bands or multiexciton states (:math:`S_1 \to S_{XX}`), producing positive absorption features.

Microscopic Formulation & State-Filling Factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In ``QDEX``, the instantaneous exciton pair populations :math:`P_{ia}(t)` (where :math:`i` denotes occupied valence orbitals and :math:`a` denotes virtual conduction orbitals) directly govern state filling.

At delay time :math:`t_k`, the single-particle fractional electron occupation in conduction orbital :math:`a` and hole occupation in valence orbital :math:`i` are obtained by tracing out the complementary carrier:

.. math::

   n_a(t_k) = \sum_{i} P_{ia}(t_k), \qquad p_i(t_k) = \sum_{a} P_{ia}(t_k)

In the unexcited ground state, each transition :math:`(i \to a)` has an intrinsic oscillator strength :math:`f_{ia}^{(0)}`. For incoherent populations the differential strength is the degeneracy-normalized state-filling signal used for quantum-dot transient absorption (Grimaldi et al., Nano Lett. 2019). A closed-shell spatial orbital holds two carriers, so :math:`g = 2`; a spinor has :math:`g = 1`. One carrier then blocks a fraction :math:`1/g` of that shell. Net optical gain requires :math:`n_a/g_a + p_i/g_i > 1`.

.. math::

   \Delta f_{ia}(t_k) = - f_{ia}^{(0)} \left( \frac{n_a(t_k)}{g_a} + \frac{p_i(t_k)}{g_i} \right)

An extra :math:`-f P_{ia}` term is available only through ``include_se`` and is not the default. Adding it on top of state filling counts stimulated emission twice.

The continuous differential absorption spectrum :math:`\Delta A(E, t_k)` is then evaluated by convolving with a Gaussian probe spectral broadening :math:`\sigma`:

.. math::

   \Delta A(E, t_k) = \sum_{ia} \Delta f_{ia}(t_k) \frac{1}{\sqrt{2\pi}\sigma} \exp\left( -\frac{(E - E_{ia}(t_k))^2}{2\sigma^2} \right)

1S Bleach Kinetic Profiling & Carrier Cooling Rates (:math:`k_C`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The hallmark signature of carrier cooling in transient absorption experiments is the **rise profile of the 1S bleach**. 

At :math:`t = 0`, photoexcitation creates carriers at high excess energy (e.g. :math:`2 \times E_g`). The initial bleach :math:`\Delta A(E, 0)` is broad and located in the high-energy window. As non-adiabatic electron-phonon scattering drives carriers toward the band edge, population accumulates in the lowest exciton state (:math:`1S`, :math:`E_{1S} \approx E_g`), causing the negative 1S bleach signal :math:`-\Delta A_{1S}(t)` to rise from zero to its plateau value.

``QDEX`` extracts the band-edge 1S bleach:

.. math::

   S_{1S}(t) = - \Delta A(E_{1S}, t)

and performs non-linear least-squares fitting to an exponential rise model:

.. math::

   S_{1S}(t) = A_0 \left( 1 - \exp\left( -\frac{t}{\tau_C} \right) \right)

This directly yields the **hot-carrier cooling time** :math:`\tau_C` (in fs/ps) and **cooling rate** :math:`k_C = 1 / \tau_C` (in :math:`\text{ps}^{-1}`), providing direct theoretical counterparts to the experimental 1S bleach rise traces reported in transient absorption spectroscopy.

Schematic Diagrams of Transient Absorption Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optical mechanisms (GSB, SE, ESA) and dynamical carrier relaxation cascade underlying transient absorption are illustrated schematically below:

**1. Three Optical Mechanisms in Nanocrystal Transient Absorption**

.. code-block:: text

   =========================================================================================================
             THREE OPTICAL MECHANISMS IN TRANSIENT ABSORPTION SPECTROSCOPY (ΔA)
   =========================================================================================================

     (a) Ground-State Bleach (GSB)            (b) Stimulated Emission (SE)             (c) Excited-State Absorption (ESA)
         ΔA < 0  (Pauli Blocking)                 ΔA < 0  (Probe Photons Added)            ΔA > 0  (Induced Absorption)

         Energy                                   Energy                                   Energy
           ▲                                        ▲                                        ▲
           │                                        │                                        │  Conduction Continuum
           │                                        │                                        │  ┌─────┐ a'
           │                                        │                                        │  │     │
           │                                        │                                        │  └─────┘
           │                                        │                                        │     ▲
           │                                        │                                        │     │  Probe Absorbed:
           │                                        │                                        │     │  hν_probe = E_a' - E_a
    CB  ┼──┼── CBM                               CB ┼──┼── CBM                            CB ┼──┼──┤  Δf_aa'^ESA > 0
           │  ┌─────┐ a                             │  ┌─────┐ a                             │  ┌──┴──┐ a (Populated)
           │  │  ●  │ n_a (Blocked!)                │  │  ●  │ P_ia (Populated pair)         │  │  ●  │ n_a
           │  └─────┘                               │  └─────┘                               │  └─────┘
           │     ▲                                  │     │
           │     │  Probe Blocked!                  │     │  Probe Stimulates:               │
           │     │  hν_probe ≈ E_ia                 │     │  hν_probe + hν_em (2 photons!)   │
           │     │  Δf_ia^GSB < 0                   │     ▼  Δf_ia^SE < 0                    │
    VB  ┼──┼── VBM                               VB ┼──┼── VBM                            VB ┼──┼── VBM
           │  ┌─────┐ i                             │  ┌─────┐ i                             │  ┌─────┐ i
           │  │  ○  │ p_i (Blocked!)                │  │  ○  │                               │  │  ●  │
           │  └─────┘                               │  └─────┘                               │  └─────┘
           └────────────────────────►               └────────────────────────►               └────────────────────────►

      Formula Mapping:                         Formula Mapping:                         Formula Mapping:
        • n_a = Σ_i P_ia                         • Default signal is state filling        • Optional, if orbital dipoles exist
        • p_i = Σ_a P_ia                         • already. include_se adds -f P_ia      • Δf_aa' = +f_aa' n_a / g
        • Δf = -f (n_a/g + p_i/g)                • and double-counts emission.           • ΔA_ESA > 0
        • g = 2 spatial, g = 1 spinor            • Net gain only if n/g + p/g > 1
        • One spatial exciton: Δf = -f

**2. Hot-Carrier Relaxation Cascade & 1S Bleach Kinetic Rise Profile**

.. code-block:: text

   =========================================================================================================
         HOT-CARRIER RELAXATION CASCADE & 1S BLEACH KINETIC RISE PROFILE
   =========================================================================================================

     1. Hot Excitation & Vibronic Cooling Cascade           2. Time-Resolved 1S Bleach Profile S_1S(t)
        Energy                                                 Bleach: S_1S(t) = -ΔA(E_1S, t)
          ▲                                                      ▲
          │  Pump: hν_pump ≈ 2 * E_g                             │
          │    ┌─────┐ e* (Hot Electron)                         │                      Plateau A_0
          │    │  ●  │                                           │               . - - - - - - - - - - - - -
          │    └─────┘                                           │           . '
          │       │                                              │        . '
          │       │ Non-adiabatic cooling                        │      . '   S_1S(t) = A_0 (1 - e^{-t/τ_C})
          │       │ d_ab = ⟨a|∂/∂t|b⟩                            │     /
          │       ▼ (Phonon emission cascade)                    │    /
     CB ──┼─── ┌─────┐ 1S_e (CBM Band Edge)                      │   /
          │    │  ●  │ State filling n_1Se(t) rises!             │  /
          │    └─────┘                                           │ /
          │                                                      │/
          │       E_1S ≈ E_g (Probe monitors 1S)                 └────────────────────────────────────────►
          │                                                      0          τ_C               t (Delay)
     VB ──┼─── ┌─────┐ 1S_h (VBM Band Edge)
          │    │  ○  │ State filling p_1Sh(t) rises!          Kinetic Parameters:
          │    └─────┘                                          • Hot carrier cooling time: τ_C (fs or ps)
          │       ▲ (Phonon emission cascade)                   • Carrier cooling rate:     k_C = 1 / τ_C (ps⁻¹)
          │       │ Non-adiabatic cooling                       • Band-edge 1S energy:      E_1S = ε_1Se - ε_1Sh - E_b
          │       │ d_ij = ⟨i|∂/∂t|j⟩                           • Bleach signal:            -ΔA(E_1S, t)
          │       │
          │    ┌─────┐
          │    │  ○  │ h* (Hot Hole)
          │    └─────┘

How QDEX Data Are Used to Compute Every Formula Term
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is the exact step-by-step mapping of how QDEX data structures evaluate every term in the transient absorption equations:

.. list-table::
   :widths: 25 35 40
   :header-rows: 1

   * - Mathematical Term
     - QDEX Data Structure / Source
     - Computational Implementation
   * - **Transition Dipoles :math:`\boldsymbol{\mu}_{ia}`**
     - ``compute_dipole_ao(shells)`` & MO transformation
     - Dipole integrals :math:`\langle \mu | e\mathbf{r} | \nu \rangle` via Libint2, transformed to MO basis: :math:`\boldsymbol{\mu}_{ia} = C_{\mathrm{occ}}^T \boldsymbol{\mu}_{\mathrm{ao}} C_{\mathrm{virt}}`.
   * - **Ground-State Strengths :math:`f_{ia}^{(0)}`**
     - ``f_pairs`` in ``frame_00000.npz``
     - Unperturbed oscillator strengths: :math:`f_{ia}^{(0)} = \frac{2}{3} \Delta E_{ia} |\boldsymbol{\mu}_{ia}|^2`.
   * - **Instantaneous Energies :math:`E_{ia}(t)`**
     - ``E_curr`` from ``step_*.npz`` or ``all_energies``
     - Time-dependent diagonal BSE / QP transition energies tracking nuclear MD motion.
   * - **Time-Dependent Populations :math:`P_{ia}(t)`**
     - ``populations[k]`` from ``run_namd_dynamics``
     - Dynamic state populations propagated via Pauli Master Equation or CPA-FSSH.
   * - **Virtual Occupation :math:`n_a(t)`**
     - ``np.bincount(a_pairs, weights=P_k)``
     - Traced electron occupation across all active virtual orbitals :math:`a`.
   * - **Occupied Occupation :math:`p_i(t)`**
     - ``np.bincount(i_pairs, weights=P_k)``
     - Traced hole occupation across all active valence orbitals :math:`i`.
   * - **Differential Strength :math:`\Delta f_{ia}(t)`**
     - ``delta_f = - f_base * (n_a[a_arr] + p_i[i_arr] + P_k)``
     - Evaluated at every time step :math:`t_k` with vector broadcasting.
   * - **2D Differential Absorbance :math:`\Delta A(E, t)`**
     - ``np.dot(delta_f, gauss_profiles)``
     - Vectorized Gaussian spectral convolution over probe energy grid :math:`E`.
   * - **1S Rise Fitting :math:`\tau_C, k_C`**
     - ``fit_bleach_rise_kinetics(times, delta_A_1s)``
     - Scipy curve-fit of :math:`-\Delta A_{1S}(t)` to exponential rise :math:`A_0 (1 - e^{-t/\tau_C})`.

Publication-Quality Visualizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Executing transient absorption generates a 3-panel publication dashboard:
* **Panel (a): 2D False-Color TA Map :math:`\Delta A(E, t)`**: Probe energy on x-axis, delay time on y-axis, using a diverging `RdBu_r` colormap (blue = negative bleach, red = positive ESA).
* **Panel (b): 1S Bleach Kinetic Rise Profile**: Tracking the negative 1S bleach with exponential rise fit and annotated cooling time :math:`\tau_C` and rate :math:`k_C`.
* **Panel (c): Differential Absorption Spectra :math:`\Delta A(E)`**: Spectral slices at selected delay times (e.g. :math:`t = 0, 100, 250, 500, 1000\text{ fs}`) showing the spectral shift from hot state filling to the sharp band-edge bleach.

---

11. CLI Flags & YAML Configuration Reference
--------------------------------------------

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--namd-precompute``
     - ``False``
     - Execute Stage 1 trajectory precomputation (cross-overlaps, NACs, phase tracking, caching).
   * - ``--namd-run``
     - ``False``
     - Execute Stage 2 NAMD carrier cooling simulation from precomputed data.
   * - ``--namd-decoherence [dir]``
     - ``None``
     - Compute and cache state-pair pure-dephasing matrices (:math:`\tau_{ij}`) into ``decoherence_times.npz``.
   * - ``--namd-compact [dir]``
     - ``None``
     - Compress precomputed directory, eliminating redundant duplicate arrays.
   * - ``--namd-soc``
     - ``False``
     - Enable relativistic Spin-Orbit Coupling across NAMD precomputation.
   * - ``--namd-ta``
     - ``False``
     - Compute ultrafast pump-probe transient absorption (TA) spectra from NAMD dynamics.
   * - ``--namd-ta-sigma <float>``
     - ``0.03``
     - Gaussian line broadening in eV for transient absorption probe spectra.
   * - ``--namd-ta-plot``
     - ``False``
     - Generate 2D false-color TA map and 1S bleach rise kinetics plot.
   * - ``--namd-ecsh-auger``
     - ``False``
     - Enable Energy-Conserving Surface Hopping (ECSH) for two-body Auger processes during NAMD.
   * - ``--namd-biexciton``
     - ``False``
     - Initialize NAMD from a biexciton state (XX) to simulate Auger annihilation dynamics.
   * - ``--namd-initial-conditions <mode>``
     - ``"single"``
     - Set initial condition sampling: ``"single"`` (start from :math:`t_0 = 0`) or ``"multiple"`` (automated ensemble sampling across uncorrelated trajectory origins).
   * - ``--namd-multi-init``
     - ``False``
     - Convenience shortcut for ``--namd-initial-conditions multiple``.
   * - ``--namd-origins <int>``
     - ``None`` (auto)
     - Explicit number of ensemble origins (when unset, calibrated automatically from :math:`\Delta t_0 \ge 2\tau_{\mathrm{corr}}`).
   * - ``--namd-window-fs <float>``
     - ``None`` (auto)
     - Simulation window duration in fs per origin (when unset, calibrated automatically from pilot cooling :math:`3\tau_{\mathrm{cool}}`).

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   namd:
     trajectory:
       dir: "./trajectory"               # Directory containing frame_* subdirectories
       dt_nuc_fs: 2.0                    # Nuclear MD time step in femtoseconds
       start_frame: 1
       end_frame: 500

     dynamics:
       method: "dish"                    # "master_equation" (PME), "dish" (DISH), or "cpa_fssh" (FSSH-EDC)
       initial_conditions: "multiple"    # "single" (t0 = 0) or "multiple" (auto-calibrated ensemble)
       temperature_k: 300.0              # Lattice temperature for detailed balance
       tau_dec_fs: "cumulant"            # "cumulant" (ab initio), "edc", or fixed float in fs
       decoherence: "edc"                # Decoherence scheme for FSSH (continuous EDC)
       n_trajectories: 1000              # Trajectory count (split evenly across origins in multi-mode)
       detailed_balance: true            # Enforce Boltzmann detailed balance factor

     integration:
       integrator: "strang"              # Unitary Strang operator splitting
       n_substeps: 2                     # Electronic sub-steps per nuclear interval

     transient_absorption:
       run: true                         # Enable pump-probe transient absorption calculation
       sigma: 0.03                       # Probe spectral broadening in eV
       plot: true                        # Generate 2D TA map and kinetics figure
       plot_file: "transient_absorption_map.png"
       csv_file: "ta_bleach_kinetics.csv"

