Part 6: Carrier Cooling Dynamics & Photoluminescence (NAMD)
===========================================================

Following the photoexcitation of a semiconductor nanocrystal or quantum dot by an ultrashort laser pulse, high-energy ("hot") electrons and holes rapidly dissipate their excess energy through electron-phonon scattering and non-adiabatic transitions. Carriers cascade down the dense ladder of excited states, cooling toward the band edges before recombining to the ground state.

``QDEX`` features an advanced, high-throughput **Non-Adiabatic Molecular Dynamics (NAMD)** engine designed to simulate carrier relaxation, phonon bottleneck phenomena, surface defect trapping/de-trapping, and photoluminescence recombination along *ab initio* molecular dynamics (AIMD) trajectories.

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
     ├─ Preserve Diabatic State Character via Hungarian Crossing Tracking
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

Unitary Matrix Exponentiation
"""""""""""""""""""""""""""""

To propagate the electronic amplitudes :math:`\mathbf{c}(\tau)` across each sub-step without any norm drift, ``QDEX`` diagonalizes :math:`\mathbf{H}_{\mathrm{eff}} = \mathbf{V} \boldsymbol{\Lambda} \mathbf{V}^\dagger` and evaluates the exact unitary matrix exponential:

.. math::

   \mathbf{c}(\tau + \delta t_{\mathrm{elec}}) = \mathbf{V} \, \exp\left( -i \boldsymbol{\Lambda} \frac{\delta t_{\mathrm{elec}}}{\hbar} \right) \mathbf{V}^\dagger \, \mathbf{c}(\tau)

This guarantees that total electronic probability is conserved to machine precision:

.. math::

   \sum_{I} |c_I(\tau)|^2 = 1.000000000000000

Tully Hopping Flux Accumulation
"""""""""""""""""""""""""""""""

Tully's fewest switches hopping probabilities are accumulated incrementally across the electronic sub-steps:

.. math::

   g_{I \to J} = \sum_{m=1}^{N_{\mathrm{sub}}} \max\left( 0, \, \frac{-2 \delta t_{\mathrm{elec}} \, \operatorname{Re}\left( c_I^*(\tau_m) c_J(\tau_m) d_{IJ} \right)}{|c_I(\tau_m)|^2} \right) \times B_{IJ}(T)

A stochastic hopping decision is then made using the net accumulated probability over the nuclear step :math:`\Delta t_{\mathrm{nuc}}`.

Electronic Sub-Stepping in the Pauli Master Equation (PME)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although the Pauli Master Equation propagates real-valued populations :math:`P_I(t)` rather than oscillating complex amplitudes :math:`c_I(t)`, sub-stepping is equally vital for numerical stability.

In dense manifolds where non-adiabatic couplings are large, single-step forward Euler integration of :math:`\frac{d\mathbf{P}}{dt} = \mathbf{R} \mathbf{P}` with a coarse step :math:`\Delta t_{\mathrm{nuc}} \sim 1\text{ fs}` can lead to **stiffness instabilities**, causing populations to oscillate or become negative (:math:`P_I < 0`).

In ``qdex.namd.master_equation``, two robust solutions are provided:

1. **Exact Matrix Exponential**:
   For modest state spaces, the population vector is propagated analytically over the nuclear step:

   .. math::

      \mathbf{P}(t + \Delta t_{\mathrm{nuc}}) = \exp\left( \mathbf{R} \, \Delta t_{\mathrm{nuc}} \right) \mathbf{P}(t)

2. **Tensorized Sub-Stepping (Diagonal BSE)**:
   For huge manifolds (:math:`10^6` exciton pairs), the tensorized rate equations are sub-stepped with :math:`N_{\mathrm{sub}} = 20 - 50` steps (:math:`\delta t = \Delta t_{\mathrm{nuc}} / N_{\mathrm{sub}} \approx 0.02 - 0.05\text{ fs}`). Within each sub-step, loss vectors :math:`\mathbf{L}_e` and :math:`\mathbf{L}_h` continuously remove population while electron and hole transition matrices inject population, preserving strict positivity (:math:`P_{ia} \ge 0`) and probability normalization (:math:`\sum_{ia} P_{ia} = 1.0`).

---

3. Pauli Master Equation (PME) vs. CPA-FSSH: When to Use Which?
---------------------------------------------------------------

A central methodological decision in non-adiabatic dynamics is choosing between a **deterministic Master Equation** and **stochastic Fewest Switches Surface Hopping (FSSH)**. Both frameworks are implemented in ``QDEX``, and each possesses distinct physical domains of applicability:

.. list-table::
   :widths: 22 38 40
   :header-rows: 1

   * - Criterion / Regime
     - Pauli Master Equation (PME)
     - Fewest Switches Surface Hopping (FSSH)
   * - **Density of States (DOS)**
     - Dense, quasi-continuous manifolds (e.g. :math:`> 10^3 - 10^6` exciton states in large QDs).
     - Sparse or discrete level manifolds (e.g. frontier :math:`1S_e, 1P_e` states, small molecules).
   * - **Quantum Coherence**
     - Fast dephasing regime: :math:`\tau_{\mathrm{dec}} \ll \tau_{\mathrm{transfer}}`. Coherences decay before population builds up.
     - Coherent regime: quantum interference, state superpositions, and phase memory persist.
   * - **Phonon Bottleneck**
     - May overestimate relaxation if multi-phonon wavepacket dynamics are approximated by simple broad rates.
     - **Essential**: captures discrete quantum transitions and coherent vibrational wavepacket motion.
   * - **Surface Traps**
     - Provides average, memoryless Markovian trapping rates; cannot capture stochastic residence times.
     - **Essential**: tracks explicit hopping and de-hopping fluctuations, barrier crossings, and bifurcation.
   * - **Computational Cost**
     - **Ultra-fast**: Tensor decomposition propagates :math:`10^6` states in :math:`< 0.2\text{ s}` per step.
     - **Heavy**: Requires averaging over :math:`10^3 - 10^4` stochastic trajectories per initial condition.

The Phonon Bottleneck Case (:math:`1P_e \to 1S_e`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In quantum-confined semiconductor nanocrystals, quantum confinement shifts atomic orbitals into discrete, shell-like states labeled by atomic-like symmetries (:math:`1S_e, 1P_e, 1D_e`). While higher-lying states form a dense, continuous manifold where cooling is ultrafast, the energy separation between the lowest unoccupied conduction state (:math:`1S_e`) and the next state (:math:`1P_e`) can be several tenths of an electronvolt:

.. math::

   \Delta E(1P_e - 1S_e) \gg \hbar \omega_{\mathrm{LO}}

Because this energy gap greatly exceeds the energy of a single longitudinal optical (LO) phonon (:math:`\hbar \omega_{\mathrm{LO}} \approx 15 - 35\text{ meV}`), relaxation cannot occur through single-phonon scattering. This phenomenon is known as the **phonon bottleneck**.

* **Why FSSH is Recommended for Bottlenecks**: The :math:`1P_e \to 1S_e` transition is mediated by rare multi-phonon wavepacket coincidences, non-adiabatic surface crossings, or Auger-type electron-hole energy exchange. CPA-FSSH explicitly evolves the time-dependent Schrödinger equation, capturing coherent quantum interference between the discrete electronic states and vibrational wavepackets, and resolving whether the bottleneck persists or is bypassed.

Surface Trap States: Hopping & De-Hopping Kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Colloidal quantum dots frequently contain under-coordinated surface atoms, halide vacancies, or localized ligand termination defects. These defects introduce discrete electronic levels located inside the fundamental band gap (surface traps).

When a hot carrier cools to the band edge:
1. **Carrier Trapping (Hopping into Trap)**: The carrier transitions from a delocalized core state (:math:`1S`) into a spatially localized defect level. This is accompanied by strong local lattice distortion (large polaron or Jahn-Teller rearrangement).
2. **Carrier De-Trapping (Hopping out of Trap)**: Thermal fluctuations from the nuclear bath can impart sufficient energy to kick the carrier back from the defect into the delocalized band states (thermally activated de-trapping).

* **Why FSSH is Recommended for Traps**: Trapping and de-trapping are stochastic, trajectory-dependent barrier-crossing events. A deterministic rate equation (PME) treats trapping as an irreversible, memoryless Markovian decay that washes out individual carrier dwell times and trapping/detrapping equilibrium fluctuations. FSSH tracks individual stochastic trajectories: some trajectories get trapped permanently, while others hop into the trap, reside there for several picoseconds, and subsequently de-hop back into the band. Capturing this physics accurately requires both **FSSH** and **extended AIMD trajectories** (typically :math:`> 10 - 50\text{ ps}`).

Summary Decision Rule
~~~~~~~~~~~~~~~~~~~~~

* Choose **``engine: "master_equation"``** when screening carrier cooling lifetimes across dense manifolds in medium-to-large quantum dots (:math:`> 500` atoms) where the density of states is high and fast dephasing dominates.
* Choose **``engine: "surface_hopping"``** when investigating discrete frontier level transitions (:math:`1P \to 1S`), quantum coherence, or stochastic hopping/de-hopping between band edges and localized surface defect traps.

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
     1 & \text{for downward transitions } (E_J \ge E_I) \\
     \exp\left( -\frac{E_I - E_J}{k_B T} \right) & \text{for upward thermal activation } (E_J < E_I)
   \end{cases}

Vectorized Tensor Decomposition for Diagonal BSE
""""""""""""""""""""""""""""""""""""""""""""""""

In a two-particle excitonic manifold with :math:`N_{\mathrm{occ}}` occupied orbitals and :math:`N_{\mathrm{virt}}` virtual orbitals, the total number of electron-hole pairs is :math:`N_{\mathrm{pairs}} = N_{\mathrm{occ}} \times N_{\mathrm{virt}}`. Constructing and multiplying an :math:`(N_{\mathrm{pairs}} \times N_{\mathrm{pairs}})` rate matrix scales as :math:`O(N_{\mathrm{pairs}}^2) = O(N_{\mathrm{occ}}^2 N_{\mathrm{virt}}^2)`. For :math:`N_{\mathrm{occ}} = N_{\mathrm{virt}} = 500`, this corresponds to an intractable :math:`250,000 \times 250,000` dense matrix (:math:`500\text{ GB}` of RAM).

Under the **Diagonal BSE** representation, the exciton state :math:`|ia\rangle` factorizes into an independent occupied hole state :math:`i` and an independent virtual electron state :math:`a`. An exciton relaxes either via an electron transition (:math:`a \to b`) with rate :math:`K_e(a \to b)` or a hole transition (:math:`i \to j`) with rate :math:`K_h(i \to j)`:

.. math::

   \frac{d P_{ia}(t)}{dt} = \sum_{b \neq a} \left[ K_e(b \to a) P_{ib} - K_e(a \to b) P_{ia} \right] + \sum_{j \neq i} \left[ K_h(j \to i) P_{ja} - K_h(i \to j) P_{ia} \right]

In matrix notation, this decomposes into an exact **BLAS Level-3 tensor product**:

.. math::

   \frac{\partial \mathbf{P}}{\partial t} = \left( \mathbf{P} \, \mathbf{K}_e^T - \mathbf{P} \operatorname{diag}(\mathbf{L}_e) \right) + \left( \mathbf{K}_h \, \mathbf{P} - \operatorname{diag}(\mathbf{L}_h) \, \mathbf{P} \right)

where :math:`\mathbf{L}_e = \sum_b K_e(a \to b)` and :math:`\mathbf{L}_h = \sum_j K_h(i \to j)` are the total state loss vectors.

This breakthrough reduces the computational scaling from :math:`O(N_{\mathrm{pairs}}^2)` to :math:`O(N_{\mathrm{occ}}^2 + N_{\mathrm{virt}}^2)`. A million exciton configurations are propagated in **less than 0.2 seconds per nuclear time step**.

2. Classical Path Approximation Surface Hopping (CPA-FSSH)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In CPA-FSSH, the classical nuclei follow a precomputed ground-state molecular dynamics trajectory :math:`\mathbf{R}(t)`. The quantum electronic wavefunction evolves according to the Time-Dependent Schrödinger Equation:

.. math::

   i\hbar \frac{d c_I(t)}{dt} = E_I(t) c_I(t) - i\hbar \sum_J d_{IJ}(t) c_J(t)

where :math:`c_I(t)` is the complex quantum amplitude of adiabatic state :math:`I`.

Tully's Fewest Switches Hopping Probability
"""""""""""""""""""""""""""""""""""""""""""

At each time step :math:`\Delta t`, an ensemble of classical trajectories is propagated. The probability for a trajectory currently residing on surface :math:`I` to switch to surface :math:`J` is given by Tully's formula:

.. math::

   g_{I \to J}(t) = \max\left( 0, \, \frac{-2 \Delta t \, \operatorname{Re}\left( c_I^*(t) c_J(t) d_{IJ}(t) \right)}{|c_I(t)|^2} \right) \times B_{IJ}(T)

A uniform random number :math:`\xi \in [0, 1]` is generated; if :math:`\sum_{K=1}^{J-1} g_{I \to K} < \xi \le \sum_{K=1}^J g_{I \to K}`, the trajectory hops to state :math:`J`.

---

5. Trajectory Precomputation & Wavefunction Tracking
-----------------------------------------------------

Numerical Non-Adiabatic Couplings (NAC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Along the classical nuclear trajectory :math:`\mathbf{R}(t)`, the non-adiabatic coupling is evaluated numerically via finite differences:

.. math::

   d_{IJ}(t + \frac{\Delta t}{2}) = \langle \psi_I(t) | \frac{\partial}{\partial t} | \psi_J(t) \rangle \approx \frac{S_{IJ}(t, t+\Delta t) - S_{JI}(t, t+\Delta t)}{2 \Delta t}

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

Hungarian Matching for Trivial Avoided Crossings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When nuclear vibrations bring two states close in energy, their adiabatic energy curves may cross. Sorting states strictly by instantaneous energy causes their physical identities to abruptly swap, introducing artificial spikes into :math:`d_{IJ}`.

``QDEX`` tracks states across time by solving the bipartite matching problem using the **Hungarian algorithm** on the cost matrix:

.. math::

   C_{IJ} = 1 - |S_{IJ}(t, t+\Delta t)|^2

This guarantees diabatic tracking and preserves the physical identity of frontier orbitals throughout the trajectory.

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

---

7. Phonon Spectral Density J(ω): Mapping Electron-Phonon Coupling
-----------------------------------------------------------------

Mathematical Definition
~~~~~~~~~~~~~~~~~~~~~~~

The **Phonon Spectral Density** :math:`J(\omega)` is the Fourier transform of the energy gap autocorrelation function:

.. math::

   J(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} C(t) \, e^{i \omega t} \, dt

In ``QDEX``, :math:`J(\omega)` is evaluated numerically using a Hann-windowed Fast Fourier Transform (FFT) of :math:`C(t)` and expressed in wavenumbers (:math:`\text{cm}^{-1}`).

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

      f_I = \frac{2}{3} \frac{m_e}{\hbar^2} E_I |\boldsymbol{\mu}_I|^2

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

3. **Trajectory Ensemble-Averaged Radiative Rate (:math:`\langle k_{\mathrm{rad}} \rangle_{\mathrm{MD}}`)**:
   When precomputed step data is available along the ab initio MD trajectory, the rate is averaged across all sampled nuclear configurations:

   .. math::

      \langle k_{\mathrm{rad}} \rangle_{\mathrm{MD}} = \frac{1}{N_{\mathrm{frames}}} \sum_{k=1}^{N_{\mathrm{frames}}} k_{\mathrm{rad}}^{\mathrm{therm}}(t_k)

Why Trajectory Averaging is Physically Crucial
""""""""""""""""""""""""""""""""""""""""""""""

* **Timescale Separation**: Carrier thermalization occurs on the femtosecond timescale (:math:`\tau_{\mathrm{cooling}} \sim 100 - 500\text{ fs}`), while radiative recombination takes nanoseconds (:math:`\tau_{\mathrm{rad}} \sim 1 - 50\text{ ns}`, over 10,000 times longer!). During this long waiting time, the nanocrystal explores its canonical thermal phase space.
* **Dynamic Symmetry Breaking & Herzberg-Teller Coupling**: At 0 K, high-symmetry nanocrystals (such as cubic perovskites or octahedral dots) often exhibit strictly dipole-forbidden dark ground excitons (:math:`f(0) = 0`). Thermal lattice vibrations dynamically break instantaneous inversion symmetry, mixing optically bright character into the ground exciton (*vibronic intensity borrowing*). Evaluating :math:`f` only at frame 0 would artificially predict an infinite radiative lifetime, whereas the trajectory ensemble average :math:`\langle f(t) \rangle_{\mathrm{MD}} > 0` correctly reproduces experimental nanosecond photoluminescence.

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

3. A Hann-windowed Fast Fourier Transform computes the power spectral density :math:`J(\omega)` in wavenumbers (:math:`\text{cm}^{-1}`):

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
3. **Electronic Coupling (:math:`V_{\mathrm{el}}`)**:
   In non-adiabatic transition theory, :math:`V_{\mathrm{el}}` is the effective off-diagonal electronic coupling between the initial excited state :math:`|1\rangle` and the ground state :math:`|0\rangle`. Via the Hellmann-Feynman theorem, the non-adiabatic coupling vector along nuclear velocity :math:`\dot{\mathbf{R}}` is:

   .. math::

      d_{10}(t) = \left\langle \psi_1(t) \middle| \frac{\partial}{\partial t} \middle| \psi_0(t) \right\rangle = \sum_A \dot{\mathbf{R}}_A \cdot \langle \psi_1 | \boldsymbol{\nabla}_A | \psi_0 \rangle

   In the time-derivative coupling representation, the effective electronic coupling matrix element is:

   .. math::

      V_{\mathrm{el}} = \hbar \, \langle |d_{10}(t)| \rangle

   where :math:`\langle |d_{10}| \rangle` is the trajectory-averaged non-adiabatic coupling magnitude between frontier orbitals, precomputed and stored in ``QDEX``'s step files.

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

10. CLI Flags & YAML Configuration Reference
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
   * - ``--namd-compact [dir]``
     - ``None``
     - Compress precomputed directory, eliminating redundant duplicate arrays.
   * - ``--namd-soc``
     - ``False``
     - Enable relativistic Spin-Orbit Coupling across NAMD precomputation.

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   namd:
     engine: "master_equation"          # "master_equation" (PME) or "surface_hopping" (CPA-FSSH)
     trajectory_dir: "./trajectory"      # Directory containing frame_* subdirectories
     dt_fs: 1.0                          # Nuclear MD time step in femtoseconds
     temperature_k: 300.0                # Lattice temperature for detailed balance
     tau_dec_fs: "cumulant"              # "cumulant" (ab initio), "edc", or fixed float in fs
     
     recombination:
       include_ground_state: true        # Couple excited manifold to ground state
       radiative: true                   # Enable Einstein spontaneous emission formula
       tau_nr_ns: 25.0                   # Defect trap non-radiative lifetime in nanoseconds
     
     storage:
       precompute_dir: "namd_precomputed"
       output_dir: "namd_results"
       compress: true
