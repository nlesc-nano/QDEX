Part 6: Carrier Cooling Dynamics & Photoluminescence (NAMD)
===========================================================

Following the photoexcitation of a semiconductor nanocrystal or quantum dot by an ultrashort laser pulse, high-energy ("hot") electrons and holes rapidly dissipate their excess energy through electron-phonon scattering and non-adiabatic transitions. Carriers cascade down the dense ladder of excited states, cooling toward the band edges before recombining to the ground state.

``miniBSE`` features an advanced, high-throughput **Non-Adiabatic Molecular Dynamics (NAMD)** engine designed to simulate carrier relaxation, phonon bottleneck phenomena, surface defect trapping/de-trapping, and photoluminescence recombination along *ab initio* molecular dynamics (AIMD) trajectories.

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

2. Pauli Master Equation (PME) vs. CPA-FSSH: When to Use Which?
---------------------------------------------------------------

A central methodological decision in non-adiabatic dynamics is choosing between a **deterministic Master Equation** and **stochastic Fewest Switches Surface Hopping (FSSH)**. Both frameworks are implemented in ``miniBSE``, and each possesses distinct physical domains of applicability:

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

3. Theoretical Foundations of the Dynamical Engines
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

4. Trajectory Precomputation & Wavefunction Tracking
-----------------------------------------------------

Numerical Non-Adiabatic Couplings (NAC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Along the classical nuclear trajectory :math:`\mathbf{R}(t)`, the non-adiabatic coupling is evaluated numerically via finite differences:

.. math::

   d_{IJ}(t + \frac{\Delta t}{2}) = \langle \psi_I(t) | \frac{\partial}{\partial t} | \psi_J(t) \rangle \approx \frac{S_{IJ}(t, t+\Delta t) - S_{JI}(t, t+\Delta t)}{2 \Delta t}

where :math:`S_{IJ}(t, t+\Delta t) = \langle \psi_I(t) | \psi_J(t+\Delta t) \rangle` is the cross-frame state overlap. In ``miniBSE``, the underlying atomic orbital cross-overlaps :math:`S_{\mu \nu}(t, t+\Delta t) = \int \chi_\mu(\mathbf{r}; \mathbf{R}(t)) \chi_\nu(\mathbf{r}; \mathbf{R}(t+\Delta t)) d\mathbf{r}` are evaluated analytically via Libint2.

Eliminating Gauge Phase Discontinuities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard electronic eigensolvers determine eigenvectors up to an arbitrary global phase factor :math:`e^{i \theta_I(t)}`. If left uncorrected, random sign flips between successive MD frames cause :math:`S_{II}(t, t+\Delta t) \approx -1`, producing spurious non-adiabatic couplings that are orders of magnitude too large.

``miniBSE`` eliminates gauge discontinuities by applying a phase rotation:

.. math::

   |\psi_I(t+\Delta t)\rangle \leftarrow |\psi_I(t+\Delta t)\rangle \, e^{-i \theta_I}

where :math:`\theta_I = \operatorname{arg}(S_{II}(t, t+\Delta t))`. This guarantees that the diagonal overlap is strictly real and positive:

.. math::

   \operatorname{Re}(S_{II}(t, t+\Delta t)) \ge 0, \quad \operatorname{Im}(S_{II}(t, t+\Delta t)) = 0

Hungarian Matching for Trivial Avoided Crossings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When nuclear vibrations bring two states close in energy, their adiabatic energy curves may cross. Sorting states strictly by instantaneous energy causes their physical identities to abruptly swap, introducing artificial spikes into :math:`d_{IJ}`.

``miniBSE`` tracks states across time by solving the bipartite matching problem using the **Hungarian algorithm** on the cost matrix:

.. math::

   C_{IJ} = 1 - |S_{IJ}(t, t+\Delta t)|^2

This guarantees diabatic tracking and preserves the physical identity of frontier orbitals throughout the trajectory.

---

5. Electronic Decoherence: Origin, Computation, and Rationale
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

Rather than guessing an empirical decoherence time (e.g. :math:`\tau_{\mathrm{dec}} = 10\text{ fs}`), ``miniBSE`` computes :math:`\tau_{\mathrm{dec}}` *ab initio* from the second-order cumulant expansion of energy gap fluctuations.

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

``miniBSE`` solves for the characteristic decoherence time :math:`\tau_{\mathrm{dec}}` satisfying:

.. math::

   D(\tau_{\mathrm{dec}}) = \frac{1}{e}

Why We Compute It This Way
""""""""""""""""""""""""""

1. **Parameter-Free**: Eliminates arbitrary empirical fitting parameters from NAMD simulations.
2. **Temperature & Lattice Sensitive**: Soft, anharmonic lattices (such as lead halide perovskites) exhibit large thermal gap fluctuations (:math:`\sigma_E \approx 50 - 100\text{ meV}`), correctly yielding short dephasing times (:math:`\tau_{\mathrm{dec}} \approx 7 - 12\text{ fs}`), whereas rigid covalent quantum dots (like InAs or Si) yield longer dephasing times (:math:`\tau_{\mathrm{dec}} \approx 20 - 40\text{ fs}`).

---

6. Phonon Spectral Density J(ω): Mapping Electron-Phonon Coupling
-----------------------------------------------------------------

Mathematical Definition
~~~~~~~~~~~~~~~~~~~~~~~

The **Phonon Spectral Density** :math:`J(\omega)` is the Fourier transform of the energy gap autocorrelation function:

.. math::

   J(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} C(t) \, e^{i \omega t} \, dt

In ``miniBSE``, :math:`J(\omega)` is evaluated numerically using a Hann-windowed Fast Fourier Transform (FFT) of :math:`C(t)` and expressed in wavenumbers (:math:`\text{cm}^{-1}`).

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

7. Radiative & Non-Radiative Recombination Mechanisms
-----------------------------------------------------

Once carriers have relaxed to the band edges (forming the lowest 1S exciton), they recombine to the ground state :math:`|S_0\rangle` through competing radiative and non-radiative channels:

.. math::

   \frac{d P_I(t)}{dt} = \sum_{J \neq I} \left[ k_{J \to I} P_J - k_{I \to J} P_I \right] - \left( k_{I \to 0}^{\mathrm{rad}} + k_{I \to 0}^{\mathrm{nr}} \right) P_I

1. Ab Initio Einstein Radiative Rate (:math:`k_{\mathrm{rad}}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spontaneous photon emission from state :math:`|I\rangle` into the vacuum radiation field inside a dielectric medium is evaluated from the **Einstein A coefficient**:

.. math::

   k_{I \to 0}^{\mathrm{rad}} = \left[ \frac{2 e^2}{4\pi \epsilon_0 m_e c^3 \hbar^2} \right] n_{\mathrm{r}} \, E_I^2 \, f_I

where:
* :math:`E_I` is the emission transition energy (in eV).
* :math:`f_I = \frac{2}{3} \frac{m_e}{\hbar^2} E_I |\boldsymbol{\mu}_I|^2` is the *ab initio* dimensionless oscillator strength.
* :math:`n_{\mathrm{r}}` is the optical refractive index of the semiconductor material, loaded from ``REFRACTIVE_INDEX_DICT`` in ``hardness.py`` (e.g. :math:`n_{\mathrm{r}} = 2.3` for :math:`\text{CsPbBr}_3`, :math:`n_{\mathrm{r}} = 3.5` for :math:`\text{InAs}`).

2. Non-Radiative Multi-Phonon Decay Across Large Gaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Direct non-radiative recombination across the fundamental band gap (:math:`E_g > 1.5\text{ eV}`) requires dissipating a massive amount of electronic energy into the nuclear bath. Because this corresponds to dozens of vibrational quanta (:math:`p = E_g / \hbar \omega_{\mathrm{LO}} \sim 50`), perturbation theory in nuclear kinetic energy yields the **Englman-Jortner Energy Gap Law**:

.. math::

   k_{\mathrm{nr}} = A_{\mathrm{nr}} \exp\left( -\gamma \frac{E_g}{\hbar \omega_{\mathrm{LO}}} \right)

where:
* :math:`\hbar \omega_{\mathrm{LO}}` is the characteristic optical phonon energy.
* :math:`\gamma = \ln\left( \frac{E_g}{S \, \hbar \omega_{\mathrm{LO}}} \right) - 1` is the electronic-vibrational coupling parameter, with :math:`S` being the Huang-Rhys factor.
* :math:`A_{\mathrm{nr}} \approx 10^{13}\text{ s}^{-1}` is the electronic pre-exponential factor.

Because :math:`k_{\mathrm{nr}}` decreases exponentially with increasing band gap, pristine, defect-free quantum dots exhibit negligible band-to-band non-radiative decay, resulting in near-unity intrinsic luminescence.

3. Intermediate & Small Gaps (Multi-Phonon Nuclear Bath)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the gap narrows (e.g. in infrared quantum dots like HgTe or PbS with :math:`E_g < 0.8\text{ eV}`), multi-phonon emission accelerates dramatically. In this regime where explicit AIMD trajectories cannot sample rare multi-phonon tunneling, a **Franck-Condon Weighted Density of States (FCWD)** model can be employed:

.. math::

   k_{\mathrm{nr}} = \frac{2\pi}{\hbar} |V_{\mathrm{el}}|^2 \, \mathrm{FCWD}(E_g)

where :math:`\mathrm{FCWD}(E_g) = \frac{1}{\sqrt{2\pi \sigma^2}} \exp\left( -\frac{(E_g - \lambda)^2}{2\sigma^2} \right)` accounts for the nuclear reorganization energy :math:`\lambda`.

4. Defect Trap-Assisted Recombination
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In real nanocrystals with unpassivated surfaces, non-radiative recombination is dominated by **Shockley-Read-Hall (SRH) defect trapping**. Deep mid-gap trap states break the large band gap into smaller sub-steps, dramatically accelerating non-radiative relaxation. In ``miniBSE``, trap-assisted recombination is modeled by specifying an effective trap lifetime:

.. code-block:: yaml

   namd:
     recombination:
       tau_nr_ns: 25.0  # Trap recombination lifetime in nanoseconds

5. Photoluminescence Quantum Yield (PLQY)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total Photoluminescence Quantum Yield is evaluated from the branching ratio between radiative emission and non-radiative loss:

.. math::

   \mathrm{PLQY} = \frac{\langle k_{\mathrm{rad}} \rangle}{\langle k_{\mathrm{rad}} \rangle + k_{\mathrm{nr}}} \times 100\%

---

8. In-Depth Analysis of NAMD Simulations
----------------------------------------

``miniBSE`` includes a dedicated analysis module (``miniBSE.namd.analysis``) that automatically processes precomputed and dynamic trajectory data.

1. Carrier Cooling Curves & Lifetime Fitting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Carrier relaxation is monitored by tracking the time-dependent excess energy of electrons (:math:`\Delta E_e`) and holes (:math:`\Delta E_h`) above their respective band edges:

.. math::

   \Delta E_e(t) = \sum_{a \in \mathrm{virt}} P_a(t) \left( \varepsilon_a(t) - \varepsilon_{\mathrm{LUMO}}(t) \right)

.. math::

   \Delta E_h(t) = \sum_{i \in \mathrm{occ}} P_i(t) \left( \varepsilon_{\mathrm{HOMO}}(t) - \varepsilon_i(t) \right)

The cooling lifetimes :math:`\tau_e` and :math:`\tau_h` are extracted by fitting the excess energy decay to an exponential function:

.. math::

   \Delta E(t) \approx \Delta E(0) \, \exp\left( -\frac{t}{\tau_{\mathrm{cooling}}} \right)

2. State-Resolved Population Kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Transient populations :math:`P_I(t)` are exported to ``carrier_cooling_populations.csv``, displaying the sequential decay of initial hot excitons into intermediate states and finally into the emitting :math:`1S` state.

3. NAC vs. Energy Gap Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To verify whether non-adiabatic transitions obey the energy-gap law, ``miniBSE`` samples pairs of states across trajectory frames and plots non-adiabatic coupling magnitudes :math:`|d_{IJ}|` against energy differences :math:`|E_J - E_I|`. This distinguishes smooth exponential decay from resonant vibronic enhancements.

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

9. CLI Flags & YAML Configuration Reference
-------------------------------------------

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
