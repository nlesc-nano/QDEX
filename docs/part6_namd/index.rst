Part 6: Carrier Cooling Dynamics & Photoluminescence (NAMD)
===========================================================

Following photoexcitation of a semiconductor nanocrystal by an ultrashort laser pulse, high-energy ("hot") electrons and holes rapidly dissipate their excess energy via electron-phonon scattering and non-adiabatic transitions, cooling down to the band edges before recombining to the ground state.

``miniBSE`` features a comprehensive, state-of-the-art **Non-Adiabatic Molecular Dynamics (NAMD)** engine designed to simulate hot carrier relaxation, non-radiative trapping, and photoluminescence across *ab initio* molecular dynamics (AIMD) trajectories.

---

1. The NAMD Pipeline Architecture
---------------------------------

The NAMD simulation workflow is decoupled into three distinct stages:

.. code-block:: text

   AIMD Trajectory Frames (frame_0001 ... frame_N)
                 │
                 ▼
   [Stage 1: Trajectory Precomputation]
     ├─ Evaluate Quasiparticle & Diagonal BSE Exciton States per frame
     ├─ Compute Cross-Frame Overlaps S(t, t+Δt) via Libint2
     ├─ Extract Non-Adiabatic Couplings (NAC) d_IJ(t)
     ├─ Correct Arbitrary Spinor Gauge Phases e^{iθ}
     ├─ Track Trivial State Inversions via the Hungarian Algorithm
     └─ Compress & Cache Precomputed Data into step_*.npz
                 │
                 ▼
   [Stage 2: Dynamics & Kinetics Simulation]
     ├─ Engine A: Pauli Master Equation (PME) [Deterministic, Tensorized]
     ├─ Engine B: CPA-FSSH [Stochastic Trajectory Surface Hopping]
     ├─ Ab Initio Cumulant Decoherence (from Energy Gap Fluctuations)
     └─ Ground-State Recombination (Einstein Emission k_rad & Energy Gap Law k_nr)
                 │
                 ▼
   [Stage 3: Analysis & Visualization]
     ├─ Carrier Cooling Lifetimes (Electron vs. Hole Relaxation)
     ├─ Transient Population Kinetics (Hot Exciton -> Band Edge 1S)
     ├─ Non-Adiabatic Coupling Matrices & Phonon Spectral Densities J(ω)
     └─ Photoluminescence Quantum Yield (PLQY) Calculation

---

2. Stage 1: Trajectory Precomputation
--------------------------------------

At each nuclear time step :math:`t_n` along the molecular dynamics trajectory, ``miniBSE`` computes the electronic or excitonic states :math:`|\psi_I(t_n)\rangle`.

Cross-Frame Overlaps & Non-Adiabatic Couplings (NAC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The non-adiabatic coupling vector :math:`\mathbf{d}_{IJ} = \langle \psi_I | \boldsymbol{\nabla}_{\mathbf{R}} | \psi_J \rangle` drives transitions between adiabatic surfaces. Along a classical nuclear trajectory :math:`\mathbf{R}(t)`, the time-derivative coupling is:

.. math::

   d_{IJ}(t) = \langle \psi_I(t) | \frac{\partial}{\partial t} | \psi_J(t) \rangle = \mathbf{d}_{IJ} \cdot \dot{\mathbf{R}}(t)

In ``miniBSE``, :math:`d_{IJ}(t)` is evaluated numerically via central finite differences from the cross-frame overlap matrix :math:`S_{IJ}(t, t+\Delta t) = \langle \psi_I(t) | \psi_J(t+\Delta t) \rangle`:

.. math::

   d_{IJ}(t + \frac{\Delta t}{2}) \approx \frac{S_{IJ}(t, t+\Delta t) - S_{JI}(t, t+\Delta t)}{2 \Delta t}

The cross-frame atomic orbital overlap matrix :math:`S_{\mu \nu}(t, t+\Delta t) = \int \chi_\mu(\mathbf{r}; \mathbf{R}(t)) \chi_\nu(\mathbf{r}; \mathbf{R}(t+\Delta t)) \, d\mathbf{r}` is evaluated analytically using the Libint2 engine.

Geometric Gauge Phase Alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because numerical eigensolvers determine electronic wavefunctions up to an arbitrary gauge phase :math:`e^{i\theta_I(t)}`, raw cross-frame overlaps exhibit random phase discontinuities that artificially inflate non-adiabatic couplings by orders of magnitude.

``miniBSE`` rigorously enforces phase continuity by rotating each state:

.. math::

   |\psi_I(t+\Delta t)\rangle \leftarrow |\psi_I(t+\Delta t)\rangle \, e^{-i \theta_I}

where the phase angle :math:`\theta_I = \operatorname{arg}\left( S_{II}(t, t+\Delta t) \right)` is chosen such that the diagonal overlap is strictly real and positive:

.. math::

   \operatorname{Re}\left( S_{II}(t, t+\Delta t) \right) \ge 0, \quad \operatorname{Im}\left( S_{II}(t, t+\Delta t) \right) = 0

Hungarian Matching for Trivial Avoided Crossings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When two states with different spatial symmetries approach each other during nuclear vibrations, their adiabatic eigenvalues may cross or exhibit a trivial avoided crossing. Sorting purely by energy causes their state identities to suddenly swap.

To maintain physical state continuity across time, ``miniBSE`` solves the linear assignment problem using the **Hungarian algorithm** on the assignment cost matrix:

.. math::

   C_{IJ} = 1 - |S_{IJ}(t, t+\Delta t)|^2

This guarantees that diabatic character is preserved across crossings.

---

3. Stage 2: Dynamical Engines
-----------------------------

``miniBSE`` provides two complementary simulation engines for propagating non-adiabatic dynamics:

Option A: Pauli Master Equation (PME)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For large quantum dots where the excitonic manifold contains hundreds of thousands of states, stochastic surface hopping is computationally impractical. ``miniBSE`` provides a deterministic Master Equation for state populations :math:`P_I(t)`:

.. math::

   \frac{d P_I(t)}{dt} = \sum_{J \neq I} \left[ k_{J \to I}(t) P_J(t) - k_{I \to J}(t) P_I(t) \right]

The instantaneous state-to-state transition rate :math:`k_{I \to J}` is formulated via Fermi's Golden Rule with a Lorentzian line shape representing quantum decoherence:

.. math::

   k_{I \to J}(t) = 2 |d_{IJ}(t)|^2 \, \left[ \frac{\tau_{\mathrm{dec}}}{1 + \left( \frac{(E_J(t) - E_I(t)) \tau_{\mathrm{dec}}}{\hbar} \right)^2} \right] \times B_{IJ}(T)

where :math:`B_{IJ}(T)` enforces thermodynamic **detailed balance** at temperature :math:`T`:

.. math::

   B_{IJ}(T) = \begin{cases}
     1 & \text{if } E_J \ge E_I \text{ (downward cooling)} \\
     \exp\left( -\frac{E_I - E_J}{k_B T} \right) & \text{if } E_J < E_I \text{ (upward thermal activation)}
   \end{cases}

Vectorized Tensor Decomposition
"""""""""""""""""""""""""""""""

Under the Diagonal BSE representation, an exciton state :math:`|ia\rangle` relaxes through independent electron (:math:`a \to b`) and hole (:math:`i \to j`) scattering channels. The rate equations decompose into a matrix product:

.. math::

   \frac{\partial \mathbf{P}}{\partial t} = \left( \mathbf{P} \mathbf{K}_e^T - \mathbf{P} \operatorname{diag}(\mathbf{L}_e) \right) + \left( \mathbf{K}_h \mathbf{P} - \operatorname{diag}(\mathbf{L}_h) \mathbf{P} \right)

where :math:`\mathbf{P} \in \mathbb{R}^{N_{\mathrm{occ}} \times N_{\mathrm{virt}}}` is the 2D exciton population matrix, :math:`\mathbf{K}_e` is the virtual electron rate matrix, and :math:`\mathbf{K}_h` is the occupied hole rate matrix.

This tensor decomposition reduces the computational complexity from :math:`O(N_{\mathrm{pairs}}^2)` down to :math:`O(N_{\mathrm{occ}}^2 + N_{\mathrm{virt}}^2)`. A system with :math:`10^6` exciton pairs propagates in **less than 0.2 seconds per time step**.

Option B: Classical Path Approximation FSSH (CPA-FSSH)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For smaller systems or when quantum coherence is essential, ``miniBSE`` supports stochastic surface hopping. The electronic wavepacket evolves along the classical path according to the Time-Dependent Schrödinger Equation:

.. math::

   i\hbar \frac{d c_I}{dt} = E_I(t) c_I(t) - i\hbar \sum_J d_{IJ}(t) c_J(t)

Trajectories hop stochastically between adiabatic states according to Tully's Fewest Switches hopping probability:

.. math::

   g_{I \to J} = \max\left( 0, \, \frac{-2 \Delta t \, \operatorname{Re}\left( c_I^* c_J d_{IJ} \right)}{|c_I|^2} \right) \times B_{IJ}(T)

---

4. Automated Ab Initio Cumulant Decoherence
--------------------------------------------

The transition rates in non-adiabatic dynamics depend sensitively on the electronic dephasing (decoherence) time :math:`\tau_{\mathrm{dec}}`. Rather than treating :math:`\tau_{\mathrm{dec}}` as an empirical fitting parameter, ``miniBSE`` calculates it *ab initio* from the second-order cumulant expansion of energy gap fluctuations.

Let :math:`\delta E_1(t) = E_1(t) - \langle E_1 \rangle` represent the fluctuation of the lowest excited state energy along the trajectory. The normalized energy gap autocorrelation function is:

.. math::

   C(t) = \frac{\langle \delta E_1(0) \delta E_1(t) \rangle}{\sigma_E^2}

where :math:`\sigma_E^2 = \langle \delta E_1^2 \rangle`. The line shape function :math:`g(t)` is obtained by double integration:

.. math::

   g(t) = \frac{\sigma_E^2}{\hbar^2} \int_0^t dt_1 \int_0^{t_1} dt_2 \, C(t_2)

The electronic dephasing decay function is:

.. math::

   D(t) = \exp\left( -g(t) \right)

``miniBSE`` numerically solves for the characteristic decoherence time :math:`\tau_{\mathrm{dec}}` satisfying:

.. math::

   D(\tau_{\mathrm{dec}}) = \frac{1}{e}

Enabling this via ``namd.tau_dec_fs: "cumulant"`` guarantees that carrier cooling rates are entirely parameter-free.

---

5. Ground-State Recombination & Photoluminescence
-------------------------------------------------

Once carriers cool to the band edge (1S exciton), they recombine to the ground state :math:`|S_0\rangle` through competing radiative and non-radiative channels:

.. math::

   \frac{d P_I}{dt} = \sum_{J \neq I} \left( k_{J \to I} P_J - k_{I \to J} P_I \right) - \left( k_{I \to 0}^{\mathrm{rad}} + k_{I \to 0}^{\mathrm{nr}} \right) P_I

1. Ab Initio Einstein Radiative Rate (:math:`k_{\mathrm{rad}}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spontaneous emission rate from excited state :math:`|I\rangle` to the ground state is computed directly from its *ab initio* transition dipole moment :math:`\boldsymbol{\mu}_I` and oscillator strength :math:`f_I`:

.. math::

   k_{I \to 0}^{\mathrm{rad}} = \left[ \frac{2 e^2}{4\pi \epsilon_0 m_e c^3 \hbar^2} \right] n_{\mathrm{r}} \, E_I^2 \, f_I

where:
* :math:`E_I` is the emission photon energy.
* :math:`n_{\mathrm{r}}` is the optical refractive index of the semiconductor crystal (e.g. :math:`n_{\mathrm{r}} = 2.3` for :math:`\text{CsPbBr}_3`), retrieved from ``REFRACTIVE_INDEX_DICT`` in ``hardness.py``.

2. Multi-Phonon Non-Radiative Rate (:math:`k_{\mathrm{nr}}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Direct non-radiative ground-state recombination across the large optical gap (:math:`> 1.5\text{ eV}`) proceeds via multi-phonon emission described by the **Englman-Jortner Energy Gap Law**:

.. math::

   k_{\mathrm{nr}} = A_{\mathrm{nr}} \exp\left( -\gamma \frac{E_g}{\hbar \omega_{\mathrm{LO}}} \right)

where :math:`\hbar \omega_{\mathrm{LO}}` is the longitudinal optical phonon energy and :math:`\gamma` is the electronic-vibrational coupling parameter. Alternatively, users can define an explicit non-radiative lifetime (e.g. ``tau_nr_ns: 20.0``) representing surface trap-assisted recombination.

3. Photoluminescence Quantum Yield (PLQY)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total Photoluminescence Quantum Yield (PLQY) is computed as the competition between ensemble-averaged radiative emission and non-radiative loss:

.. math::

   \mathrm{PLQY} = \frac{\langle k_{\mathrm{rad}} \rangle}{\langle k_{\mathrm{rad}} \rangle + k_{\mathrm{nr}}} \times 100\%

---

6. CLI Flags & YAML Configuration Reference
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
     - Run Stage 1 NAMD precomputation (cross-frame overlaps, NACs, phase tracking, and caching).
   * - ``--namd-run``
     - ``False``
     - Run Stage 2 NAMD carrier cooling dynamics from precomputed cache.
   * - ``--namd-compact [dir]``
     - ``None``
     - Compress and clean redundant arrays in the precomputed directory.
   * - ``--namd-soc``
     - ``False``
     - Enable relativistic Spin-Orbit Coupling for NAMD precomputation.

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   namd:
     engine: "master_equation"        # "master_equation" (PME) or "surface_hopping" (FSSH)
     dt_fs: 1.0                       # MD time step in femtoseconds
     temperature_k: 300.0             # Lattice temperature for detailed balance
     tau_dec_fs: "cumulant"           # "cumulant" (ab initio) or fixed float in fs
     
     recombination:
       include_ground_state: true
       radiative: true                # Uses Einstein spontaneous emission formula
       tau_nr_ns: 25.0                # Non-radiative lifetime in nanoseconds
     
     storage:
       precompute_dir: "namd_precomputed"
       output_dir: "namd_results"
