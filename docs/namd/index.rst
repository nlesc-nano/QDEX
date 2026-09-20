Non-Adiabatic Dynamics (NAMD)
==============================

The `miniBSE` NAMD engine simulates non-adiabatic excited-state relaxation, carrier cooling, and photoluminescence recombination in semiconductor nanostructures across molecular dynamics trajectories.

NAMD Pipeline Architecture
--------------------------

.. code-block:: text

   Trajectory MD Frames (frame_*)
                 │
                 ▼
   [Stage 1: Precomputation]
     ├─ Quasiparticle & Diagonal BSE Exciton States
     ├─ Cross-Frame Overlaps S(t, t+Δt) & Non-Adiabatic Couplings (NAC)
     ├─ Spinor Phase Correction & Hungarian Crossing Tracking
     └─ Compact Caching (step_*.npz)
                 │
                 ▼
   [Stage 2: Dynamics & Kinetics]
     ├─ Option A: Pauli Master Equation (PME) [Deterministic, Tensorized]
     ├─ Option B: Classical Path Approximation FSSH (CPA-FSSH) [Stochastic]
     ├─ Ab Initio Cumulant Decoherence from Energy Gap Fluctuations
     └─ Recombination to Ground State (Einstein Spontaneous Emission & PLQY)
                 │
                 ▼
   [Stage 3: Analysis & Visualization]
     ├─ Carrier Cooling Curves (Electron vs Hole Lifetimes)
     ├─ State-Resolved Populations (1Se, 1Sh, 1S Exciton)
     ├─ NAC vs Energy-Gap Law & Phonon Spectral Density J(ω)
     └─ 6-Panel Publication Figures & Interactive Plotly Dashboards

Stage 1: Trajectory Precomputation
----------------------------------

At each nuclear time step $\Delta t$, the non-adiabatic coupling (NAC) $d_{IJ} = \langle \psi_I | \frac{\partial}{\partial t} | \psi_J \rangle$ is evaluated numerically from the cross-frame overlap matrix:

.. math::

   d_{IJ}(t) \approx \frac{S_{IJ}(t, t+\Delta t) - S_{JI}(t, t+\Delta t)}{2 \Delta t}

Phase Alignment & Hungarian Tracking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because electronic eigensolvers determine wavefunctions up to an arbitrary gauge phase :math:`e^{i\theta}`, phase discontinuities must be eliminated. `miniBSE` enforces:

1. **Phase Correction**: Rotates each state by :math:`e^{-i\theta_k}` such that :math:`\operatorname{Re}(S_{kk}) \ge 0`.
2. **Hungarian Tracking**: At trivial crossings where adiabatic states invert, states are tracked across time using the Hungarian matching algorithm on the cost matrix :math:`C_{ij} = 1 - |S_{ij}|^2`.

Stage 2: Dynamics Frameworks
----------------------------

1. Pauli Master Equation (PME)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For dense manifolds containing millions of states, `miniBSE` uses a deterministic Master Equation:

.. math::

   \frac{d P_I(t)}{dt} = \sum_{J \neq I} \left[ k_{J \to I}(t) P_J(t) - k_{I \to J}(t) P_I(t) \right]

The transition rate $k_{I \to J}$ is given by Fermi's Golden Rule with Lorentzian decoherence broadening:

.. math::

   k_{I \to J} = 2 |d_{IJ}|^2 \, \left[ \frac{\tau_{\mathrm{dec}}}{1 + \left( \frac{\Delta E_{IJ} \tau_{\mathrm{dec}}}{\hbar} \right)^2} \right] \times B_{IJ}(T)

where :math:`B_{IJ}(T) = \min\left(1, e^{-(E_J - E_I)/k_B T}\right)` enforces detailed balance.

**Vectorized Tensor Decomposition**:
In diagonal BSE, exciton states :math:`|ia\rangle` relax via independent electron (:math:`a \to b`) and hole (:math:`i \to j`) channels:

.. math::

   \frac{\partial \mathbf{P}}{\partial t} = \left( \mathbf{P} \mathbf{K}_e^T - \mathbf{P} \operatorname{diag}(\mathbf{L}_e) \right) + \left( \mathbf{K}_h \mathbf{P} - \operatorname{diag}(\mathbf{L}_h) \mathbf{P} \right)

This reduces the complexity from :math:`O(N_{\mathrm{pairs}}^2)` to :math:`O(N_{\mathrm{occ}}^2 + N_{\mathrm{virt}}^2)`, propagating millions of exciton pairs in $<0.2\text{ s}$ per step.

2. Classical Path Approximation FSSH (CPA-FSSH)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, `miniBSE` supports stochastic trajectory surface hopping. The electronic wavepacket evolves via the Time-Dependent Schrödinger Equation:

.. math::

   i\hbar \frac{d c_I}{dt} = E_I c_I - i\hbar \sum_J d_{IJ} c_J

and hops between surfaces with Tully's fewest-switches probability:

.. math::

   g_{I \to J} = \max\left( 0, \, \frac{-2 \Delta t \, \operatorname{Re}(c_I^* c_J d_{IJ})}{|c_I|^2} \right) \times B_{IJ}(T)

Automated Cumulant Decoherence
------------------------------

Rather than guessing an empirical decoherence time :math:`\tau_{\mathrm{dec}}`, `miniBSE` computes it *ab initio* from the second-order cumulant expansion of the lowest excited state gap fluctuations :math:`\delta E_1(t) = E_1(t) - \langle E_1 \rangle`:

.. math::

   C(t) = \frac{\langle \delta E_1(0) \delta E_1(t) \rangle}{\sigma^2}, \quad g(t) = \frac{\sigma^2}{\hbar^2} \int_0^t dt_1 \int_0^{t_1} dt_2 \, C(t_2)

The dephasing function is :math:`D(t) = \exp(-g(t))`, and :math:`\tau_{\mathrm{dec}}` is solved such that :math:`D(\tau_{\mathrm{dec}}) = 1/e`. Set `tau_dec_fs: "cumulant"` to enable.

Recombination & Photoluminescence
---------------------------------

`miniBSE` unifies femtosecond cooling with nanosecond photoluminescence by coupling the excited states to the ground state :math:`|S_0\rangle`:

.. math::

   \frac{d P_I}{dt} = \sum_{J \neq I} (k_{J \to I} P_J - k_{I \to J} P_I) - \left( k_{I \to 0}^{\mathrm{rad}} + k_{I \to 0}^{\mathrm{nr}} \right) P_I

1. **Radiative Rate ($k_{\mathrm{rad}}$)**:
   Computed from *ab initio* transition dipoles via Einstein's spontaneous emission formula:

   .. math::

      k_{I \to 0}^{\mathrm{rad}} = \left[ \frac{2 e^4}{4\pi\epsilon_0 m_e c^3 \hbar^2} \right] n_{\mathrm{r}} \, E_I^2 \, f_I

   where $n_{\mathrm{r}}$ is the material's optical refractive index from `REFRACTIVE_INDEX_DICT`.
2. **Non-Radiative Rate ($k_{\mathrm{nr}}$)**:
   Evaluates multi-phonon emission via the Englman-Jortner Energy Gap Law:

   .. math::

      k_{\mathrm{nr}} = A_{\mathrm{nr}} \exp\left( -\gamma \frac{E_g}{\hbar \omega_{\mathrm{LO}}} \right)

   or incorporates surface defect trap lifetimes via `tau_nr_ns`.
3. **Photoluminescence Quantum Yield (PLQY)**:

   .. math::

      \mathrm{PLQY} = \frac{\langle k_{\mathrm{rad}} \rangle}{\langle k_{\mathrm{rad}} \rangle + k_{\mathrm{nr}}} \times 100\%
