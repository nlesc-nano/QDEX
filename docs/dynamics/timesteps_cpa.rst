Timesteps cpa
=============

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.precompute``
* Callable: ``qdex.namd.precompute.precompute_namd_data``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   precompute_namd_data(config)


2. Multi-Timescale Integration: Separating Nuclear and Electronic Time Steps
----------------------------------------------------------------------------

A fundamental challenge in simulating non-adiabatic carrier dynamics is the dramatic **timescale mismatch** between nuclear vibrations and electronic phase oscillations.


The Timescale Mismatch
~~~~~~~~~~~~~~~~~~~~~~

* **Nuclear Motion** (:math:`\sim 1 - 2\text{ fs}`): 
  Atomic nuclei are thousands of times heavier than electrons (:math:`M_{\mathrm{Pb}} / m_e \approx 3.8 \times 10^5`). Nuclear motion is governed by acoustic and optical phonon frequencies (:math:`\omega_{\mathrm{ph}} \approx 50 - 300\text{ cm}^{-1}`), corresponding to vibrational periods of :math:`T_{\mathrm{vib}} \approx 100 - 600\text{ fs}`. A nuclear time step of :math:`\Delta t_{\mathrm{nuc}} \approx 1.0 - 2.0\text{ fs}` is therefore fully sufficient to integrate classical Newton's equations of motion with energy conservation.

* **Electronic Wavefunction Oscillations** (:math:`\sim 0.005 - 0.05\text{ fs}`):
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

