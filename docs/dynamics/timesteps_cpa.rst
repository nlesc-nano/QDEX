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


Nuclear steps and electronic sub-steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Nuclear step.** The nuclear time step :math:`\Delta t` is the spacing of the MD frames
(``namd.trajectory.dt_nuc_fs``, 2 fs for the CsPbX\ :sub:`3` trajectories). Everything that depends on
the geometry is known only at the frames: the pair energies :math:`E(t_k)` and the step coupling
:math:`d_k` (from the overlap of frames :math:`k` and :math:`k+1`, :doc:`nacs_tracking`). Hops, DISH
decoherence events, the FSSH decoherence damping and the PME rates act once per nuclear step.

**Electronic sub-steps.** The amplitudes of the electron and of the hole channel are propagated through
the step with :math:`N_{\mathrm{sub}}` sub-steps (``namd.integration.n_substeps``),

.. math::

   \delta t = \frac{\Delta t}{N_{\mathrm{sub}}}, \qquad
   \tau_m = \left(m + \tfrac{1}{2}\right)\delta t, \quad m = 0 \ldots N_{\mathrm{sub}} - 1 .

Inside the step the energies are interpolated linearly between the frames and the coupling is held at its
step value (the generator of the step for ``nac_scheme: logm``):

.. math::

   E_I(\tau_m) = E_I(t_k) + \frac{\tau_m}{\Delta t}\left[E_I(t_{k+1}) - E_I(t_k)\right], \qquad
   \mathbf{H}_{\mathrm{eff}}(\tau_m) = \operatorname{diag}\, \mathbf{E}(\tau_m) - i\hbar\, \mathbf{d}_k .

:math:`\mathbf{d}` is anti-Hermitian, so :math:`\mathbf{H}_{\mathrm{eff}}` is Hermitian. Each sub-step is a
second-order Strang splitting (``integrator: strang``),

.. math::

   \mathbf{c}(\tau + \delta t) = D_m\, O\, D_m\, \mathbf{c}(\tau), \qquad
   D_m = e^{-i\, \mathbf{E}(\tau_m)\, \delta t / 2\hbar}, \qquad O = e^{-\mathbf{d}_k\, \delta t},

exactly unitary; :math:`O` is computed once per nuclear step and applied with one matrix product per
sub-step to all trajectories at once (columns of :math:`\mathbf{c}`). With :math:`\mathbf{E}` constant the
sub-steps compose to :math:`e^{-\mathbf{d}_k \Delta t} = U^\dagger` for the logm coupling: the populations
move exactly as the overlap of the two frames says. The splitting error comes from the energy phases and
the coupling not commuting, :math:`\mathcal{O}\big((\Delta E\,\delta t/\hbar)\,(|d|\,\delta t)\big)` per
sub-step: with :math:`|d|\,\Delta t` up to 0.9 in the valence band of the CsPbX\ :sub:`3` dots and pair-energy
spreads of a few eV, :math:`N_{\mathrm{sub}} = 20` (δt = 0.1 fs) is used for production (default 2;
convergence: :ref:`namd-substeps`).

``step_unitary_matrix_exp`` (single wavefunction, exact :math:`e^{-i\mathbf{H}_{\mathrm{eff}}\delta t/\hbar}`
from a diagonalization) is not used on the ensemble path; ``rk4`` is available for tests.

Tully hopping flux (CPA-FSSH)
"""""""""""""""""""""""""""""

The fewest-switches probability of the active state :math:`I` is accumulated over the sub-steps,

.. math::

   g_{I \to J} = \sum_{m} \max\left(0,\; \frac{2\,\delta t\, \operatorname{Re}\left(c_I^*(\tau_m)\, c_J(\tau_m)\, d_{IJ}\right)}{|c_I(\tau_m)|^2}\right),

and an upward hop is accepted with the Boltzmann factor :math:`B_{IJ} = \min(1, e^{-(E_J - E_I)/k_BT})`,
the classical-path stand-in for velocity rescaling (Parandekar & Tully, J. Chem. Phys. 2005). One hop
(electron or hole) per nuclear step; a hop collapses the amplitudes onto the new state. Trajectories that
do not hop are damped once per nuclear step, :math:`c_J \leftarrow c_J\, e^{-\Delta t/\tau_{IJ}}` for
:math:`J \ne I` with the state-pair dephasing times :math:`\tau_{IJ}` (``decoherence_times.npz``), and
renormalised. The amplitudes are **not** reset to the active state at every step.

DISH
""""

DISH (``method: dish``, Jaeger, Fischer & Prezhdo, J. Chem. Phys. 137, 22A545 (2012)) propagates the same
amplitudes without flux hops. Once per nuclear step a decoherence event for each inactive state :math:`J`
occurs with probability :math:`1 - e^{-\Delta t/\tau_{IJ}}`; on an event the trajectory hops to :math:`J`
with probability :math:`|c_J|^2 B_{IJ}` (only to stored pairs) and collapses onto it, otherwise
:math:`c_J` is set to zero and the amplitudes renormalised (:doc:`dish`).


Electronic sub-stepping in the Pauli master equation (PME)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The PME propagates populations with golden-rule rates built once per nuclear step from :math:`|d_k|^2`,
the energy gaps and the line widths (state-pair dephasing times :math:`\tau_{IJ}` with
``pme_tau: pairs``, :doc:`pme`). For the diagonal-BSE pair manifold each of 20 sub-steps applies a
row-stochastic matrix :math:`T = I + \delta t\, K` of the electron channel and of the hole channel
(:math:`\mathbf{P} \leftarrow \mathbf{P} T_e`, then :math:`\mathbf{P} \leftarrow T_h^\mathsf{T} \mathbf{P}`);
a row whose leaving probability would exceed 1 is renormalised onto its outgoing transitions, which keeps
every population non-negative and conserves the total. Recombination :math:`e^{-k_{\mathrm{loss}}\Delta t}`
is applied once per nuclear step.


.. _namd-substeps:

Sub-step convergence
~~~~~~~~~~~~~~~~~~~~

Tested on the spin-free CsPbBr\ :sub:`3` cube (1300 × 800 window, 2 fs frames, logm couplings, one
origin, 200 fs, 1000 trajectories, pump 2 E\ :sub:`g`): see the table below (filled in from
``namd/CsPbBr3/conv`` on NHR).
