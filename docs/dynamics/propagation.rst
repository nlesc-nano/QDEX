Propagating the electrons: nuclear steps and electronic sub-steps
=================================================================

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

* Module: ``qdex.namd.integrator`` (``propagate_channel_batch_strang``, ``accumulate_fssh_flux``)
* YAML: ``namd.integration.integrator``, ``namd.integration.n_substeps``

Two time scales
---------------

The nuclei move on the scale of vibrational periods (50–300 cm\ :sup:`−1`: 100–600 fs), and the MD
trajectory stores a frame every :math:`\Delta t` = 1–2 fs. The electronic amplitudes oscillate with the
Bohr frequencies of the energy differences in the window, :math:`2\pi\hbar/\Delta E` ≈ 1–4 fs for
:math:`\Delta E` = 1–4 eV, and they rotate into each other at the rate :math:`|d|`, up to
0.45 fs\ :sup:`−1` in the valence band of the CsPbX\ :sub:`3` dots. The electronic equation is therefore
integrated with sub-steps inside each nuclear step.

Nuclear step
------------

The nuclear time step :math:`\Delta t` is the spacing of the MD frames (``namd.trajectory.dt_nuc_fs``,
2 fs for the CsPbX\ :sub:`3` trajectories). Everything that depends on the geometry is known only at the
frames: the pair energies :math:`E(t_k)` and the step coupling :math:`d_k` (from the overlap of frames
:math:`k` and :math:`k+1`, :doc:`states_couplings`). Hops, DISH decoherence events, the FSSH decoherence
damping and the PME rates act once per nuclear step.

Electronic sub-steps
--------------------

The amplitudes of the electron and of the hole channel are propagated through the step with
:math:`N_{\mathrm{sub}}` sub-steps (``namd.integration.n_substeps``),

.. math::

   \delta t = \frac{\Delta t}{N_{\mathrm{sub}}}, \qquad
   \tau_m = \left(m + \tfrac{1}{2}\right)\delta t, \quad m = 0 \ldots N_{\mathrm{sub}} - 1 .

Inside the step the energies are interpolated linearly between the frames and the coupling is held at its
step value (the generator of the step for ``nac_scheme: logm``):

.. math::

   E_I(\tau_m) = E_I(t_k) + \frac{\tau_m}{\Delta t}\left[E_I(t_{k+1}) - E_I(t_k)\right], \qquad
   \mathbf{H}_{\mathrm{eff}}(\tau_m) = \operatorname{diag}\, \mathbf{E}(\tau_m) - i\hbar\, \mathbf{d}_k .

:math:`\mathbf{d}` is anti-Hermitian, so :math:`\mathbf{H}_{\mathrm{eff}}` is Hermitian. In the electron
channel :math:`E_I` are the pair energies at the hole of the trajectory and :math:`\mathbf{d}` the virtual
couplings; in the hole channel the pair energies at the electron and :math:`\mathbf{d}^{\mathrm{occ}*}`.

Each sub-step is a second-order Strang splitting (``integrator: strang``),

.. math::

   \mathbf{c}(\tau + \delta t) = D_m\, O\, D_m\, \mathbf{c}(\tau), \qquad
   D_m = e^{-i\, \mathbf{E}(\tau_m)\, \delta t / 2\hbar}, \qquad O = e^{-\mathbf{d}_k\, \delta t},

exactly unitary. :math:`O` is computed once per nuclear step and applied with one matrix product per
sub-step to all trajectories at once (columns of :math:`\mathbf{c}`). With :math:`\mathbf{E}` constant the
sub-steps compose to :math:`e^{-\mathbf{d}_k \Delta t} = U^\dagger` for the logm coupling: the populations
move exactly as the overlap of the two frames says. The splitting error comes from the energy phases and
the coupling not commuting, :math:`\mathcal{O}\big((\Delta E\,\delta t/\hbar)(|d|\,\delta t)\big)` per
sub-step; production runs use :math:`N_{\mathrm{sub}} = 20` (δt = 0.1 fs; default 2, see the
convergence test below).

``rk4`` (fourth-order Runge–Kutta with renormalisation) is available for tests;
``step_unitary_matrix_exp`` (exact :math:`e^{-i\mathbf{H}_{\mathrm{eff}}\delta t/\hbar}` of one wavefunction)
is not used on the ensemble path.

What happens at the end of each nuclear step
--------------------------------------------

* **FSSH** (:doc:`fssh`): the hopping probabilities accumulated over the sub-steps decide a hop; a hop
  collapses the amplitudes, otherwise the inactive amplitudes are damped with the dephasing times.
* **DISH** (:doc:`dish`): decoherence events are drawn for every inactive state; an event collapses onto
  that state or removes it.
* **PME** (:doc:`pme`): no amplitudes; the populations take 20 row-stochastic sub-steps with the rates of
  the step.

.. _namd-substeps:

Sub-step convergence
--------------------

Tested on the spin-free CsPbBr\ :sub:`3` cube (1300 × 800 window, 2 fs frames, logm couplings, one
origin, 200 fs, 1000 trajectories, pump 2 E\ :sub:`g`, same random seed). Change of the mean exciton
energy (eV):

.. list-table::
   :header-rows: 1

   * - :math:`N_{\mathrm{sub}}` (δt)
     - DISH, 50 fs
     - DISH, 200 fs
     - FSSH, 50 fs
     - FSSH, 200 fs
     - wall time (FSSH, 16 threads)
   * - 5 (0.4 fs)
     - −0.0227
     - −0.1895
     - −0.0332
     - −0.1338
     - 1.9 min
   * - 20 (0.1 fs)
     - −0.0227
     - −0.1895
     - −0.0333
     - −0.1324
     - 4.0 min
   * - 50 (0.04 fs)
     - −0.0227
     - −0.1895
     - −0.0282
     - −0.1298
     - 8.8 min

DISH does not change with the sub-steps (its hops are decided once per nuclear step from
:math:`|c_J|^2`). FSSH changes by 3–4 meV of 130 meV at 200 fs, within the scatter of its hops (the
amplitudes differ slightly, so the same random numbers give different hops). :math:`N_{\mathrm{sub}}` = 20
is converged; 10 is enough for the more expensive SOC runs.
