Fewest-switches surface hopping (CPA-FSSH)
==========================================

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

* Module: ``qdex.namd.surface_hopping`` (``run_namd_dynamics``), ``qdex.namd.integrator``
* CLI: ``--namd-run`` (``--namd-method cpa_fssh``)
* YAML: ``namd.dynamics.method: cpa_fssh``, ``n_trajectories``, ``temperature_k``, ``detailed_balance``,
  ``nac_scheme``, ``seed``

Idea
----

Tully's fewest-switches surface hopping (J. Chem. Phys. 93, 1061 (1990)) replaces the superposition of
Ehrenfest dynamics by an ensemble of trajectories, each in one **active** state at a time. The amplitudes
:math:`c_I` are propagated with the time-dependent Schrödinger equation (:doc:`propagation`); a trajectory
jumps from its active state :math:`K` to :math:`J` with the probability that just keeps the fraction of
trajectories in each state equal to :math:`|c_J|^2`, with as few jumps as possible. The fraction of
trajectories in each state is the population.

In the classical path approximation (:doc:`introduction`) the hop does not change the nuclear trajectory:
upward hops are weighted with the Boltzmann factor instead of a velocity rescaling, and decoherence is
added explicitly.

In QDEX the active state of a trajectory is an exciton :math:`(i, a)`; the hole amplitudes
:math:`c^h_j` and the electron amplitudes :math:`c^e_b` are propagated in their own channels with the pair
energies at the active partner carrier, :math:`E_{jb}` with :math:`b = a` (hole channel) or :math:`j = i`
(electron channel), and the couplings :math:`d^{\mathrm{occ}*}` and :math:`d^{\mathrm{virt}}`
(:doc:`states_couplings`).

One nuclear step
----------------

1. **Propagation.** Both channels are propagated over the step with :math:`N_{\mathrm{sub}}` Strang
   sub-steps (:doc:`propagation`).

2. **Hopping probabilities.** Over the sub-steps the fewest-switches flux out of the active state is
   accumulated,

   .. math::

      g_{K \to J} = \sum_{m=1}^{N_{\mathrm{sub}}}
      \max\!\left(0,\; \frac{2\,\delta t\; \operatorname{Re}\!\big(c_K^*(\tau_m)\, c_J(\tau_m)\, d_{KJ}\big)}{|c_K(\tau_m)|^2}\right),

   the rate at which population flows from :math:`K` to :math:`J`
   (:math:`\dot{|c_K|^2} = -\sum_J 2\operatorname{Re}(c_K^* c_J d_{KJ})`) divided by the population of
   :math:`K`. Targets that are not pairs of the stored basis get zero; upward targets are multiplied by
   :math:`B_{KJ} = \min(1, e^{-(E_J - E_K)/k_BT})` with the pair energies at :math:`t_{k+1}` (detailed
   balance in place of velocity rescaling: Parandekar & Tully, J. Chem. Phys. 122, 094102 (2005)).

3. **Hop.** One random number decides between no hop, an electron hop (:math:`a \to b`, probabilities
   :math:`g^e`) and a hole hop (:math:`i \to j`, :math:`g^h`); if :math:`\sum g^e + \sum g^h > 1` both are
   scaled to a total of 1. A hop moves the active pair and **collapses** the amplitudes of that channel
   onto the new state (the trajectory has chosen its branch).

4. **Decoherence.** In a trajectory that did not hop, the inactive amplitudes are damped once per nuclear
   step with the pair dephasing times (:doc:`decoherence`),

   .. math::

      c_J \leftarrow c_J\, e^{-\Delta t/\tau_{KJ}} \quad (J \neq K), \qquad
      c_K \leftarrow \frac{c_K}{|c_K|} \sqrt{1 - \sum_{J \neq K} |c_J|^2},

   so that every coherence :math:`c_K c_J^*` decays at :math:`1/\tau_{KJ}` and the norm is kept. The
   amplitudes are **not** reset to the active state at every step: resetting would project the
   wavefunction every 2 fs and freeze the transitions (quantum Zeno effect), while keeping them without
   damping lets the trajectories stay coherent for ever (over-coherence) and hop back and forth.

Without ``decoherence_times.npz`` a warning is printed and all pairs get one time: ``tau_dec_fs`` (with
``decoherence: cumulant`` the lowest-exciton cumulant time) or 20 fs. ``method: cpa_fssh_gdc`` damps with :math:`e^{-(\Delta t/\tau)^2/2}` per step instead; applied
step after step without memory of the last collapse it is an exponential with a step-dependent rate, so
it is kept only for comparison.

Ensemble and initial conditions
-------------------------------

All trajectories of an origin are propagated together (amplitudes as columns of one matrix). The initial
active pairs are drawn from the bright pairs (oscillator strength) within the pump line
(``initial_excitation``: ``ratio_eg`` × lowest exciton, Gaussian of ``pulse_fwhm_ev``), and the amplitudes
start on them. With ``initial_conditions: multiple`` the run is repeated from ``n_origins`` frames along
the trajectory, ``window_fs`` each, and the results are averaged; ``seed`` makes the random numbers
reproducible.

Output: the populations of the pairs, the mean exciton energy and the excess energies of the electron
and of the hole above the instantaneous band edges in time (:doc:`analysis`).
