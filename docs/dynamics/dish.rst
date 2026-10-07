Decoherence-induced surface hopping (DISH)
==========================================

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

* Callable: ``qdex.namd.integrator.step_dish_batch``
* CLI: ``--namd-run`` (``--namd-method dish``)
* YAML: ``namd.dynamics.method: dish`` (needs ``decoherence_times.npz``, :doc:`decoherence`)

Idea
----

In DISH (Jaeger, Fischer & Prezhdo, J. Chem. Phys. 137, 22A545 (2012)) the hops are not driven by the
coupling flux but by decoherence. The amplitudes evolve with the time-dependent Schrödinger equation,
without damping; at random times, set by the dephasing time of each pair, a state decoheres from the
active one. At that moment the wavefunction is measured: it collapses onto that state with probability
:math:`|c_J|^2` (a hop), or that state is removed from the superposition. In a dense manifold with fast
dephasing, this is the regime of carrier cooling in dots.

One nuclear step
----------------

1. **Propagation** of both channels over the step (:doc:`propagation`), with no damping and no flux.

2. **Decoherence events.** For every inactive state :math:`J` of each channel an event occurs with
   probability

   .. math::

      P^{\mathrm{dec}}_J = 1 - e^{-\Delta t/\tau_{KJ}} ,

   with :math:`K` the active state and :math:`\tau_{KJ}` the pair dephasing time.

3. **Collapse or removal.** For each decohered :math:`J` the trajectory hops to :math:`J` with probability

   .. math::

      P^{\mathrm{hop}}_J = |c_J|^2\; B_{KJ}, \qquad B_{KJ} = \min\!\big(1, e^{-(E_J - E_K)/k_BT}\big),

   only to pairs of the stored basis. On a hop the amplitudes collapse onto :math:`J` (if several states
   qualify, one is drawn with weights :math:`P^{\mathrm{hop}}`); otherwise :math:`c_J \leftarrow 0` and the
   amplitudes are renormalised. The electron channel is decided first (at the current hole), then the hole
   channel at the new electron.

Relation to the master equation
-------------------------------

For two states starting in :math:`K`, the amplitude in :math:`J` builds up between decoherence events,
:math:`|c_J(s)|^2 \approx |d_{KJ}|^2 s^2/[1 + (\Delta E\, s/2\hbar)^2]` for a coherence of age :math:`s`, and
an event converts it into a hop. Averaged over the Poisson-distributed ages
(:math:`\langle s^2 \rangle = 2\tau^2`) the rate is

.. math::

   k_{K \to J}^{\mathrm{DISH}} = \frac{2\,|d_{KJ}|^2\, \tau_{KJ}}{1 + (\Delta E_{KJ}\, \tau_{KJ}/\hbar)^2}\; B_{KJ},

the golden-rule rate with a Lorentzian of width :math:`\hbar/\tau_{KJ}` used by the PME (:doc:`pme`). Two-
level tests of ``step_dish_batch`` reproduce it within 5–8 % (4000 trajectories, :math:`\Delta E` =
0–0.1 eV, :math:`\tau` = 5–10 fs) and detailed balance within 1 % (:math:`p_1/p_0` = 0.310 against 0.313).
DISH therefore gives the PME rates when both use the same :math:`\tau_{KJ}` (``pme_tau: pairs``), but
keeps single trajectories: dwell times, branching, trapping and de-trapping events.

Differences from the original formulation
-----------------------------------------

* The decoherence time of state :math:`J` is the pair time :math:`\tau_{KJ}` with the active state, not
  the amplitude-weighted :math:`1/\tau_J = \sum_k |c_k|^2/\tau_{Jk}` of Jaeger *et al.*; the two agree when
  the active amplitude is close to 1, which frequent decoherence events ensure.
* The active state itself does not decohere.
* Several states can decohere in one step (each with its own probability).
