Dish
====

Part of :doc:`/dynamics/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.namd.integrator``
* Callable: ``qdex.namd.integrator.step_dish_batch``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   step_dish_batch(C, active_surfaces, E_batch, dt_fs, tau_mat=None, beta=None, detailed_balance=True, min_tau_fs=1.0)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part6_namd/index.rst:414-522``

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

* **Case A: Hop Accepted** (:math:`R_2 < P_{\mathrm{hop}, J}`):
  The trajectory successfully transitions to state :math:`J`. The active surface switches :math:`K \leftarrow J`, and the wavepacket undergoes complete projective collapse onto state :math:`J`:

  .. math::

     c_J \leftarrow 1.0, \qquad c_{L \neq J} \leftarrow 0.0

  If multiple states simultaneously qualify for a hop in a single time step, one state is selected with probability proportional to :math:`P_{\mathrm{hop}, J}`.

* **Case B: Hop Rejected** (:math:`R_2 \ge P_{\mathrm{hop}, J}`):
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
