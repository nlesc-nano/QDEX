Fssh edc
========

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.surface_hopping``
* Callable: ``qdex.namd.surface_hopping.run_namd_dynamics``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   run_namd_dynamics(config)


2. Classical Path Approximation Surface Hopping with Energy-Based Decoherence (CPA-FSSH-EDC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In CPA-FSSH, the classical nuclei follow a precomputed ground-state molecular dynamics trajectory :math:`\mathbf{R}(t)`. The quantum electronic wavefunction :math:`|\Psi(t)\rangle` evolves according to the Time-Dependent Schrödinger Equation (TDSE):

.. math::

   i\hbar \frac{\partial |\Psi(t)\rangle}{\partial t} = \hat{H}_{\mathrm{elec}}(\mathbf{r}; \mathbf{R}(t)) |\Psi(t)\rangle

Expanding the electronic state in the instantaneous adiabatic Kohn-Sham or diagonal-BSE pair basis :math:`|\psi_I(\mathbf{R}(t))\rangle`:

.. math::

   |\Psi(t)\rangle = \sum_I c_I(t) |\psi_I(\mathbf{R}(t))\rangle

Substituting this expansion into the TDSE and projecting onto :math:`\langle \psi_I|` yields the coupled equations of motion for the complex amplitudes :math:`c_I(t)`:

.. math::

   i\hbar \frac{d c_I(t)}{dt} = E_I(t) c_I(t) - i\hbar \sum_J d_{IJ}(t) c_J(t)

where :math:`E_I(t) = \langle \psi_I | \hat{H}_{\mathrm{elec}} | \psi_I \rangle` is the instantaneous adiabatic energy and :math:`d_{IJ}(t) = \langle \psi_I | \frac{\partial}{\partial t} | \psi_J \rangle` is the non-adiabatic coupling matrix element.

In matrix notation, this forms the effective non-adiabatic Schrödinger equation:

.. math::

   \frac{d \mathbf{c}(t)}{dt} = -\frac{i}{\hbar} \mathbf{H}_{\mathrm{eff}}(t) \mathbf{c}(t), \qquad \mathbf{H}_{\mathrm{eff}}(t) = \operatorname{diag}\left(\mathbf{E}(t)\right) - i\hbar \mathbf{d}(t)

Because :math:`\mathbf{d}` is anti-Hermitian (:math:`d_{IJ} = -d_{JI}^*`), the off-diagonal coupling term :math:`-i\hbar \mathbf{d}` is **strictly Hermitian**, ensuring that :math:`\mathbf{H}_{\mathrm{eff}}` is Hermitian and the total electronic probability is conserved:

.. math::

   \sum_I |c_I(t)|^2 = 1

Unitary Time Evolution via Second-Order Strang Splitting
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To propagate the electronic coefficients across the nuclear time step :math:`[t_k, t_{k+1}]`, ``QDEX`` divides the interval into fine sub-steps :math:`\delta t = \Delta t_{\mathrm{nuc}} / N_{\mathrm{sub}}`. Within each sub-step, the propagator is evaluated using second-order Strang splitting:

.. math::

   \mathbf{U}(\delta t) = \exp\left( -i \mathbf{E}(\tau_m) \frac{\delta t}{2\hbar} \right) \exp\left( -\mathbf{d}(t_k) \, \delta t \right) \exp\left( -i \mathbf{E}(\tau_m) \frac{\delta t}{2\hbar} \right)

where :math:`\mathbf{E}(\tau_m)` is linearly interpolated between :math:`\mathbf{E}(t_k)` and :math:`\mathbf{E}(t_{k+1})`. Because :math:`\mathbf{d}` is anti-Hermitian, :math:`\exp(-\mathbf{d}\,\delta t)` is an exact unitary rotation, preserving norm conservation to machine precision.

Tully's Fewest Switches Hopping Probability
"""""""""""""""""""""""""""""""""""""""""""

In Fewest Switches Surface Hopping (Tully, *J. Chem. Phys.* 93, 1061, 1990), an ensemble of :math:`N_{\mathrm{traj}}` classical trajectories is propagated. Each trajectory :math:`tr` resides on a specific "active surface" :math:`K(t)`.

At each time step, the probability for a trajectory on active surface :math:`K` to switch to an inactive surface :math:`J \neq K` is determined by the outward probability flux accumulated across the electronic sub-steps:

.. math::

   g_{K \to J} = \sum_{m=1}^{N_{\mathrm{sub}}} \max\left( 0, \, \frac{2 \delta t \, \operatorname{Re}\left( c_K^*(\tau_m) c_J(\tau_m) d_{KJ}(t_k) \right)}{|c_K(\tau_m)|^2} \right) \times B_{KJ}(T)

where :math:`B_{KJ}(T)` enforces detailed balance at lattice temperature :math:`T`:

.. math::

   B_{KJ}(T) = \min\left( 1, \, \exp\left( -\frac{\max(E_J - E_K, 0)}{k_B T} \right) \right)

In standard molecular dynamics, detailed balance is maintained through momentum rescaling along the non-adiabatic coupling vector :math:`\mathbf{d}_{KJ}` (rejecting upward hops if nuclear kinetic energy is insufficient). In the Classical Path Approximation, nuclear velocities follow a fixed ground-state trajectory; multiplying upward hops by :math:`B_{KJ}(T)` provides the rigorous classical-path counterpart to velocity rescaling (Parandekar & Tully, *J. Chem. Phys.* 122, 094102, 2005; Jain, Alguire, & Subotnik, *J. Chem. Phys.* 144, 214110, 2016).

Hopping Decision & Wavepacket Collapse
""""""""""""""""""""""""""""""""""""""

A uniform random number :math:`\xi \in [0, 1)` is generated:

* If :math:`\sum_{L=1}^{J-1} g_{K \to L} < \xi \le \sum_{L=1}^J g_{K \to L}`, trajectory :math:`tr` hops to state :math:`J`.
* The active surface is reassigned: :math:`K \leftarrow J`.
* The electronic wavepacket undergoes **projective collapse** onto the new active state:

  .. math::

     c_J \leftarrow 1.0, \qquad c_{L \neq J} \leftarrow 0.0

This projective collapse models the immediate decoherence of the trajectory as it separates along the newly occupied potential surface.

The Overcoherence Problem & Continuous Energy-Based Decoherence (EDC)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In standard FSSH, when a trajectory does **not** hop, electronic amplitudes :math:`c_J` remain populated indefinitely because classical nuclei in the CPA follow a shared trajectory. In dense quantum-dot manifolds with hundreds of states, this causes artificial state superpositions to linger, producing unphysical hopping loops (**overcoherence**).

Conversely, if the wavepacket is reset to a pure state at every nuclear step (:math:`c_K = 1, c_{J \neq K} = 0`), frequent projective measurements collapse the wavepacket before transitions can develop, completely freezing carrier relaxation—a pathological numerical artifact known as the **Quantum Zeno effect**.

To resolve both problems, ``QDEX`` implements the **Energy-Based Decoherence (EDC)** scheme (Granucci & Persico, *J. Chem. Phys.* 126, 134114, 2007). 

In EDC, electronic wavepacket amplitudes persist continuously across nuclear steps. When no hop occurs, inactive states :math:`J \neq K` are continuously damped over the nuclear step :math:`\Delta t`:

.. math::

   c_J(t + \Delta t) \leftarrow c_J(t + \Delta t) \, \exp\left( -\frac{\Delta t}{\tau_{KJ}} \right), \qquad \forall J \neq K

where :math:`\tau_{KJ}` is the state-pair pure-dephasing time. To strictly preserve total probability conservation (:math:`\sum_L |c_L|^2 = 1`) without altering the quantum phase of the active state, the active amplitude :math:`c_K` is renormalized:

.. math::

   c_K(t + \Delta t) \leftarrow \sqrt{ 1 - \sum_{J \neq K} |c_J(t + \Delta t)|^2 } \; \frac{c_K(t + \Delta t)}{|c_K(t + \Delta t)|}

This continuous damping smoothly suppresses off-diagonal coherences on physical dephasing timescales (:math:`10 - 40\text{ fs}`) while allowing short-time quantum interference to drive hopping.

State-Pair Dephasing Times: Ab Initio Covariance vs. Granucci-Persico
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``QDEX`` supports two rigorous evaluations of :math:`\tau_{KJ}`:

1. **Ab Initio Pure-Dephasing Covariance Matrix**:
   Evaluated from the trajectory energy fluctuations:

   .. math::

      \tau_{KJ} = \min\left( \tau_{\max}, \, \frac{\hbar \sqrt{2}}{\sigma_{KJ}} \right), \qquad \sigma_{KJ}^2 = \operatorname{Var}\left( E_K(t) - E_J(t) \right)

   Diagonal elements are capped at :math:`\tau_{\max} = 500\text{ fs}`. Precomputed across all orbital pairs and cached in ``decoherence_times.npz``.

2. **Granucci-Persico Energy-Gap Formula**:
   Evaluated instantaneously from the energy gap and nuclear kinetic energy:

   .. math::

      \tau_{KJ} = \frac{\hbar}{|E_K - E_J|} \left( 1 + \frac{C}{E_{\mathrm{kin}}} \right)

   where :math:`C = 0.1\text{ Hartree} = 2.721\text{ eV}`, and :math:`E_{\mathrm{kin}} = \frac{3}{2} N_{\mathrm{atoms}} k_B T` is the classical thermal kinetic energy of the nuclear lattice.

