Decoherence
===========

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.integrator``
* Callable: ``qdex.namd.integrator.apply_edc_decoherence``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   apply_edc_decoherence(c, active_surface, E_vec, dt_elec, c_param=0.1, n_atoms=1, temp_k=300.0, decay_type='exponential')


6. Electronic Decoherence: Origin, Computation, and Rationale
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

Rather than guessing an empirical decoherence time (e.g. :math:`\tau_{\mathrm{dec}} = 10\text{ fs}`), ``QDEX`` computes :math:`\tau_{\mathrm{dec}}` *ab initio* from the second-order cumulant expansion of energy gap fluctuations.

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

``QDEX`` solves for the characteristic decoherence time :math:`\tau_{\mathrm{dec}}` satisfying:

.. math::

   D(\tau_{\mathrm{dec}}) = \frac{1}{e}

Why We Compute It This Way
""""""""""""""""""""""""""

1. **Parameter-Free**: Eliminates arbitrary empirical fitting parameters from NAMD simulations.
2. **Temperature & Lattice Sensitive**: Soft, anharmonic lattices (such as lead halide perovskites) exhibit large thermal gap fluctuations (:math:`\sigma_E \approx 50 - 100\text{ meV}`), correctly yielding short dephasing times (:math:`\tau_{\mathrm{dec}} \approx 7 - 12\text{ fs}`), whereas rigid covalent quantum dots (like InAs or Si) yield longer dephasing times (:math:`\tau_{\mathrm{dec}} \approx 20 - 40\text{ fs}`).


State-Pair Pure-Dephasing Matrices (tau_ij)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In multi-state quantum dynamics, electronic dephasing times vary substantially across orbital pairs depending on energy gaps, orbital spatial localization, and lattice phonon coupling. Rather than assuming an arbitrary uniform dephasing time, ``QDEX`` precomputes the complete state-pair pure-dephasing matrices :math:`\boldsymbol{\tau}_{\mathrm{occ}} \in \mathbb{R}^{N_{\mathrm{occ}} \times N_{\mathrm{occ}}}` and :math:`\boldsymbol{\tau}_{\mathrm{virt}} \in \mathbb{R}^{N_{\mathrm{virt}} \times N_{\mathrm{virt}}}` directly from the AIMD trajectory energy fluctuations.

Microscopic Derivation from Short-Time Cumulant Expansion
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For any pair of adiabatic states :math:`(i, j)`, the instantaneous energy difference along the trajectory is:

.. math::

   \Delta \varepsilon_{ij}(t) = \varepsilon_i(t) - \varepsilon_j(t)

The thermal fluctuation around its equilibrium trajectory mean is:

.. math::

   \delta \Delta \varepsilon_{ij}(t) = \Delta \varepsilon_{ij}(t) - \langle \Delta \varepsilon_{ij} \rangle

with variance:

.. math::

   \sigma_{ij}^2 = \operatorname{Var}\left( \Delta \varepsilon_{ij} \right) = \langle \left( \delta \Delta \varepsilon_{ij}(t) \right)^2 \rangle

In second-order cumulant expansion of the reduced density operator, the pure-dephasing decay function :math:`D_{ij}(t)` is driven by the bath line-shape function :math:`g_{ij}(t)`:

.. math::

   D_{ij}(t) = \exp\left( -g_{ij}(t) \right), \qquad g_{ij}(t) = \frac{1}{\hbar^2} \int_0^t dt_1 \int_0^{t_1} dt_2 \, \langle \delta \Delta \varepsilon_{ij}(0) \delta \Delta \varepsilon_{ij}(t_2) \rangle

On timescales shorter than the characteristic nuclear phonon correlation time (:math:`t \ll \tau_{\mathrm{bath}} \sim 50\text{ fs}`), the energy fluctuation autocorrelation is essentially static: :math:`\langle \delta \Delta \varepsilon_{ij}(0) \delta \Delta \varepsilon_{ij}(t_2) \rangle \approx \sigma_{ij}^2`.

The double time integral simplifies analytically to:

.. math::

   g_{ij}(t) \approx \frac{\sigma_{ij}^2}{\hbar^2} \int_0^t dt_1 \, t_1 = \frac{\sigma_{ij}^2 \, t^2}{2 \hbar^2}

yielding a Gaussian dephasing profile:

.. math::

   D_{ij}(t) \approx \exp\left( -\frac{\sigma_{ij}^2 \, t^2}{2 \hbar^2} \right) = \exp\left( - \left(\frac{t}{\tau_{ij}}\right)^2 \right)

Defining the characteristic dephasing time :math:`\tau_{ij}` at the standard :math:`1/e` decay threshold (:math:`D_{ij}(\tau_{ij}) = 1/e`):

.. math::

   \frac{\sigma_{ij}^2 \, \tau_{ij}^2}{2 \hbar^2} = 1 \implies \tau_{ij} = \frac{\hbar \sqrt{2}}{\sigma_{ij}}

To prevent numerical singularities along the diagonal (:math:`i = j` where :math:`\sigma_{ii} = 0`), the diagonal elements are set to a maximum physical cutoff :math:`\tau_{\max} = 500\text{ fs}`:

.. math::

   \tau_{ij} = \min\left( \tau_{\max}, \, \frac{\hbar \sqrt{2}}{\sigma_{ij}} \right)

The resulting pairwise dephasing matrices are precomputed via ``qdex --namd-decoherence`` and cached in ``decoherence_times.npz``, where both **CPA-FSSH-EDC** and **DISH** dynamically consume them at run time.


Surface Hopping Schemes: FSSH-EDC vs. DISH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To prevent the Quantum Zeno effect (where resetting wavepackets at every nuclear step freezes carrier cooling) while curing FSSH overcoherence, ``QDEX`` provides two physically rigorous stochastic algorithms:

1. **Continuous CPA-FSSH with Energy-Based Decoherence Correction (EDC)**
   (Granucci & Persico, *J. Chem. Phys.* 126, 134114, 2007)
   Electronic amplitudes :math:`\mathbf{C}(t)` persist continuously across nuclear steps. At each nuclear step:
   * Electronic TDSE is integrated over the nuclear interval to accumulate Tully's Fewest Switches flux :math:`g_{I \to J}`.
   * If a hop occurs: the active state switches to :math:`J` and the wavepacket collapses onto state :math:`J` (:math:`c_J = 1`, :math:`c_{K \neq J} = 0`).
   * If no hop occurs: off-diagonal amplitudes are continuously damped by the state-pair dephasing factor:
     
     .. math::
     
        c_J \leftarrow c_J \exp\left( -\frac{\Delta t}{\tau_{IJ}} \right), \quad \forall J \neq I
        
     and the active amplitude :math:`c_I` is renormalized to conserve total probability :math:`\sum_K |c_K|^2 = 1`.

2. **Decoherence-Induced Surface Hopping (DISH)**
   (Jaeger, Fischer, Prezhdo, *J. Chem. Phys.* 137, 22A545, 2012; Akimov & Prezhdo, *J. Chem. Phys.* 138, 124102, 2013)
   In DISH, electronic propagation is strictly unitary between stochastic collapse events (no continuous artificial damping). At each nuclear step:
   * For each state :math:`J \neq I`, a Poisson dephasing event occurs with probability:
     
     .. math::
     
        P_{\mathrm{dec}, J} = 1 - \exp\left( -\frac{\Delta t}{\tau_{IJ}} \right)
        
   * If a dephasing event occurs for state :math:`J`, the system attempts a collapse (hop) to state :math:`J` with probability:
     
     .. math::
     
        P_{\mathrm{hop}, J} = |c_J|^2 \times \min\left( 1, \, e^{-\beta \max(E_J - E_I, 0)} \right)
        
   * If accepted, the trajectory hops to state :math:`J` and collapses into a pure state (:math:`c_J = 1`).
   * If rejected, the coherence with state :math:`J` is quenched (:math:`c_J \to 0`) and the remaining states are renormalized.

In dense semiconductor nanocrystals, DISH analytically recovers the Pauli Master Equation cooling rates while maintaining individual trajectory statistics.

