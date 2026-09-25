Comparison
==========

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.precompute``
* Callable: ``qdex.namd.precompute.precompute_namd_data``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   precompute_namd_data(config)


4. Comparative Synthesis: PME vs. CPA-FSSH-EDC vs. DISH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following comprehensive comparison synthesizes the mathematical foundations, computational characteristics, and physical applicability of the three dynamical engines in ``QDEX``:

.. list-table::
   :widths: 20 26 27 27
   :header-rows: 1

   * - Theoretical Dimension
     - Pauli Master Equation (PME)
     - CPA-FSSH with EDC
     - Decoherence-Induced Surface Hopping (DISH)
   * - **Primary Dynamical Variable**
     - Real-valued population vector :math:`\mathbf{P}(t) \in \mathbb{R}^{N}`.
     - Complex amplitude vector :math:`\mathbf{c}(t) \in \mathbb{C}^{N}` plus active surface index :math:`K(t)`.
     - Complex amplitude vector :math:`\mathbf{c}(t) \in \mathbb{C}^{N}` plus active surface index :math:`K(t)`.
   * - **Electronic Equation of Motion**
     - Markovian rate equation: :math:`\frac{d\mathbf{P}}{dt} = \mathbf{K} \mathbf{P}`.
     - Continuous TDSE: :math:`i\hbar \dot{\mathbf{c}} = \mathbf{H}_{\mathrm{eff}} \mathbf{c}` with continuous EDC damping.
     - Piecewise unitary TDSE: :math:`\mathbf{c}(t+\Delta t) = \mathbf{U} \mathbf{c}(t)` (no continuous damping).
   * - **Hopping / Transition Mechanism**
     - Continuous probability flux between populations via FGR rate matrix.
     - Tully's Fewest Switches non-adiabatic coupling flux: :math:`g_{K \to J} \propto \operatorname{Re}(c_K^* c_J d_{KJ})`.
     - Environment-driven stochastic branching: Poisson dephasing followed by Born-rule collapse :math:`P_{\mathrm{hop}} \propto |c_J|^2`.
   * - **Decoherence Treatment**
     - Implicit: Lorentzian broadening of energy conservation delta function.
     - Continuous exponential damping of inactive states: :math:`c_J \leftarrow c_J e^{-\Delta t / \tau_{KJ}}`.
     - Discrete stochastic quenching of inactive states: :math:`c_J \leftarrow 0` upon rejected dephasing.
   * - **Detailed Balance Enforcement**
     - Transition rates scaled by Boltzmann factor :math:`B_{IJ}(T) = \exp(-\Delta E / k_B T)`.
     - Hopping fluxes scaled by Boltzmann factor :math:`B_{IJ}(T)` (classical-path velocity rescaling stand-in).
     - Hopping probabilities scaled by Boltzmann factor :math:`B_{IJ}(T)`.
   * - **Quantum Zeno Vulnerability**
     - **Immune**: Propagates macroscopic populations; no wavepacket resetting.
     - **Immune**: Wavepacket amplitudes persist across steps; EDC damps smoothly without abrupt collapse.
     - **Immune**: Decoherence events are stochastic and Poisson-distributed; active state remains untouched.
   * - **Equivalence in Dense Limit**
     - Canonical baseline (Lorentzian FGR).
     - Diverges if dephasing times are miscalibrated or couplings are stiff.
     - **Analytically proves equivalence** to Lorentzian FGR in the fast-dephasing limit.
   * - **Computational Complexity**
     - :math:`O(N_{\mathrm{occ}}^2 + N_{\mathrm{virt}}^2)` via BLAS-3 tensor contraction (< 0.2 s per step for :math:`10^6` pairs).
     - :math:`O(N_{\mathrm{traj}} \cdot N_{\mathrm{sub}} \cdot N_{\mathrm{dyn}}^2)` (heavy; requires fine sub-stepping).
     - :math:`O(N_{\mathrm{traj}} \cdot N_{\mathrm{sub}} \cdot N_{\mathrm{dyn}}^2)` (vectorized batching across trajectories).
   * - **Recommended Regime**
     - Ultrafast screening across huge state manifolds (:math:`> 1,000` states), high DOS.
     - Discrete frontier states (:math:`1P \to 1S`), persistent quantum coherence, small systems.
     - Nanocrystals with dense bands, avoided crossings, defect trapping, and single-carrier dwell times.


