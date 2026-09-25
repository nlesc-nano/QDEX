Auger matrix elements
=====================

Part of :doc:`/recombination/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.auger``
* Callable: ``qdex.auger.compute_auger_matrix_element``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   compute_auger_matrix_element(q_recomb: np.ndarray, q_eject_1: np.ndarray, q_eject_2: np.ndarray, W_resta: np.ndarray, q_recomb_alt: Optional[np.ndarray]=None)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part8_auger/index.rst:73-81``

2. Many-Body Formulation & Matrix Elements
------------------------------------------

Within first-order time-dependent perturbation theory (**Fermi's Golden Rule**), the transition rate from an initial multi-carrier state :math:`|i\rangle` to a continuum manifold of final states :math:`|f\rangle` is:

.. math::

   \Gamma_{\mathrm{Auger}} = \frac{2\pi}{\hbar} \sum_f \left| M_{if} \right|^2 \delta(E_i - E_f)


.. rubric:: From ``docs/part8_auger/index.rst:82-110``

The Two-Body Screened Coulomb Operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The interaction driving the transition is the screened two-body Coulomb operator:

.. math::

   \hat{W} = \frac{1}{2} \sum_{pqrs} \langle pq | \hat{W} | rs \rangle \, c_p^\dagger c_q^\dagger c_s c_r

For the negative trion / biexciton :math:`eeh` channel:
* Initial 3-carrier state: :math:`|i\rangle = c_{e_1}^\dagger c_{e_2}^\dagger c_h |0\rangle`
* Final 1-carrier state: :math:`|f\rangle = c_{e'}^\dagger |0\rangle`

Evaluating the matrix element :math:`M_{if} = \langle f | \hat{W} | i \rangle` yields an anti-symmetrized pair of direct and exchange two-electron integrals:

.. math::

   M_{if} = V_{e' e_1, h e_2} - V_{e' e_2, h e_1}

where each two-electron orbital integral is defined as:

.. math::

   V_{e' e_1, h e_2} = \iint \phi_{e'}^*(\mathbf{r}_1) \, \phi_{e_1}(\mathbf{r}_1) \; W(\mathbf{r}_1, \mathbf{r}_2) \; \phi_h^*(\mathbf{r}_2) \, \phi_{e_2}(\mathbf{r}_2) \; d\mathbf{r}_1 \, d\mathbf{r}_2

Notice the physical factorization of the two coordinate spaces:
1. :math:`\rho_{\mathrm{recomb}}(\mathbf{r}_2) = \phi_h^*(\mathbf{r}_2) \phi_{e_2}(\mathbf{r}_2)` is the **recombination transition density** of the annihilating electron-hole pair.
2. :math:`\rho_{\mathrm{eject}}(\mathbf{r}_1) = \phi_{e'}^*(\mathbf{r}_1) \phi_{e_1}(\mathbf{r}_1)` is the **ejection transition density** describing excitation of the spectator carrier into the continuum.


.. rubric:: From ``docs/part8_auger/index.rst:111-123``

Spin Summation & Multiplicity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a closed-shell spatial orbital representation where conduction band edge electrons have paired spins, the two electrons can occupy parallel or antiparallel spin states:
* Parallel spins: only the antisymmetric channel :math:`|V_{\mathrm{dir}} - V_{\mathrm{exch}}|^2` is spin-allowed.
* Antiparallel spins: exchange vanishes due to orthogonal spin factors, leaving :math:`|V_{\mathrm{dir}}|^2`.

Summing over all spin-allowed final configurations yields the transition probability weight:

.. math::

   \sum_{\sigma_f} |M_{\mathrm{Auger}}|^2 = |V_{\mathrm{dir}} - V_{\mathrm{exch}}|^2 + |V_{\mathrm{dir}}|^2


.. rubric:: From ``docs/part8_auger/index.rst:124-138``

Universal Multiexciton Statistical Scaling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a neutral biexciton containing :math:`2e + 2h`:
* Either of the 2 conduction electrons can recombine with either of the 2 valence holes, giving :math:`2 \times 2 = 4` independent electron-hole recombination channels.
* Energy can be transferred to either the remaining spectator electron (:math:`eeh`) or the remaining spectator hole (:math:`hhe`).

The relation tested against single-dot and ensemble lifetimes is the superposition principle (Park, Lim, Klimov et al., ACS Nano 2017; Hou, Peng et al., Nat. Commun. 2019):

.. math::

   k_{XX} = 2 k_{X^-} + 2 k_{X^+} \qquad \Longleftrightarrow \qquad \frac{1}{\tau_{XX}} = 2\left(\frac{1}{\tau_{X^-}} + \frac{1}{\tau_{X^+}}\right)

Here :math:`k_{X^-}` and :math:`k_{X^+}` are the physical trion rates. For a twofold 1S shell each of those rates is twice the elementary three-carrier pathway, because either of the two identical carriers can recombine. When the two trion rates are equal this is :math:`\tau_{\mathrm{trion}} = 4 \tau_{XX}`. It is not an identity if one channel dominates.


.. rubric:: From ``docs/part8_auger/index.rst:139-148``

Schematic Diagrams of the Auger Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The energy-level configurations, particle movements, and coupling formulas for the :math:`eeh`, :math:`hhe`, and biexciton :math:`XX` processes are summarized schematically below:

**1. Negative Trion / Biexciton eeh Channel (Electron Ejected to Continuum)**

.. code-block:: text

   =========================================================================================================

                             NEGATIVE TRION / BIEXCITON: eeh (Electron Ejected)
   =========================================================================================================

          Energy
            ▲
            │                                                              e' (Hot Continuum Electron)
            │                                                            ┌─────┐
            │                                                            │  ●  │ ε_e' = ε_e1 + ΔE_recomb
            │                                                            └─────┘
            │                                                               ▲
            │                                                               │  ΔE_eject = ε_e' - ε_e1
            │                                                               │  (Absorbs Coulomb energy)
            │                                                               │
     CB     ┼─── CBM ───   ┌─────┐             ┌─────┐                  ┌─────┐
            │              │ e_1 │(Spectator)  │ e_2 │(Recombines)      │     │
            │              │  ●  │             │  ●  │                  │     │
            │              └─────┘             └─────┘                  └─────┘
            │                 │                   │
            │                 │                   │ ΔE_recomb = ε_e2 - ε_h ≈ E_g
            │                 │                   ▼
     VB     ┼─── VBM ───   ┌─────┐             ┌─────┐                  ┌─────┐
            │              │     │             │  ○  │(Hole h)          │     │ (Recombined: 0 holes)
            │              │     │             │     │                  │     │
            │              └─────┘             └─────┘                  └─────┘
            │
            │             INITIAL STATE (e_1, e_2, h)                  FINAL STATE (e')
            └────────────────────────────────────────────────────────────────────────────────────────►
                                                         Reaction Coordinate

     Formula Mapping:
       • Recombination energy:    ΔE_recomb = ε_e2 - ε_h ≈ E_g
       • Recombination charge:    q_A^{e2 h} = Σ_{λ∈A, σ} C_{λ e2} S_{λσ} C_{σ h}
       • Ejection charge:         q_B^{e' e1} = Σ_{μ∈B, ν} C_{μ e'} S_{μν} C_{ν e1}
       • Direct Coulomb integral: V_dir  = (q^{e2 h})^T  W^Resta  q^{e' e1}
       • Exchange integral:       V_exch = (q^{e1 h})^T  W^Resta  q^{e' e2}
       • Net matrix element:      |M_if|² = |V_dir - V_exch|² + |V_dir|²
       • Recombination rate:      Γ_eeh = (2π / ħ) Σ_e' |M_if|² ρ(ε_e' - ε_e1 - ΔE_recomb)

**2. Positive Trion / Biexciton hhe Channel (Hole Pushed Deep into Valence)**

.. code-block:: text

   =========================================================================================================

                             POSITIVE TRION / BIEXCITON: hhe (Hole Ejected)
   =========================================================================================================

          Energy
            ▲
            │
     CB     ┼─── CBM ───   ┌─────┐             ┌─────┐                  ┌─────┐
            │              │  e  │(Recombines) │     │                  │     │ (Recombined: 0 electrons)
            │              │  ●  │             │     │                  │     │
            │              └─────┘             └─────┘                  └─────┘
            │                 │
            │                 │ ΔE_recomb = ε_e - ε_h2 ≈ E_g
            │                 ▼
     VB     ┼─── VBM ───   ┌─────┐             ┌─────┐                  ┌─────┐
            │              │  ○  │(Hole h_2)   │  ○  │(Spectator h_1)   │     │
            │              │     │             │     │                  │     │
            │              └─────┘             └─────┘                  └─────┘
            │                                     │
            │                                     │ ΔE_eject = ε_h1 - ε_h'
            │                                     │ (Pushed deep into valence)
            │                                     ▼
            │                                                           ┌─────┐
            │                                                           │  ○  │ h' (Deep Valence Hole)
            │                                                           └─────┘ ε_h' = ε_h1 - ΔE_recomb
            │
            │             INITIAL STATE (e, h_1, h_2)                  FINAL STATE (h')
            └────────────────────────────────────────────────────────────────────────────────────────►
                                                         Reaction Coordinate

     Formula Mapping:
       • Recombination energy:    ΔE_recomb = ε_e - ε_h2 ≈ E_g
       • Recombination charge:    q_A^{e h2} = Σ_{λ∈A, σ} C_{λ e} S_{λσ} C_{σ h2}
       • Ejection charge:         q_B^{h1 h'} = Σ_{μ∈B, ν} C_{μ h1} S_{μν} C_{ν h'}
       • Direct Coulomb integral: V_dir  = (q^{e h2})^T  W^Resta  q^{h1 h'}
       • Exchange integral:       V_exch = (q^{e h1})^T  W^Resta  q^{h2 h'}
       • Net matrix element:      |M_if|² = |V_dir - V_exch|² + |V_dir|²
       • Recombination rate:      Γ_hhe = (2π / ħ) Σ_h' |M_if|² ρ(ε_h1 - ε_h' - ΔE_recomb)

**3. Neutral Biexciton Annihilation (XX) & ECSH Dynamical Relaxation Cascade**

.. code-block:: text

   =========================================================================================================

         NEUTRAL BIEXCITON (XX) ANNIHILATION & ENERGY-CONSERVING SURFACE HOPPING (ECSH) CASCADE
   =========================================================================================================

     Step 1: Initial Biexciton State (2e + 2h)
             Energy: E_XX ≈ 2 * E_g
             Conduction: [ e_1 , e_2 ]      Valence: [ h_1 , h_2 ]
                    │
                    │  Auger Annihilation Hop: ΔN_orb = 2 (Two-Particle Coulomb Operator)
                    │  Statistical Scaling: Γ_XX = 4 * Γ_eeh + 4 * Γ_hhe
                    │  ECSH Rule: Energy conserved internally within electrons (|E_XX - E_X*| <= k_B*T)
                    │  NO phonon loss! NO Boltzmann damping! (Zero nuclear velocity rescaling)
                    ▼
     Step 2: Hot Single Exciton State (X*)
             One (e, h) pair annihilated; excess energy absorbed by spectator carrier:
             • Either Hot Electron State: (e', h_2) with ε_e' ≈ CBM + E_g
             • Or Hot Hole State:         (e_1, h') with ε_h' ≈ VBM - E_g
                    │
                    │  Non-Adiabatic Electron-Phonon Cascades: ΔN_orb = 1 (Single-Particle NACs)
                    │  Coupling: d_ij = <φ_i | ∂/∂t | φ_j> along AIMD trajectory
                    │  Dissipates electronic energy into lattice vibrations (heat/phonons)
                    │  Upward hops Boltzmann-scaled: P * exp(-ΔE / k_B*T)
                    ▼
     Step 3: Thermalized Band-Edge Exciton (1S)
             Cold electron at CBM (1S_e), Cold hole at VBM (1S_h)
             Energy: E_1S ≈ E_g  -->  Subsequent slow radiative emission (τ_rad ~ 10 - 50 ns)

---
