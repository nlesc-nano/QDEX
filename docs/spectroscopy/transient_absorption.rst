Transient absorption
====================

Part of :doc:`/spectroscopy/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.namd.transient_absorption``
* Callable: ``qdex.namd.transient_absorption.compute_transient_absorption``
* CLI: ``--sigma``
* YAML: ``namd.transient_absorption``

.. code-block:: python

   compute_transient_absorption(times_fs: np.ndarray, populations: np.ndarray, E_pairs: np.ndarray, f_pairs: np.ndarray, i_pairs: np.ndarray, a_pairs: np.ndarray, sigma_ev: float=0.03, e_range: Optional[Tuple[float, float]]=None, n_e_points: int=300, include_se: bool=False, all_energies: Optional[np.ndarray]=None, state_degeneracy: float=2.0, include_esa: bool=False, f_elec_esa: Optional[np.ndarray]=None, e_elec_esa: Optional[np.ndarray]=None)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part6_namd/index.rst:1158-1168``

10. Ultrafast Pump-Probe Transient Absorption (TA) Spectroscopy
---------------------------------------------------------------

In ultrafast optical experiments on semiconductor nanocrystals, **pump-probe transient absorption (TA)** spectroscopy tracks the non-equilibrium evolution of photoexcited carriers. A high-energy pump pulse creates an initial non-thermal carrier distribution, while a broadband, time-delayed probe pulse monitors the differential absorbance:

.. math::

   \Delta A(E, t) = A_{\mathrm{pump-on}}(E, t) - A_{\mathrm{pump-off}}(E)

as a function of probe photon energy :math:`E` and pump–probe delay time :math:`t`.


.. rubric:: From ``docs/part6_namd/index.rst:1169-1182``

Physical Mechanisms in Nanocrystal Transient Absorption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transient absorption signal :math:`\Delta A(E, t)` comprises three distinct physical mechanisms:

1. **Ground-State Bleach (GSB) & State-Filling (Pauli Blocking)** (:math:`\Delta A < 0`):
   When electrons occupy virtual conduction states :math:`a` and holes occupy valence states :math:`i`, optical transitions between these occupied states are blocked by the Pauli exclusion principle. As carriers cascade down the non-adiabatic ladder toward the band edge, absorption into the lowest exciton state (:math:`1S`) is progressively quenched, forming the prominent negative **1S Bleach**.

2. **Stimulated Emission (SE)** (:math:`\Delta A < 0`):
   Probe photons passing through the sample stimulate coherent radiative transitions from populated excited levels back to the ground state. Because stimulated emission adds photons to the transmitted probe beam, it appears as a negative differential absorption feature at the same transition frequencies as the ground-state bleach.

3. **Excited-State Absorption (ESA)** (:math:`\Delta A > 0`):
   A populated excited carrier can absorb a second probe photon to transition into higher-lying continuum bands or multiexciton states (:math:`S_1 \to S_{XX}`), producing positive absorption features.


.. rubric:: From ``docs/part6_namd/index.rst:1183-1207``

Microscopic Formulation & State-Filling Factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In ``QDEX``, the instantaneous exciton pair populations :math:`P_{ia}(t)` (where :math:`i` denotes occupied valence orbitals and :math:`a` denotes virtual conduction orbitals) directly govern state filling.

At delay time :math:`t_k`, the single-particle fractional electron occupation in conduction orbital :math:`a` and hole occupation in valence orbital :math:`i` are obtained by tracing out the complementary carrier:

.. math::

   n_a(t_k) = \sum_{i} P_{ia}(t_k), \qquad p_i(t_k) = \sum_{a} P_{ia}(t_k)

In the unexcited ground state, each transition :math:`(i \to a)` has an intrinsic oscillator strength :math:`f_{ia}^{(0)}`. For incoherent populations the differential strength is the degeneracy-normalized state-filling signal used for quantum-dot transient absorption (Grimaldi et al., Nano Lett. 2019). A closed-shell spatial orbital holds two carriers, so :math:`g = 2`; a spinor has :math:`g = 1`. One carrier then blocks a fraction :math:`1/g` of that shell. Net optical gain requires :math:`n_a/g_a + p_i/g_i > 1`.

.. math::

   \Delta f_{ia}(t_k) = - f_{ia}^{(0)} \left( \frac{n_a(t_k)}{g_a} + \frac{p_i(t_k)}{g_i} \right)

An extra :math:`-f P_{ia}` term is available only through ``include_se`` and is not the default. Adding it on top of state filling counts stimulated emission twice.

The continuous differential absorption spectrum :math:`\Delta A(E, t_k)` is then evaluated by convolving with a Gaussian probe spectral broadening :math:`\sigma`:

.. math::

   \Delta A(E, t_k) = \sum_{ia} \Delta f_{ia}(t_k) \frac{1}{\sqrt{2\pi}\sigma} \exp\left( -\frac{(E - E_{ia}(t_k))^2}{2\sigma^2} \right)


.. rubric:: From ``docs/part6_namd/index.rst:1208-1228``

1S Bleach Kinetic Profiling & Carrier Cooling Rates (:math:`k_C`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The hallmark signature of carrier cooling in transient absorption experiments is the **rise profile of the 1S bleach**. 

At :math:`t = 0`, photoexcitation creates carriers at high excess energy (e.g. :math:`2 \times E_g`). The initial bleach :math:`\Delta A(E, 0)` is broad and located in the high-energy window. As non-adiabatic electron-phonon scattering drives carriers toward the band edge, population accumulates in the lowest exciton state (:math:`1S`, :math:`E_{1S} \approx E_g`), causing the negative 1S bleach signal :math:`-\Delta A_{1S}(t)` to rise from zero to its plateau value.

``QDEX`` extracts the band-edge 1S bleach:

.. math::

   S_{1S}(t) = - \Delta A(E_{1S}, t)

and performs non-linear least-squares fitting to an exponential rise model:

.. math::

   S_{1S}(t) = A_0 \left( 1 - \exp\left( -\frac{t}{\tau_C} \right) \right)

This directly yields the **hot-carrier cooling time** :math:`\tau_C` (in fs/ps) and **cooling rate** :math:`k_C = 1 / \tau_C` (in :math:`\text{ps}^{-1}`), providing direct theoretical counterparts to the experimental 1S bleach rise traces reported in transient absorption spectroscopy.


.. rubric:: From ``docs/part6_namd/index.rst:1229-1238``

Schematic Diagrams of Transient Absorption Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optical mechanisms (GSB, SE, ESA) and dynamical carrier relaxation cascade underlying transient absorption are illustrated schematically below:

**1. Three Optical Mechanisms in Nanocrystal Transient Absorption**

.. code-block:: text

   =========================================================================================================

             THREE OPTICAL MECHANISMS IN TRANSIENT ABSORPTION SPECTROSCOPY (ΔA)
   =========================================================================================================

     (a) Ground-State Bleach (GSB)            (b) Stimulated Emission (SE)             (c) Excited-State Absorption (ESA)
         ΔA < 0  (Pauli Blocking)                 ΔA < 0  (Probe Photons Added)            ΔA > 0  (Induced Absorption)

         Energy                                   Energy                                   Energy
           ▲                                        ▲                                        ▲
           │                                        │                                        │  Conduction Continuum
           │                                        │                                        │  ┌─────┐ a'
           │                                        │                                        │  │     │
           │                                        │                                        │  └─────┘
           │                                        │                                        │     ▲
           │                                        │                                        │     │  Probe Absorbed:
           │                                        │                                        │     │  hν_probe = E_a' - E_a
    CB  ┼──┼── CBM                               CB ┼──┼── CBM                            CB ┼──┼──┤  Δf_aa'^ESA > 0
           │  ┌─────┐ a                             │  ┌─────┐ a                             │  ┌──┴──┐ a (Populated)
           │  │  ●  │ n_a (Blocked!)                │  │  ●  │ P_ia (Populated pair)         │  │  ●  │ n_a
           │  └─────┘                               │  └─────┘                               │  └─────┘
           │     ▲                                  │     │
           │     │  Probe Blocked!                  │     │  Probe Stimulates:               │
           │     │  hν_probe ≈ E_ia                 │     │  hν_probe + hν_em (2 photons!)   │
           │     │  Δf_ia^GSB < 0                   │     ▼  Δf_ia^SE < 0                    │
    VB  ┼──┼── VBM                               VB ┼──┼── VBM                            VB ┼──┼── VBM
           │  ┌─────┐ i                             │  ┌─────┐ i                             │  ┌─────┐ i
           │  │  ○  │ p_i (Blocked!)                │  │  ○  │                               │  │  ●  │
           │  └─────┘                               │  └─────┘                               │  └─────┘
           └────────────────────────►               └────────────────────────►               └────────────────────────►

      Formula Mapping:                         Formula Mapping:                         Formula Mapping:
        • n_a = Σ_i P_ia                         • Default signal is state filling        • Optional, if orbital dipoles exist
        • p_i = Σ_a P_ia                         • already. include_se adds -f P_ia      • Δf_aa' = +f_aa' n_a / g
        • Δf = -f (n_a/g + p_i/g)                • and double-counts emission.           • ΔA_ESA > 0
        • g = 2 spatial, g = 1 spinor            • Net gain only if n/g + p/g > 1
        • One spatial exciton: Δf = -f

**2. Hot-Carrier Relaxation Cascade & 1S Bleach Kinetic Rise Profile**

.. code-block:: text

   =========================================================================================================

         HOT-CARRIER RELAXATION CASCADE & 1S BLEACH KINETIC RISE PROFILE
   =========================================================================================================

     1. Hot Excitation & Vibronic Cooling Cascade           2. Time-Resolved 1S Bleach Profile S_1S(t)
        Energy                                                 Bleach: S_1S(t) = -ΔA(E_1S, t)
          ▲                                                      ▲
          │  Pump: hν_pump ≈ 2 * E_g                             │
          │    ┌─────┐ e* (Hot Electron)                         │                      Plateau A_0
          │    │  ●  │                                           │               . - - - - - - - - - - - - -
          │    └─────┘                                           │           . '
          │       │                                              │        . '
          │       │ Non-adiabatic cooling                        │      . '   S_1S(t) = A_0 (1 - e^{-t/τ_C})
          │       │ d_ab = ⟨a|∂/∂t|b⟩                            │     /
          │       ▼ (Phonon emission cascade)                    │    /
     CB ──┼─── ┌─────┐ 1S_e (CBM Band Edge)                      │   /
          │    │  ●  │ State filling n_1Se(t) rises!             │  /
          │    └─────┘                                           │ /
          │                                                      │/
          │       E_1S ≈ E_g (Probe monitors 1S)                 └────────────────────────────────────────►
          │                                                      0          τ_C               t (Delay)
     VB ──┼─── ┌─────┐ 1S_h (VBM Band Edge)
          │    │  ○  │ State filling p_1Sh(t) rises!          Kinetic Parameters:
          │    └─────┘                                          • Hot carrier cooling time: τ_C (fs or ps)
          │       ▲ (Phonon emission cascade)                   • Carrier cooling rate:     k_C = 1 / τ_C (ps⁻¹)
          │       │ Non-adiabatic cooling                       • Band-edge 1S energy:      E_1S = ε_1Se - ε_1Sh - E_b
          │       │ d_ij = ⟨i|∂/∂t|j⟩                           • Bleach signal:            -ΔA(E_1S, t)
          │       │
          │    ┌─────┐
          │    │  ○  │ h* (Hot Hole)
          │    └─────┘


.. rubric:: From ``docs/part6_namd/index.rst:1311-1350``

How QDEX Data Are Used to Compute Every Formula Term
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is the exact step-by-step mapping of how QDEX data structures evaluate every term in the transient absorption equations:

.. list-table::
   :widths: 25 35 40
   :header-rows: 1

   * - Mathematical Term
     - QDEX Data Structure / Source
     - Computational Implementation
   * - **Transition Dipoles** :math:`\boldsymbol{\mu}_{ia}`
     - ``compute_dipole_ao(shells)`` & MO transformation
     - Dipole integrals :math:`\langle \mu | e\mathbf{r} | \nu \rangle` via Libint2, transformed to MO basis: :math:`\boldsymbol{\mu}_{ia} = C_{\mathrm{occ}}^T \boldsymbol{\mu}_{\mathrm{ao}} C_{\mathrm{virt}}`.
   * - **Ground-State Strengths** :math:`f_{ia}^{(0)}`
     - ``f_pairs`` in ``frame_00000.npz``
     - Unperturbed oscillator strengths: :math:`f_{ia}^{(0)} = \frac{2}{3} \Delta E_{ia} |\boldsymbol{\mu}_{ia}|^2`.
   * - **Instantaneous Energies** :math:`E_{ia}(t)`
     - ``E_curr`` from ``step_*.npz`` or ``all_energies``
     - Time-dependent diagonal BSE / QP transition energies tracking nuclear MD motion.
   * - **Time-Dependent Populations** :math:`P_{ia}(t)`
     - ``populations[k]`` from ``run_namd_dynamics``
     - Dynamic state populations propagated via Pauli Master Equation or CPA-FSSH.
   * - **Virtual Occupation** :math:`n_a(t)`
     - ``np.bincount(a_pairs, weights=P_k)``
     - Traced electron occupation across all active virtual orbitals :math:`a`.
   * - **Occupied Occupation** :math:`p_i(t)`
     - ``np.bincount(i_pairs, weights=P_k)``
     - Traced hole occupation across all active valence orbitals :math:`i`.
   * - **Differential Strength** :math:`\Delta f_{ia}(t)`
     - ``delta_f = - f_base * (n_a[a_arr] + p_i[i_arr] + P_k)``
     - Evaluated at every time step :math:`t_k` with vector broadcasting.
   * - **2D Differential Absorbance** :math:`\Delta A(E, t)`
     - ``np.dot(delta_f, gauss_profiles)``
     - Vectorized Gaussian spectral convolution over probe energy grid :math:`E`.
   * - **1S Rise Fitting** :math:`\tau_C, k_C`
     - ``fit_bleach_rise_kinetics(times, delta_A_1s)``
     - Scipy curve-fit of :math:`-\Delta A_{1S}(t)` to exponential rise :math:`A_0 (1 - e^{-t/\tau_C})`.


.. rubric:: From ``docs/part6_namd/index.rst:1351-1360``

Publication-Quality Visualizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Executing transient absorption generates a 3-panel publication dashboard:
* **Panel (a): 2D False-Color TA Map** :math:`\Delta A(E, t)`: Probe energy on x-axis, delay time on y-axis, using a diverging `RdBu_r` colormap (blue = negative bleach, red = positive ESA).
* **Panel (b): 1S Bleach Kinetic Rise Profile**: Tracking the negative 1S bleach with exponential rise fit and annotated cooling time :math:`\tau_C` and rate :math:`k_C`.
* **Panel (c): Differential Absorption Spectra** :math:`\Delta A(E)`: Spectral slices at selected delay times (e.g. :math:`t = 0, 100, 250, 500, 1000\text{ fs}`) showing the spectral shift from hot state filling to the sharp band-edge bleach.

---
