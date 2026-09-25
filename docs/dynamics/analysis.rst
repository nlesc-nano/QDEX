Analysis
========

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.precompute``
* Callable: ``qdex.namd.precompute.precompute_namd_data``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   precompute_namd_data(config)


9. In-Depth Analysis of NAMD Simulations
----------------------------------------

``QDEX`` includes a dedicated analysis module (``qdex.namd.analysis``) that automatically processes precomputed and dynamic trajectory data.


1. Carrier Cooling Curves, Lifetimes, and Band Edge Arrival Times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Carrier relaxation is monitored by tracking the time-dependent excess energy of electrons (:math:`\Delta E_e`), holes (:math:`\Delta E_h`), and total exciton (:math:`\Delta E_{\mathrm{exc}}`) above their respective band edges:

.. math::

   \Delta E_e(t) = \sum_{a \in \mathrm{virt}} P_a(t) \left( \varepsilon_a(t) - \varepsilon_{\mathrm{LUMO}}(t) \right)

.. math::

   \Delta E_h(t) = \sum_{i \in \mathrm{occ}} P_i(t) \left( \varepsilon_{\mathrm{HOMO}}(t) - \varepsilon_i(t) \right)

Cooling Rates and Exponential Lifetimes
"""""""""""""""""""""""""""""""""""""""

The cooling lifetimes :math:`\tau_{\mathrm{exc}}`, :math:`\tau_e`, and :math:`\tau_h` are extracted by linear regression on the logarithmic excess energy decay:

.. math::

   \Delta E(t) \approx \Delta E(0) \, \exp\left( -\frac{t}{\tau_{\mathrm{cooling}}} \right)

The corresponding carrier cooling rate is:

.. math::

   k_{\mathrm{cool}} = \frac{1}{\tau_{\mathrm{cooling}}} \quad [\text{ps}^{-1}]

Time to Reach the Band Edge: Analytical Estimates vs. Actual Trajectory
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In experiments and device modeling, researchers require not only the initial relaxation slope :math:`\tau`, but also the **total time required for hot carriers to fully reach the band edge**:

1. **Analytical Estimates from Exponential Fit**:
   * **95% Excess Energy Dissipated**:

     .. math::

        t_{\mathrm{est}}^{95\%} = -\ln(0.05) \, \tau \approx 3.0 \, \tau

   * **99% Excess Energy Dissipated (Complete Thermalization)**:

     .. math::

        t_{\mathrm{est}}^{99\%} = -\ln(0.01) \, \tau \approx 4.6 \, \tau

2. **Actual Numerical Arrival Times from Simulation Trajectory**:
   * :math:`t_{\mathrm{act}}^{95\%}`: The first simulation timestamp :math:`t` where :math:`\Delta E(t) \le 0.05 \, \Delta E(0)`.
   * :math:`t_{\mathrm{act}}^{99\%}`: The first simulation timestamp :math:`t` where :math:`\Delta E(t) \le 0.01 \, \Delta E(0)`.
   * :math:`t_{\mathrm{act}}^{\mathrm{therm}}`: The first timestamp where excess energy drops below the thermal energy of the lattice bath:

     .. math::

        \Delta E(t) \le k_B T \approx 25.8\text{ meV at } 300\text{ K}

**Physical Insight**: While an ideal exponential decay satisfies :math:`t_{\mathrm{act}} \approx t_{\mathrm{est}}`, realistic atomistic trajectories often exhibit non-exponential behaviors—such as an initial **phonon bottleneck** across discrete sub-bands or delayed cascades through intermediate surface states. Comparing :math:`t_{\mathrm{est}}` with :math:`t_{\mathrm{act}}` immediately diagnostics whether carrier cooling proceeds smoothly or is delayed by bottlenecks.


2. State-Resolved Population Kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Transient populations :math:`P_I(t)` are exported to ``carrier_cooling_populations.csv``, displaying the sequential decay of initial hot excitons into intermediate states and finally into the emitting :math:`1S` state.


3. NAC vs. Energy Gap Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To verify whether non-adiabatic transitions obey the energy-gap law, ``QDEX`` samples pairs of states across trajectory frames and plots non-adiabatic coupling magnitudes :math:`|d_{IJ}|` against energy differences :math:`|E_J - E_I|`. This distinguishes smooth exponential decay from resonant vibronic enhancements.


4. 6-Panel Publication Figures & Dashboards
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Executing the analysis workflow generates a comprehensive 6-panel summary figure (``namd_analysis_6panel.png``) and an interactive Plotly HTML dashboard:

* **Panel A**: Carrier cooling curves (:math:`\Delta E_e(t)` vs. :math:`\Delta E_h(t)` with fitted lifetimes).
* **Panel B**: Time-dependent populations of frontier exciton states.
* **Panel C**: Band-gap thermal fluctuation trajectory :math:`E_g(t)`.
* **Panel D**: Energy gap autocorrelation function :math:`C(t)` and cumulant dephasing decay :math:`D(t)`.
* **Panel E**: Phonon Spectral Density :math:`J(\omega)` in :math:`\text{cm}^{-1}`.
* **Panel F**: Non-adiabatic coupling distribution :math:`|d_{IJ}|` vs. :math:`\Delta E_{IJ}`.

