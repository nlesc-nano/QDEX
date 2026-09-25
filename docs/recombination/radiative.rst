Radiative
=========

Part of :doc:`/recombination/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.compute_radiative_rates``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   compute_radiative_rates(E_ev, f_osc, refractive_index=2.0)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part6_namd/index.rst:891-918``

2. Einstein Radiative Rate: Single-Frame vs. NAMD Trajectory Averaging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spontaneous photon emission into the vacuum radiation field inside a dielectric medium of refractive index :math:`n_{\mathrm{r}}` is given by the Einstein A coefficient:

.. math::

   k_{I \to 0}^{\mathrm{rad}} = \left[ \frac{2 e^2}{4\pi \epsilon_0 m_e c^3 \hbar^2} \right] n_{\mathrm{r}} \, E_I^2 \, f_I = C_{\mathrm{rad}} \, n_{\mathrm{r}} \, E_I^2 \, f_I

where :math:`C_{\mathrm{rad}} = 4.3391988 \times 10^7\text{ s}^{-1}\text{ eV}^{-2}` (:math:`4.3391988 \times 10^{-8}\text{ fs}^{-1}\text{ eV}^{-2}`) and :math:`n_{\mathrm{r}}` is loaded from ``REFRACTIVE_INDEX_DICT`` in ``qdex.hardness`` (e.g. :math:`n_{\mathrm{r}} = 2.19` for :math:`\text{CsPbBr}_3`, :math:`n_{\mathrm{r}} = 3.30` for :math:`\text{GaAs}`).

Single-Frame vs. Trajectory Ensemble Averaging
""""""""""""""""""""""""""""""""""""""""""""""

A critical question is whether :math:`f_I` and :math:`E_I` should be taken from a single static snapshot (frame 0) or averaged along the NAMD trajectory:

1. **Static / Single-Frame Rate** (:math:`k_{\mathrm{rad}}(t=0)`):
   Evaluates :math:`E_I(0)` and :math:`f_I(0)` at the relaxed ground-state equilibrium geometry.
2. **Thermalized Band-Edge Rate at Frame 0**:
   Because the fine-structure splitting between band-edge exciton states (:math:`1-10\text{ meV}`) is much smaller than thermal energy (:math:`k_B T \approx 25.8\text{ meV}` at :math:`300\text{ K}`), carriers rapidly reach thermal equilibrium among low-lying states before radiating:

   .. math::

      k_{\mathrm{rad}}^{\mathrm{therm}}(t=0) = \frac{\sum_I k_{\mathrm{rad}, I}(0) \, \exp\left( -\frac{E_I(0) - E_0(0)}{k_B T} \right)}{\sum_I \exp\left( -\frac{E_I(0) - E_0(0)}{k_B T} \right)}

3. **What the dynamics actually use**:
   Dipoles are computed on frame 0. The photoluminescence yield uses the frame-0 thermal average above, with the lowest-exciton lifetime taken from :math:`\arg\min E_I`, not from pair index 0. A trajectory average :math:`\langle k_{\mathrm{rad}}\rangle_{\mathrm{MD}}` would require the dipole on every frame. That average is not computed. A dark band-edge exciton at the first geometry therefore stays dark for the reported yield. Herzberg–Teller intensity borrowing along the trajectory is a real physical effect and is not yet in the rate.
