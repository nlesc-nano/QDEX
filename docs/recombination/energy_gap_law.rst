Energy gap law
==============

Part of :doc:`/recombination/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.compute_energy_gap_law_rate``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   compute_energy_gap_law_rate(E_gap_ev, E_LO_ev=0.018, S_hr=1.0, A_nr=10000000000000.0)


3. Non-Radiative Decay Across Large Gaps: Englman-Jortner Energy Gap Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Direct non-radiative recombination across a wide semiconductor band gap (:math:`E_g > 1.0\text{ eV}`) requires dissipating enormous electronic energy into the nuclear lattice. Because this involves dozens of vibrational quanta (:math:`p = E_g / \hbar \omega_{\mathrm{LO}} \sim 30 - 80`), multi-phonon perturbation theory (Englman & Jortner, 1970) yields the **Energy Gap Law**:

.. math::

   k_{\mathrm{nr}} = A_{\mathrm{nr}} \exp\left( -\gamma \frac{E_g}{\hbar \omega_{\mathrm{LO}}} \right)

where:

* :math:`\hbar \omega_{\mathrm{LO}}` is the dominant accepting optical phonon energy.
* :math:`\gamma` is the multi-phonon coupling parameter:

  .. math::

     \gamma = \ln\left( \frac{E_g}{S \, \hbar \omega_{\mathrm{LO}}} \right) - 1 = \ln\left( \frac{E_g}{\lambda} \right) - 1

  where :math:`S` is the dimensionless Huang-Rhys factor and :math:`\lambda = S \hbar \omega_{\mathrm{LO}}` is the nuclear reorganization energy.
* :math:`A_{\mathrm{nr}}` is the electronic prefactor:

  .. math::

     A_{\mathrm{nr}} = \frac{C_{\mathrm{el}}^2}{\hbar} \sqrt{\frac{2\pi}{\hbar \omega_{\mathrm{LO}} E_g}} \sim 10^{12} - 10^{13}\text{ s}^{-1}

  where :math:`C_{\mathrm{el}} \approx V_{\mathrm{el}}` is the non-adiabatic coupling matrix element between excited and ground states.
