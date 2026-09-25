Localization
============

Part of :doc:`/electronic_structure/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.pdos_coop``
* Callable: ``qdex.pdos_coop.compute_pdos_and_coop``
* CLI: ``--charge_type, --run_fuzzy``
* YAML: ``physics.charge_type, fuzzy.run``

.. code-block:: python

   compute_pdos_and_coop(C, S, eps_eV, shells, pdos_atoms, coop_pairs, ewin, sigma=0.03, is_soc=False, prefix='sf', pops=None, population_bars=None, device='numpy')


3. Inverse Participation Ratio (IPR)
------------------------------------

The **Inverse Participation Ratio (IPR)** provides a quantitative metric of wavefunction spatial localization. 


Mathematical Formulation
~~~~~~~~~~~~~~~~~~~~~~~~

In ``QDEX``, the orbital IPR for molecular orbital :math:`m` is defined in terms of its Mulliken atomic orbital weights :math:`P_{\mu, m}`:

.. math::

   \mathrm{IPR}_m = \sum_{\mu=1}^{N_{\mathrm{ao}}} \left( P_{\mu, m} \right)^2


Physical Interpretation
~~~~~~~~~~~~~~~~~~~~~~~

* **Delocalized Band States**: If an orbital is uniformly delocalized over :math:`N_{\mathrm{ao}}` atomic orbitals, each :math:`P_{\mu, m} \approx 1 / N_{\mathrm{ao}}`, resulting in:

  .. math::

     \mathrm{IPR}_m \approx N_{\mathrm{ao}} \left( \frac{1}{N_{\mathrm{ao}}} \right)^2 = \frac{1}{N_{\mathrm{ao}}} \ll 1

* **Localized Trap or Defect States**: If an orbital is localized entirely on a single dangling bond or surface atom, :math:`P_{\mu, m} \approx 1` for that site and 0 elsewhere, resulting in:

  .. math::

     \mathrm{IPR}_m \approx 1.0

Plotting IPR alongside the PDOS immediately identifies localized trap states within the fundamental gap or resonant near the band edges.

