Qp analysis
===========

Part of :doc:`/electronic_structure/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.pdos_coop``
* Callable: ``qdex.pdos_coop.compute_pdos_and_coop``
* CLI: ``--charge_type, --run_fuzzy``
* YAML: ``physics.charge_type, fuzzy.run``

.. code-block:: python

   compute_pdos_and_coop(C, S, eps_eV, shells, pdos_atoms, coop_pairs, ewin, sigma=0.03, is_soc=False, prefix='sf', pops=None, population_bars=None, device='numpy')


14. Combining Electronic Structure Analysis with QP Shifts
----------------------------------------------------------

``QDEX`` seamlessly integrates QP corrections into ground-state electronic structure analyses:

* **QP-Corrected PDOS**: Convolves the projected density of states on the corrected quasiparticle energy axis :math:`\varepsilon^{\mathrm{QP}}`, opening the gap to experimental values while preserving Mulliken orbital weights.
* **QP-Corrected Fuzzy Bands**: Unfolds nanocrystal orbitals onto bulk :math:`k`-paths using the quasiparticle dispersion:

  .. math::

     P_{\mathbf{k}}(\varepsilon) = \sum_m |\langle e^{i \mathbf{k} \cdot \mathbf{r}} | \phi_m \rangle|^2 \, \delta(\varepsilon - \varepsilon_m^{\mathrm{QP}}).

Controlled via ``--dashboard_energy_mode {dft, qp, both}``.

