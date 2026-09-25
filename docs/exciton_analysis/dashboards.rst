Dashboards
==========

Part of :doc:`/exciton_analysis/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.nto``
* Callable: ``qdex.nto.analyze_nto_state``
* CLI: ``--nto, --nto-states``
* YAML: ``analysis.nto, analysis.nto_states``

.. code-block:: python

   analyze_nto_state(solver, vec, energy_ev, f_osc, state_index, coords, symbols, mu_ia=None, soc_U=None, top_n=3, context=None)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part5_exciton_analysis/index.rst:209-222``

5. Interactive Plotly 6-Panel Dashboard
---------------------------------------

When running with ``--plot``, ``QDEX`` exports a self-contained, interactive HTML dashboard (``exciton_analysis.html``) featuring six coordinated subplots:

1. **Participation Ratio (PR)** vs. Energy (identifying multiconfigurational states).
2. **True Exciton Size** (:math:`d_{eh}`) vs. Energy (characterizing Bohr radius scaling).
3. **Particle Spread** (:math:`\sigma_h, \sigma_e`) (comparing electron vs. hole delocalization).
4. **Charge-Transfer Distance** (:math:`d_{\mathrm{CT}}`) color-coded by CT ratio.
5. **Spatial Pearson Correlation** (:math:`R_{eh}`) (distinguishing bound from dissociated excitons).
6. **Simulated UV-Vis Absorption Spectrum** with oscillator strength stick spectra and Gaussian broadening.

---
