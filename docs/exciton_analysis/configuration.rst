Configuration
=============

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

.. rubric:: From ``docs/part5_exciton_analysis/index.rst:223-225``

6. CLI Flags & YAML Configuration Reference
-------------------------------------------


.. rubric:: From ``docs/part5_exciton_analysis/index.rst:226-260``

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--nto``
     - ``False``
     - Run Natural Transition Orbital (NTO) analysis after solving BSE.
   * - ``--nto-states <list>``
     - ``1 2 3``
     - Specific exciton states to analyze with NTOs (1-indexed, e.g. ``--nto-states 1 5 10``).
   * - ``--nto-top <int>``
     - ``3``
     - Number of dominant NTO pairs to report per state.
   * - ``--nto-csv``
     - ``False``
     - Export detailed NTO descriptors and weights to ``nto_results.csv``.
   * - ``--bse_states <list>``
     - ``1 2 3``
     - Specific exciton roots to export as 3D volumetric ``.cube`` files.
   * - ``--nbse <int>``
     - ``3``
     - Number of lowest exciton states to export as ``.cube`` files if ``--bse_states`` is omitted.
   * - ``--plot``
     - ``False``
     - Generate publication-ready figures and interactive Plotly HTML dashboards.
   * - ``--show``
     - ``False``
     - Display interactive plots in the web browser upon calculation completion.


.. rubric:: From ``docs/part5_exciton_analysis/index.rst:261-276``

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   analysis:
     nto: true
     nto_states: [1, 2, 3]
     nto_top: 3
     nto_csv: true
     plot: true

   cube:
     export: true
     bse_states: [1, 2]
     spacing_ang: 0.5
