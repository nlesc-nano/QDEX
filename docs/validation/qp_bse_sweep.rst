Reproducing the results: the QP × BSE sweep
===========================================

Part of :doc:`/validation/index`.

``benchmarks/qp_bse_sweep.py`` runs every consistent combination of QP model, kernel and integral
representation for one dot. Copy it into a folder that holds ``config.yaml`` and its input files:

.. code-block:: bash

   cp benchmarks/qp_bse_sweep.py tests/CdSe/2.0nm/
   cd tests/CdSe/2.0nm
   python qp_bse_sweep.py --list                  # 70 cases
   python qp_bse_sweep.py --jobs 2 --nthreads 4   # run (resumable)
   python qp_bse_sweep.py --collect               # sweep/summary.md, summary.csv

For dots with thousands of basis functions use the targeted profile (fixed 25 × 25 active space,
Davidson, 36 cases, results in ``sweep_large/``):

.. code-block:: bash

   python qp_bse_sweep.py --profile large --nthreads 16
   python qp_bse_sweep.py --profile large --no-xs --no-qsgw     # cheapest subset

Each case folder holds its complete ``config.yaml`` in the current layout, so any case can be rerun
with ``qdex --config config.yaml``.

Models
------

.. list-table::
   :header-rows: 1

   * - Case name
     - ``quasiparticles.model``
     - ``excitations``
   * - ``sbse-resta``, ``sbse-dim``
     - ``bulk``
     - ``mode: sbse``, ``kernel: resta`` / ``dim``
   * - ``brus``
     - ``brus``
     - ``kernel: resta``
   * - ``sgw-resta``, ``evgw-resta``, ``qsgw-resta``
     - same
     - ``kernel: qp`` (the model's W)
   * - ``sgw-dim``, ``evgw-dim``, ``qsgw-dim``
     - same
     - ``kernel: qp``

Groups
------

.. list-table::
   :header-rows: 1

   * - Group
     - Question
   * - A
     - every model, {mnok, xs} × {vacuum, solvent}
   * - B
     - Z = 1 and Z = 0.8 against the plasmon-pole Z
   * - C
     - classical self-energy vs ΔCOHSEX; Born vs sphere environment term; rigid vs orbital levels
   * - E
     - SOC in the solvent for ``sbse-resta``, ``sgw-resta``, ``sgw-dim``, ``qsgw-dim``
   * - F
     - active space 50 × 50 and 100 × 100 (Davidson)
   * - G
     - Löwdin instead of Mulliken charges
   * - H
     - triplets; diagonal solvers

The tables of :doc:`cdse_experiment` correspond to groups A and E with mnok. Compare the bright SOC
state in the solvent with experiment; ``benchmarks/compare_models.py --exp-ref aubert-hens-2022-zb``
evaluates the sizing curves at the cluster's core diameter.
