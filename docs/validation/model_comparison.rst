Comparing QP and BSE options
============================

Part of :doc:`/validation/index`.

QDEX has three largely independent switches for the excited-state calculation. They combine into
many possible runs, and it is easy to lose track of which difference comes from which switch. This
page explains what each switch changes and describes a script that runs the combinations one factor
at a time and checks them.

.. figure:: /_static/figures/qp_hierarchy.svg
   :width: 100%
   :alt: qp hierarchy

   The ``qp_gap`` options and what each one corrects.

The three switches
------------------

.. list-table::
   :header-rows: 1
   :widths: 18 30 52

   * - switch
     - what it changes
     - options
   * - ``qp_gap``
     - the single-particle energies, i.e. the QP gap and level spacing
     - ``pbe`` (none). ``brus`` (EMA kinetic term). ``gw`` (anchor-scaled scissor). ``sgw-dim``,
       ``sgw-resta``, ``sgw`` (static charging models). ``evgw-*`` (gap iteration). ``qsgw-*``
       (orbital relaxation). ``env`` (bulk GW + dielectric-sphere polarization, shared with the BSE).
   * - ``kernel``
     - the screened interaction :math:`W` in the direct e–h term, which sets the exciton binding
     - ``resta`` (default: bulk :math:`\epsilon_\infty`, Thomas–Fermi-like profile). ``dim`` (Thole
       dipoles). ``sbse`` (monopole RPA). ``xs-resta`` (AO density-pair integrals + Resta). ``bse``
       (unscreened MNOK: an upper bound for the binding, not a physical choice).
   * - ``excitation_mode``
     - how the transition-space matrix is solved
     - ``independent_dft``, ``independent_qp`` (no kernel), ``diagonal_bse`` (kernels on the diagonal
       only), ``bse`` (full TDA BSE)

Two more settings matter:

* ``eps_out``: the solvent, which should be the optical :math:`n^2`.
* The active-space size ``nhomos``/``nlumos``: a convergence parameter, not a physical choice.

Which combination to use
------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - goal
     - recommended
   * - optical gap / absorption onset in a solvent
     - ``qp_gap: env`` + ``kernel: resta`` + ``excitation_mode: bse``, with ``eps_out`` = solvent
       :math:`n^2`. The optical gap is then solvent-consistent.
       Add ``env_anchor_residual: true`` to include the size dependence calibrated on the monomer GW
       anchor.
   * - QP gap / IP–EA trend with size
     - ``qp_gap: env`` or ``gw``. Report which, and give ``eps_out``.
   * - exciton binding energy
     - ``qp_gap: env`` (binding = QP − S\ :sub:`1` is environment-consistent), active space
       ≥ 100 × 100 with Davidson (``solver.full_diag: false``).
   * - fast screening of many structures
     - ``excitation_mode: diagonal_bse`` for peak positions near the edge (it recovers about 98 % of
       the S\ :sub:`1` shift for CdSe, but not the oscillator-strength distribution).
   * - absolute S\ :sub:`1` with the older scissor models (``gw``, ``sgw-*``)
     - only in the solvent they were calibrated for. Their optical gap moves with ``eps_out`` as much as
       the QP gap does (:doc:`environment_cancellation`).

The comparison script
---------------------

``benchmarks/compare_models.py`` runs the normal ``qdex`` CLI for each case, in its own directory,
starting from your ``config.yaml``. It then writes a report ``summary.md`` (plus ``summary.csv`` and
``summary.png``). The cases are grouped by the question they answer, so that only one factor changes
inside a group:

.. list-table::
   :header-rows: 1
   :widths: 10 40 50

   * - group
     - question
     - cases
   * - A
     - What QP gap does each ``qp_gap`` option give?
     - all QP options, vacuum, Resta kernel
   * - B
     - Is the optical gap insensitive to the solvent?
     - ``gw``, ``sgw-resta``, ``sgw-dim``, ``env`` at ``eps_out`` = 1 and the solvent value;
       ``env`` in a matched medium
   * - C
     - How much binding does each :math:`W` give?
     - ``resta``, ``dim``, ``xs-resta``, ``sbse``, bare MNOK at the same QP gap
   * - D
     - Do the solver modes obey their exact relations?
     - ``independent_dft``, ``independent_qp``, ``diagonal_bse``, triplet, Löwdin, no exchange
   * - E
     - Is S\ :sub:`1` converged in the active space?
     - 25 × 25, 50 × 50, 100 × 100 (Davidson)

Run it with:

.. code-block:: bash

   # quick: 12 runs, about 2 minutes for the 149-atom CdSe cluster on 4 cores
   python benchmarks/compare_models.py --system tests/CdSe --profile quick

   # standard (about 25 runs) or full (32 runs), with an experimental window
   python benchmarks/compare_models.py --system tests/CdSe --profile full --exp-gap 2.70 2.95

   # your own system, only groups A and B, hexane
   python benchmarks/compare_models.py --system path/to/dir --only A,B --eps-solvent 1.89

   # see which runs a profile contains
   python benchmarks/compare_models.py --system tests/CdSe --profile full --list

Finished cases are reused, so an interrupted run can be restarted; ``--force`` recomputes them. The
default is spin-free; ``--soc`` keeps SOC on.

Checks and how to read them
---------------------------

Every check is reported as PASS, WARN or FAIL.

**Exact relations.** A FAIL here means a bug.

* ``independent_dft`` lowest excitation = DFT gap, and ``independent_qp`` lowest = QP gap, to within
  2 meV.
* Full BSE S\ :sub:`1` ≤ lowest diagonal element (``diagonal_bse``). This follows from the
  variational principle.
* T\ :sub:`1` ≤ S\ :sub:`1`, and a singlet without exchange equals the triplet.
* QP gap − S\ :sub:`1` > 0 for every bound exciton.
* ``env`` in a matched medium (``eps_out`` = :math:`\epsilon_\infty`): QP gap = DFT gap + bulk
  GW–PBE opening, because the reaction field must vanish.

**Physics checks.** A WARN marks a known model limitation or an unconverged setting.

* *Optical gap insensitive to solvent*: the ratio :math:`|\Delta S_1|/|\Delta E_{QP}|` between vacuum
  and solvent. It is < 0.2 for ``env``. It is ≈ 1 for the scissor models, a known limitation rather
  than a bug.
* *Binding robust to the W model*: the spread of QP − S\ :sub:`1` over kernels. For CdSe this spread is
  more than 1 eV. It shows how much the choice of :math:`W` matters, and it is why the ``resta``
  kernel is the documented default.
* *Converged in the active space*: the change of S\ :sub:`1` between the two largest windows is
  < 10 meV.
* *Within experiment*: only when ``--exp-gap`` is given.

Reference results (CdSe 2 nm, full profile)
-------------------------------------------

The full profile on ``tests/CdSe`` produced the tables in
``audit/AUDIT_2026-09-25_QP_EXCITED_STATES.md`` §7. All exact relations pass. The WARNs are the
environment check of the scissor models, the kernel spread, and the slow active-space convergence
(binding still +8 meV per doubling at 100 × 100).
