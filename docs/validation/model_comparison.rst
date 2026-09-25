QP and BSE compatibility
========================

Part of :doc:`/validation/index`.

.. figure:: /_static/figures/qp_hierarchy.svg
   :width: 100%
   :alt: qp hierarchy

   The ``qp_gap`` options and what each one corrects.

The rule: the same W in GW and BSE
----------------------------------

In GW–BSE, the screened interaction :math:`W` that corrects the quasiparticle energies is the same
:math:`W` that binds the electron and hole in the BSE direct term. A Delta-W QP model computes a
self-energy correction from a model :math:`W_{QD}` (including its solvent term). Its excitations are
consistent only if the BSE kernel is that same :math:`W_{QD}`. Otherwise:

* the charged excitation (QP gap) sees the solvent and the reduced interior screening, while
* the neutral excitation (BSE) does not.

The classical polarization energies then fail to cancel in the optical gap, and S\ :sub:`1` follows
the QP gap one to one when the solvent changes.

QDEX therefore enforces the rule. Each Delta-W QP model exposes the :math:`W` it built, and the BSE
uses it through ``kernel: qp``.

Compatibility table
-------------------

.. list-table::
   :header-rows: 1
   :widths: 22 38 20 20

   * - ``qp_gap``
     - :math:`W` of the QP model
     - allowed ``kernel``
     - default ``kernel``
   * - ``sgw-resta``, ``evgw-resta``, ``qsgw-resta``
     - Resta profile with the size-scaled :math:`\epsilon_{\mathrm{eff}}` (final iteration for ``evgw``/``qsgw``) plus solvent term
     - ``qp``
     - ``qp``
   * - ``sgw-resta-pure``
     - Resta profile with bulk :math:`\epsilon_\infty` plus solvent term
     - ``qp``
     - ``qp``
   * - ``sgw-dim``, ``evgw-dim``, ``qsgw-dim``
     - DIM/Thole-scaled profile (final iteration for ``evgw``/``qsgw``) plus solvent term
     - ``qp``
     - ``qp``
   * - ``sgw`` (atom mode)
     - sBSE RPA :math:`W`, solvent inside :math:`J`
     - ``qp`` (``sbse`` is accepted as the same matrix)
     - ``qp``
   * - ``gw`` / ``sgw-anchor``, ``brus``, ``pbe``, numeric gap
     - none (gap-only models)
     - ``resta``, ``xs-resta``, ``dim``, ``sbse``, ``bse``
     - ``bse`` (legacy)

Rules enforced by the CLI:

* A Delta-W model combined with any kernel other than ``qp`` is an error. The escape hatch
  ``--allow-inconsistent-kernel`` exists only to reproduce old results.
* ``kernel: qp`` with a gap-only model is an error, because there is no :math:`W` to share.
* ``kernel: qp`` is atom-resolved and requires ``two_electron_integrals: mnok``.
* ``qp_z`` / ``dynamic_z`` with a gap-only model is an error; these models have no Z.

For the gap-only models the kernel is an independent modelling choice. The gap model and the kernel do
not share a :math:`W`, and the solvent term of ``gw`` does not enter the BSE.

The quasiparticle weight Z
--------------------------

The Delta-W correction is :math:`\Delta\varepsilon_p = Z\cdot\tfrac12 q_p^T\Delta W q_p`. The BSE uses
the full :math:`W`. ``qp_z`` selects Z:

* ``qp_z: 1.0``: the static limit. The classical polarization terms then cancel between QP gap and
  BSE, and S\ :sub:`1` is independent of the solvent. **Recommended with a static BSE.**
* ``qp_z: 0.8`` (default of the one-shot models): scales the QP shift but not the BSE attraction. A
  fraction :math:`1-Z` of the polarization energy is left over, so S\ :sub:`1` shifts with the solvent.
  In the opposite direction to the old mismatched combination, and 5× smaller.
* ``qp_z: derived``: empirical state-dependent :math:`Z_p` (0.94–0.97 for CdSe). Not a computed
  self-energy derivative.

A reduced Z paired with a static BSE is not a consistent pair. In full GW–BSE, the renormalization of
the QP energies is largely compensated by the dynamical screening of the BSE kernel, which a static
kernel omits.

What consistency does to the results (CdSe 2 nm)
------------------------------------------------

Cd\ :sub:`68`\ Se\ :sub:`55`\ Cl\ :sub:`26`, spin-free, 25 × 25 active space. S\ :sub:`1` in eV.

.. list-table::
   :header-rows: 1

   * - model
     - vacuum QP / S\ :sub:`1`
     - toluene (:math:`\epsilon_{\mathrm{out}}=2.24`) QP / S\ :sub:`1`
   * - ``sgw-resta``, old mismatched kernel
     - 3.762 / 3.535
     - 3.253 / 3.026
   * - ``sgw-resta``, shared W, Z = 0.8
     - 3.762 / 2.240
     - 3.253 / 2.383
   * - ``sgw-resta``, shared W, Z = 1
     - 4.021 / 2.499
     - 3.385 / 2.515
   * - ``sgw-dim``, shared W, Z = 1
     - 3.809 / 2.529
     - 3.173 / 2.531
   * - ``qsgw-dim``, shared W
     - 3.955 / 2.689
     -
   * - ``qsgw-resta``, shared W
     - 4.234 / 2.738
     -

The results show four things:

* **Solvent.** With the shared :math:`W` and Z = 1, S\ :sub:`1` changes by 2–16 meV between vacuum and
  toluene, compared with 0.5 eV for the mismatched combination.
* **W model.** Resta and DIM give the same S\ :sub:`1` to within 30 meV. The interior screening
  contrast raises the QP gap and the binding by almost the same amount, just as the solvent term does.
* **What sets S₁.** With a consistent static :math:`W`, S\ :sub:`1` ≈ DFT gap + bulk GW–PBE opening − binding
  with bulk screening. That is about 2.5 eV spin-free, about 2.4 eV for the first bright state with SOC.
  This is 0.3–0.5 eV below the experimental window of 2.70–2.95 eV (Yu et al. 2003 sizing curve). The
  earlier agreement of the mismatched ``sgw-*`` runs in toluene came from the missing electron–hole
  polarization.
* **Orbital relaxation.** ``qsgw-*`` is higher by about 0.2 eV. Its screened-exchange operator
  :math:`-\tfrac12 P\circ\Delta W` relaxes the orbitals and does not reduce to a classical charging energy,
  so it is not cancelled.

The comparison script
---------------------

``benchmarks/compare_models.py`` runs the ``qdex`` CLI once per case, starting from your
``config.yaml``, and writes ``summary.md``, ``summary.csv`` and ``summary.png``:

.. list-table::
   :header-rows: 1
   :widths: 8 92

   * - group
     - question
   * - A
     - QP gap and S\ :sub:`1` of every ``qp_gap`` option, each with its compatible kernel (vacuum)
   * - B
     - Is S\ :sub:`1` insensitive to the solvent, while the QP gap moves?
   * - C
     - Effect of Z (0.8, 1.0, derived) with the shared :math:`W`
   * - D
     - Independent kernels for the gap-only ``gw`` model and the binding each gives
   * - E
     - Exact relations of the excitation modes, spin and charges
   * - F
     - Active-space convergence
   * - L
     - The legacy mismatched combination, for comparison

.. code-block:: bash

   python benchmarks/compare_models.py --system tests/CdSe --profile quick           # ~12 runs
   python benchmarks/compare_models.py --system tests/CdSe --profile full --soc --exp-gap 2.70 2.95
   python benchmarks/compare_models.py --system path/to/dir --only A,B --eps-solvent 1.89
   python benchmarks/compare_models.py --system tests/CdSe --profile full --list    # list the runs

Finished cases are reused (``--force`` recomputes them). Use ``--soc`` to compare the first bright
state with experiment. The checks are reported as follows:

* **FAIL**: a violated exact relation, which means a bug. Examples: ``independent_qp`` lowest ≠ QP gap;
  full BSE S\ :sub:`1` above the lowest diagonal element; T\ :sub:`1` above S\ :sub:`1`; singlet
  without exchange ≠ triplet; negative binding.
* **WARN**: a model limitation or an unconverged setting. Examples: S\ :sub:`1` moving with the solvent
  (gap-only models, or Z < 1); a large spread of the independent kernels for ``gw``; active-space
  convergence; outside the experimental window.
