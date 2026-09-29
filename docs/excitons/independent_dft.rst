Independent KS transitions (``independent_dft``)
================================================

Part of :doc:`/excitons/index`.

.. code-block:: yaml

   excitations:
     mode: independent_dft

**Energies.** Differences of Kohn–Sham eigenvalues, with no QP correction and no electron–hole
interaction:

.. math::

   \Omega_{ia} = \varepsilon_a - \varepsilon_i .

**Intensities.** From the transition dipoles of the orbital pair,

.. math::

   \boldsymbol\mu_{ia} = \langle\psi_i|\mathbf r|\psi_a\rangle,\qquad
   f_{ia} = \tfrac23\,\Omega_{ia}\,|\boldsymbol\mu_{ia}|^2 .

**What it is for.**

* The reference against which every correction is measured.
* The orbital character of the spectrum: which transitions are bright, how dense the manifold is.

**Why it is not the optical gap.** The KS gap is neither the QP gap nor the optical gap. It
underestimates the QP gap (missing self-energy) and has no binding. In small dots the two errors
partly cancel, so the KS gap can land near S₁, but this is fortuitous and does not carry over across
sizes or solvents.

**Cost.** Negligible: no two-electron integrals are built.
