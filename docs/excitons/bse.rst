Coupled BSE (``bse``)
=====================

Part of :doc:`/excitons/index`.

.. code-block:: yaml

   excitations:
     mode: bse
     nhomos: 25
     nlumos: 25
     nroots: 40
     full_diag: false          # Davidson; true = dense diagonalization

**Energies.** The full Tamm–Dancoff Bethe–Salpeter matrix is diagonalized:

.. math::

   A_{ia,jb} = (\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}})\,\delta_{ij}\delta_{ab}
   + 2K^x_{ia,jb} - K^d_{ia,jb},\qquad \mathbf A\,\mathbf X_S = \Omega_S\,\mathbf X_S ,

with :math:`-K^d` only for triplets and :math:`K^x - K^d` in the spinor basis with SOC. The intensities
follow from the exciton transition dipole :math:`\boldsymbol\mu_S = \sum_{ia} X^S_{ia}\boldsymbol\mu_{ia}`.

**What mixing adds over** :doc:`diagonal_bse`.

* **Correlation.** The off-diagonal K\ :sup:`d` lets the electron and hole correlate their positions,
  which lowers S₁ further.
* **Intensity.** The off-diagonal K\ :sup:`x` collects oscillator strength into the bright state.
* **Fine structure.** Near-degenerate transitions split into bright and dark states; with SOC this
  gives the band-edge exciton fine structure.

**Tamm–Dancoff approximation.** The full BSE also couples excitations to de-excitations (the B
block). The TDA neglects it. For band-edge excitons of semiconductor dots the error is below about
0.1 eV; the problem becomes Hermitian and half the size.

**Solver.** For the lowest roots QDEX uses a block-Davidson method: the matrix is never stored, only
products :math:`\mathbf A\mathbf v` are formed; the subspace is expanded with the residuals
preconditioned by :math:`(\mathrm{diag}\,\mathbf A - \omega)^{-1}` until the residual norm is below
``tol`` (10⁻⁵). Memory is :math:`\mathcal O(N_{\mathrm{pairs}}\,n_{\mathrm{roots}})` instead of
:math:`\mathcal O(N_{\mathrm{pairs}}^2)`. ``full_diag: true`` diagonalizes the dense matrix and returns
all states.

**Active space.** S₁ converges slowly with the number of transitions: 2.499 → 2.443 eV at 2 nm from
25 × 25 to 100 × 100. Keep ``nhomos``/``nlumos`` fixed within a comparison.

Perturbative transition selection
---------------------------------

.. code-block:: yaml

   excitations:
     nhomos: 50
     nlumos: 50
     selection: perturbative
     selection_energy: 3.3     # E_thr in eV, on the diagonal A_ia,ia
     # selection_pt: 1.0e-4    # t in hartree

Grimme's selection (J. Chem. Phys. 138, 244104 (2013); ``ptselect`` in std2), for the coupled solvers
``bse``, ``sbse`` and ``stda`` in the spin-free calculation:

1. **Primary transitions**: every ia in the active space with
   :math:`A_{ia,ia} = \Delta E_{ia} + k_x K^x_{ia,ia} - K^d_{ia,ia} \le E_{\mathrm{thr}}`.
2. **Added transitions**: a higher transition u is kept if its second-order coupling to the primaries
   exceeds t,

   .. math::

      E_u^{(2)} = \sum_{v\in P}\frac{|A_{uv}|^2}{A_{uu} - A_{vv}} > t .

3. **Rejected transitions** lower the diagonal of the primaries by their second-order contributions,
   :math:`A_{vv} \to A_{vv} - \sum_{u\,\mathrm{rejected}} |A_{uv}|^2/(A_{uu}-A_{vv})`.

The active space (``nhomos``, ``nlumos``, ``e_thresh``) is the pool from which the selection draws; it
replaces Grimme's orbital window. The selection lets a large pool be used at the cost of a small
matrix. Defaults are those of std2 (E_thr = 7 eV, t = 10⁻⁴ E\ :sub:`h`); for band-edge excitons of dots
E_thr should be set just above the lowest states.

**Cd₆₈Se₅₅Cl₂₆** (sBSE, bulk Resta W, vacuum), pool 50 × 50 = 2,500 transitions:

.. list-table::
   :header-rows: 1

   * - Setting
     - transitions
     - S₁ (eV)
     - f(S₁)
   * - full 25 × 25
     - 625
     - 2.8001
     - 0.303
   * - full 50 × 50
     - 2,500
     - 2.7923
     - 0.241
   * - E_thr = 3.3 eV
     - 37 + 36
     - 2.7855
     - 0.390
   * - E_thr = 3.6 eV
     - 96 + 133
     - 2.7863
     - 0.336
   * - E_thr = 3.3 eV, t = 10⁻⁵
     - 37 + 447
     - 2.7884
     - 0.288

S₁ is within 4–7 meV of the full 50 × 50 result with 3–20 % of the transitions, approaching it from
below (the second-order lowering of the primaries slightly overshoots). Oscillator strengths converge
more slowly: they need the larger selections. The diagonal solvers and the SOC calculation do not use
the selection.
