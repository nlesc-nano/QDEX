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

**Active space.** S₁ converges slowly with the number of transitions: 2.779 → 2.760 → 2.720 eV at 2 nm
(sBSE, vacuum) from 25 × 25 to 50 × 50 to 100 × 100. Keep ``nhomos``/``nlumos`` fixed within a comparison.

Perturbative transition selection
---------------------------------

.. code-block:: yaml

   excitations:
     nhomos: 50
     nlumos: 50
     selection: perturbative
     selection_energy: 3.3     # E_thr in eV, on the diagonal A_ia,ia
     # selection_pt: 1.0e-4    # t in hartree
     # selection_shift: on     # second-order lowering of the primaries (std2)

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

**Cd₆₈Se₅₅Cl₂₆** (sBSE, bulk Resta W, vacuum, default MNOK integrals), pool 50 × 50 = 2,500
transitions:

.. list-table::
   :header-rows: 1

   * - Setting
     - transitions
     - S₁ (eV)
     - f(S₁)
   * - full 25 × 25
     - 625
     - 2.7786
     - 0.254
   * - full 50 × 50
     - 2,500
     - 2.7599
     - 0.195
   * - E_thr = 3.3 eV
     - 43 + 116
     - 2.7461
     - 0.296
   * - E_thr = 3.6 eV
     - 103 + 314
     - 2.7484
     - 0.249
   * - E_thr = 3.3 eV, t = 10⁻⁵
     - 43 + 1,293
     - 2.7576
     - 0.200

With a 100 × 100 pool (10,000 transitions; the full dense calculation takes 245 s, the selected ones
6–18 s):

.. list-table::
   :header-rows: 1

   * - Setting
     - transitions
     - S₁ (eV), with lowering
     - S₁ (eV), ``selection_shift: off``
   * - full 100 × 100
     - 10,000
     - 2.7195
     - 2.7195
   * - E_thr = 3.3 eV
     - 43 + 122
     - 2.6928
     - 2.7889
   * - E_thr = 3.6 eV
     - 103 + 393
     - 2.6959
     - 2.7769
   * - E_thr = 3.3 eV, t = 10⁻⁵
     - 43 + 2,159
     - 2.7110
     - 2.7525

* **50 × 50 pool.** S₁ is within 11–14 meV of the full result with 6–17 % of the transitions, and
  within 2 meV with t = 10⁻⁵ (53 %). It is approached from below: the second-order lowering of the
  primaries slightly overshoots.
* **100 × 100 pool.** The added transitions are nearly the same as with the 50 × 50 pool; the extra
  ones enter only through the second-order lowering (mean 40–50 meV, up to 110 meV per primary). With
  the lowering S₁ is 24–27 meV too low, 9 meV with t = 10⁻⁵ (22 % of the transitions, 18 s instead
  of 245 s). Without it S₁ is 33–69 meV too high. The two bracket the full result; the std2 default
  (lowering on) is clearly the closer one.
* **Oscillator strengths** converge slowly, both with the pool (0.25 → 0.20 → 0.10 from 25 × 25 to
  100 × 100) and with the selection; they need the larger selections.
* The diagonal solvers and the SOC calculation do not use the selection.

**Dynamics.** The selection is meant for optical spectra near the band edge. Carrier cooling needs
every state between the pump energy and the band edge, and the same transitions in every frame, so the
NAMD precompute does not use it: it fixes the transition set on the first frame with an energy window
(``namd.storage.active_energy_window_ev``) and propagates within ``[gap − 0.2, pump + 0.3]`` eV by default.
A perturbative selection would change with geometry from frame to frame, breaking the state tracking
and the couplings, and its second-order shifts would add noise to the energies.
