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
