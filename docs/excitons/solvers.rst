Solvers
=======

Part of :doc:`/excitons/index`.

.. important::

   ``bse`` solves the resonant/Tamm-Dancoff matrix. Davidson returns selected roots; ``diagonal_bse`` uses the Hamiltonian diagonal without an iterative eigensolve.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.davidson``
* Callable: ``qdex.davidson.davidson``
* CLI: ``--excitation-mode, --include-direct-eh, --include-exchange``
* YAML: ``physics.excitation_mode, physics.include_direct_eh, physics.include_exchange``

.. code-block:: python

   davidson(matvec, diag, nroots, max_iter=500, tol=1e-06, max_subspace=None, device='numpy')


7. Full BSE & The Davidson Iterative Solver
-------------------------------------------

To solve for the lowest :math:`k` roots without dense :math:`O(N_{\mathrm{pairs}}^3)` matrix diagonalization, ``QDEX`` implements an optimized **Davidson iterative subspace solver**:

1. Projects :math:`\mathbf{A}` into a small trial subspace :math:`\mathbf{V} = [\mathbf{v}_1, \dots, \mathbf{v}_m]`.
2. Computes the matrix-vector product :math:`\mathbf{w}_j = \mathbf{A} \mathbf{v}_j` on-the-fly.
3. Solves the projected eigenvalue problem :math:`\mathbf{V}^\dagger \mathbf{A} \mathbf{V} \mathbf{y} = \omega \mathbf{y}`.
4. Computes the residual :math:`\mathbf{r} = \mathbf{A} \mathbf{x} - \omega \mathbf{x}` and preconditioner :math:`\boldsymbol{\delta} = (\operatorname{diag}(\mathbf{A}) - \omega)^{-1} \mathbf{r}`.
5. Expands the subspace until the norm :math:`||\mathbf{r}|| < \text{tol}` (default :math:`10^{-5}`).

This reduces memory requirements from :math:`O(N_{\mathrm{pairs}}^2)` to :math:`O(N_{\mathrm{pairs}} \times k)`.

