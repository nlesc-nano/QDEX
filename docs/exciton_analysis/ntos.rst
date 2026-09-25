Ntos
====

Part of :doc:`/exciton_analysis/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.nto``
* Callable: ``qdex.nto.compute_nto_pairs``
* CLI: ``--nto, --nto-states``
* YAML: ``analysis.nto, analysis.nto_states``

.. code-block:: python

   compute_nto_pairs(x_mat)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part5_exciton_analysis/index.rst:161-184``

3. Natural Transition Orbitals (NTOs)
-------------------------------------

In large systems, the BSE expansion vector :math:`X_{ia}^S` often contains hundreds of small non-zero coefficients, obscuring the primary orbital character of the excitation.

**Natural Transition Orbitals (NTOs)** solve this problem by applying a **Singular Value Decomposition (SVD)** to the transition density matrix :math:`\mathbf{X}^S`:

.. math::

   \mathbf{X}^S = \mathbf{U} \, \boldsymbol{\Lambda} \, \mathbf{V}^\dagger = \sum_{k=1}^{K_{\max}} \lambda_k \, \mathbf{u}_k \, \mathbf{v}_k^\dagger

where:
* :math:`\lambda_k \ge 0` are the singular values, normalized such that :math:`\sum_k \lambda_k^2 = 1`.
* :math:`\mathbf{u}_k` defines the :math:`k`-th **NTO hole orbital**: :math:`\tilde{\phi}_k^h(\mathbf{r}) = \sum_i U_{ik} \phi_i(\mathbf{r})`.
* :math:`\mathbf{v}_k` defines the :math:`k`-th **NTO electron orbital**: :math:`\tilde{\phi}_k^e(\mathbf{r}) = \sum_a V_{ak} \phi_a(\mathbf{r})`.

The complex BSE excitation is transformed into a compact sum of uncoupled transitions:

.. math::

   |\Psi_S\rangle = \sum_{k} \lambda_k \, |\tilde{\phi}_k^h \to \tilde{\phi}_k^e\rangle

In most semiconductor quantum dots, a single NTO pair accounts for :math:`> 85\%` of the total transition weight (:math:`\lambda_1^2 > 0.85`), providing an immediate visual interpretation of the excitation.


.. rubric:: From ``docs/part5_exciton_analysis/index.rst:185-195``

NTO Compactness Metrics
~~~~~~~~~~~~~~~~~~~~~~~

``QDEX`` outputs comprehensive NTO statistics:
* **Lead Weight** (:math:`\lambda_1^2`): Fraction of the transition described by the dominant NTO pair.
* **NTO Participation Ratio**: :math:`\mathrm{PR}_{\mathrm{NTO}} = \frac{1}{\sum_k \lambda_k^4}`.
* **NTO Shannon Entropy**: :math:`S_{\mathrm{NTO}} = -\sum_k \lambda_k^2 \ln(\lambda_k^2)`.
* :math:`N_{90}` and :math:`N_{99}`: Number of NTO pairs required to recover 90% and 99% of the transition weight.

---
