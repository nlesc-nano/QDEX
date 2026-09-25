Qp screening
============

Part of :doc:`/interactions/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.build_dim_screening_factors``
* CLI: ``--2e-integrals, --kernel, --eps-out``
* YAML: ``physics.2e-integrals, physics.kernel, physics.eps_out``

.. code-block:: python

   build_dim_screening_factors(coords, atom_symbols, material_name=None, eps_out=2.4, alpha=1.0)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:318-355``

Screening Formulations for :math:`W^{\mathrm{QD}}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``QDEX`` implements three distinct physical methods for evaluating :math:`W^{\mathrm{QD}}`:

1. **Discrete Dipole Interaction Model (``sgw-dim``)**:
   Models the semiconductor lattice as an array of atom-centered polarizable point dipoles :math:`\boldsymbol{\mu}_A = \alpha_A \mathbf{E}^{\mathrm{tot}}(\mathbf{R}_A)`.
   The total local field includes external fields and self-consistent dipolar interactions:

   .. math::

      \boldsymbol{\mu}_A = \alpha_A \left[ \mathbf{E}_0(\mathbf{R}_A) - \sum_{B \neq A} \mathbf{T}_{AB} \boldsymbol{\mu}_B \right]

   where :math:`\mathbf{T}_{AB} = \frac{3 \mathbf{R}_{AB} \otimes \mathbf{R}_{AB} - R_{AB}^2 \mathbf{I}}{R_{AB}^5} f_{\mathrm{Thole}}(R_{AB})` is the Thole-damped dipole interaction tensor.
   Solving the linear system :math:`(\mathbf{I} + \boldsymbol{\alpha}\mathbf{T}) \boldsymbol{\mu} = \boldsymbol{\alpha}\mathbf{E}_0` yields the microscopic screened interaction :math:`W_{AB}^{\mathrm{QD}}`.

2. **Resta-Penn Screened Dielectric Model (``sgw-resta``)**:
   Evaluates distance-dependent electronic screening based on the Thomas-Fermi model of the valence electron gas:

   .. math::

      W_{AB}^{\mathrm{QD}} = \frac{1}{\epsilon_{\mathrm{in}}(R) r_{AB}} + \frac{1 - \epsilon_{\mathrm{in}}(R)^{-1}}{r_{AB}} \exp\left( -\frac{r_{AB}}{\lambda_{\mathrm{TF}}} \right)

   where quantum confinement suppresses the core dielectric constant according to the Penn model:

   .. math::

      \epsilon_{\mathrm{in}}(R) = 1 + (\epsilon_\infty - 1) \frac{1}{1 + (R_p / R)^2}, \quad R_p = \frac{\pi}{2} \left( \frac{13.6\text{ eV}}{E_g^{\mathrm{bulk}}} \right) a_0.

3. **Simplified BSE Screened Interaction (``sgw``)**:
   Implements the Cho, Bintrim, and Berkelbach [JCTC 18, 3054 (2022)] polarizability kernel:

   .. math::

      W = (\mathbf{I} + \mathbf{J}_{\mathrm{solv}} \boldsymbol{\Pi}^0)^{-1} \mathbf{J}_{\mathrm{solv}}

   where :math:`\boldsymbol{\Pi}^0` is the non-interacting transition polarizability matrix.
