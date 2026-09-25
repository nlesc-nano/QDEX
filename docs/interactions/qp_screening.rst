Qp screening
============

Part of :doc:`/interactions/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.build_dim_screening_factors``
* CLI: ``--2e-integrals, --kernel, --eps-out``
* YAML: ``physics.2e-integrals, physics.kernel, physics.eps_out``

.. code-block:: python

   build_dim_screening_factors(coords, atom_symbols, material_name=None, eps_out=2.4, alpha=1.0)


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

   where :math:`1/r_{AB}` is regularized by the MNOK damping and
   :math:`\lambda_{\mathrm{TF}} = d_{\mathrm{NN}}/\sqrt{\epsilon_{\mathrm{in}}-1}`. The implemented
   size dependence of :math:`\epsilon_{\mathrm{in}}` is an interpolation in the confinement energy
   (``estimate_sgw_resta_qp_gap``, one-shot path):

   .. math::

      \epsilon_{\mathrm{in}} = 1 + \frac{\epsilon_\infty - 1}{1 + \left(\Delta E_{\mathrm{conf}} / E_g^{\mathrm{PBE,bulk}}\right)^2},
      \qquad \Delta E_{\mathrm{conf}} = E_g^{\mathrm{DFT,QD}} - E_g^{\mathrm{PBE,bulk}}

   (the iterated ``evgw-resta``/``qsgw-resta`` paths use the current QP gap and the bulk GW gap
   instead).

   .. note::

      Despite its name, this is not Penn's scaling. In the Penn model
      :math:`\epsilon-1\propto(\hbar\omega_p/E_P)^2`, where :math:`E_P` is the average (Penn) gap,
      about 4–5 eV for CdSe, not the fundamental band gap. A Penn-consistent size correction would be
      :math:`\epsilon_{\mathrm{in}}-1=(\epsilon_\infty-1)\,[E_P/(E_P+\Delta E_{\mathrm{conf}})]^2`.
      Using the small PBE bulk gap as the reference energy makes the reduction much stronger. For the
      CdSe test cluster (:math:`\Delta E_{\mathrm{conf}}=0.82` eV) the implemented expression gives
      :math:`\epsilon_{\mathrm{in}}\approx3.0`, compared with :math:`\approx4.7` from the Penn form with
      :math:`E_P=4.5` eV. Atomistic calculations of ~2 nm dots (Wang and Zunger, PRL 73, 1039 (1994);
      Delerue, Lannoo and Allan, PRB 68, 115411 (2003)) find a reduction of tens of percent, not a halving.

3. **Simplified BSE Screened Interaction (``sgw``)**:
   Implements the Cho, Bintrim, and Berkelbach [J. Chem. Theory Comput. 18, 3438 (2022)] polarizability kernel:

   .. math::

      W = (\mathbf{I} + \mathbf{J}_{\mathrm{solv}} \boldsymbol{\Pi}^0)^{-1} \mathbf{J}_{\mathrm{solv}}

   where :math:`\boldsymbol{\Pi}^0` is the non-interacting transition polarizability matrix.
