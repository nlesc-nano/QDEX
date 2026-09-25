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
   :math:`\lambda_{\mathrm{TF}} = d_{\mathrm{NN}}/\sqrt{\epsilon_{\mathrm{in}}-1}`. The size
   dependence of :math:`\epsilon_{\mathrm{in}}` follows the Penn model (Penn, Phys. Rev. 128, 2093
   (1962); Tsu, Babić and Ioriatti, J. Appl. Phys. 82, 1327 (1997)). With one oscillator,
   :math:`\epsilon-1=(\hbar\omega_p/E_P)^2`, and confinement opens the average gap
   :math:`E_P\to E_P+\Delta E_{\mathrm{conf}}`:

   .. math::

      \epsilon_{\mathrm{in}}(R) = 1 + (\epsilon_\infty - 1)\left[\frac{E_P}{E_P+\Delta E_{\mathrm{conf}}}\right]^2,
      \qquad E_P = \frac{\hbar\omega_p}{\sqrt{\epsilon_\infty-1}},
      \qquad \Delta E_{\mathrm{conf}} = E_g^{\mathrm{DFT,QD}} - E_g^{\mathrm{PBE,bulk}}.

   :math:`\hbar\omega_p = \sqrt{4\pi n_v e^2/m}` is the free-electron plasmon of the bulk valence
   (s, p) electron density (``qdex.hardness.valence_plasmon_ev``; 14.1 eV and :math:`E_P` = 6.2 eV for
   CdSe). The iterated ``evgw-resta``/``qsgw-resta`` paths use the current QP gap minus the bulk GW
   gap as :math:`\Delta E_{\mathrm{conf}}`. For CdSe this gives :math:`\epsilon_{\mathrm{in}}` = 3.97
   for Cd₁₆Se₁₃Cl₆ (1.2 nm) and 5.1 for the 2 nm test cluster.

   The same :math:`\hbar\omega_p` and :math:`\epsilon_{\mathrm{in}}` fix the plasmon-pole frequency of
   the derived Z (:doc:`/quasiparticles/dynamic_z`), so one oscillator parameterizes the screening and
   its dynamics.

   .. note::

      Earlier versions used
      :math:`\epsilon_{\mathrm{in}} = 1 + (\epsilon_\infty-1)/[1+(\Delta E_{\mathrm{conf}}/E_g^{\mathrm{PBE,bulk}})^2]`.
      Its reference energy was the 0.64 eV PBE band gap instead of the Penn gap. That gave
      :math:`\epsilon_{\mathrm{in}}` = 1.48 at 1.2 nm and 2.98 at 2 nm, a much stronger reduction than
      atomistic calculations (Wang and Zunger, PRL 73, 1039 (1994); Delerue, Lannoo and Allan, PRB 68,
      115411 (2003)).

3. **Simplified BSE Screened Interaction (``sgw``)**:
   Implements the Cho, Bintrim, and Berkelbach [J. Chem. Theory Comput. 18, 3438 (2022)] polarizability kernel:

   .. math::

      W = (\mathbf{I} + \mathbf{J}_{\mathrm{solv}} \boldsymbol{\Pi}^0)^{-1} \mathbf{J}_{\mathrm{solv}}

   where :math:`\boldsymbol{\Pi}^0` is the non-interacting transition polarizability matrix.
