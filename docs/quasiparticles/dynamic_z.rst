Dynamic z
=========

Part of :doc:`/quasiparticles/index`.

.. important::

   ``compute_dynamic_z`` is an empirical scalar damping formula. It does not evaluate the derivative of a computed frequency-dependent self-energy or enforce an f-sum rule. The default pole energy is 15 eV; material-specific plasmon energies are not read from a database.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.compute_dynamic_z``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   compute_dynamic_z(delta_sigma_stat_ev, gap_ev, eps_eff, material_name=None, omega_p_ev=15.0)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:469-471``

8. Avenue 3: Dynamic Renormalization & Self-Consistency
-------------------------------------------------------


.. rubric:: From ``docs/part3_gw_scissor/index.rst:472-488``

Dynamic Renormalization Factor :math:`Z_p` (``--dynamic_z``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than assuming a fixed empirical renormalization factor (such as :math:`Z \approx 0.8`), ``QDEX`` can compute state-dependent dynamic weights :math:`Z_p` from an **empirical plasmon-energy damping ansatz** (no frequency-dependent self-energy or :math:`f`-sum-rule constraint is evaluated):

.. math::

   Z_p^{\mathrm{model}} = \left( 1 + \frac{\max(0,\Delta \Sigma_p^{\mathrm{stat}})}{\tilde{\Omega}_p} \right)^{-1}

where :math:`\Delta \Sigma_p^{\mathrm{stat}}` is the static COHSEX self-energy shift and :math:`\tilde{\Omega}_p` is the screened plasmon frequency:

.. math::

   \tilde{\Omega}_p = \sqrt{ \frac{\Omega_p^2}{\max(1.0, \epsilon_{\mathrm{eff}} - 1.0)} + E_g^2 }.

Here :math:`\Omega_p = \sqrt{4\pi n_v e^2 / m_e}` is the valence electron plasmon energy (:math:`\approx 15 - 20\text{ eV}`). In large quantum dots with strong screening (:math:`\epsilon_{\mathrm{eff}} \gg 1`), :math:`Z_p \to 0.80 - 0.85`; in ultra-small clusters with suppressed screening, :math:`Z_p \to 0.90 - 0.95`.


.. rubric:: From ``docs/part3_gw_scissor/index.rst:489-498``

Implementation in QDEX (``--dynamic_z``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dynamic state-dependent renormalization factors are evaluated by :func:`qdex.hardness.compute_dynamic_z`:

- Uses the optional ``omega_p_ev`` argument (15 eV by default); no material plasmon-energy table is queried.
- Evaluates the screened plasmon frequency :math:`\tilde{\Omega}_p = \sqrt{\frac{\Omega_p^2}{\max(1.0, \epsilon_{\mathrm{eff}} - 1.0)} + E_g^2}` from the current gap and effective screening constant :math:`\epsilon_{\mathrm{eff}}`.
- Computes state renormalization factors :math:`Z_H = (1 + \sigma_H^{\mathrm{stat}} / \tilde{\Omega}_p)^{-1}` and :math:`Z_L = (1 + \sigma_L^{\mathrm{stat}} / \tilde{\Omega}_p)^{-1}`.
- Clamps values to the physically robust interval :math:`Z \in [0.50, 0.99]` to safeguard against numerical unphysicalities in ultra-small clusters.
