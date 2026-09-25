Environmental screening and optical cancellation
================================================

The Resta and DIM quasiparticle paths include an exterior dielectric reaction
term through ``eps_out``. Their BSE direct kernels currently use the interior
screening model and do not receive that same exterior reaction operator. This
is an explicit approximation of the implemented model, not a general
consequence of the Delerue cancellation argument.

For a common linear reaction operator :math:`\delta W`, the environmental
contribution to the energy of a neutral electron-hole pair can be written,
within a static charge-density approximation, as

.. math::

   \delta E_{\mathrm{opt}}
   = \tfrac12(\rho_e-\rho_h)^T\delta W(\rho_e-\rho_h).

The quasiparticle gap contains the two separate charging contributions;
the electron-hole interaction supplies the cross term. Complete cancellation
requires the electron and hole densities to be sufficiently similar under the
same reaction operator. Spatially separated carriers, surface-localized
states, and different solvent or cavity models generally leave a residual
optical shift. The Si-nanocrystal result of Delerue, Lannoo and Allan is
evidence for strong cancellation in the systems they studied, not a universal
identity for CdSe or every QDEX kernel.

For an environment-dependent benchmark, report the QP gap, direct binding,
and optical excitation separately at each ``eps_out``. Keep the same
screening/reference convention in the charged and neutral calculations.
The currently unchanged BSE direct kernel does not by itself verify a
compensating change in binding.

Implementation: ``qdex.hardness.estimate_sgw_dim_qp_gap`` and
``estimate_sgw_resta_qp_gap`` contain the exterior reaction term;
``qdex.hardness.build_dim_screening_factors`` and ``build_xs_kernel`` supply
interior BSE screening. The relevant CLI option is ``--eps-out`` and the YAML
key is ``physics.eps_out``.

See also :doc:`/interactions/environment` and :doc:`/interactions/screening`.
