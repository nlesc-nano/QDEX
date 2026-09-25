Foundations
===========

Part of :doc:`/quasiparticles/index`.

.. figure:: /_static/figures/qp_hierarchy.svg
   :width: 100%
   :alt: qp hierarchy

   Map of the ``qp_gap`` options and what each one corrects.


.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)


Standard semi-local Kohn-Sham Density Functional Theory (DFT) using functionals like PBE severely underestimates the fundamental band gap of semiconductor nanostructures. For instance, PBE predicts a band gap of :math:`\approx 1.5\text{ eV}` for bulk :math:`\text{CsPbBr}_3`, whereas the experimental quasiparticle gap is :math:`\approx 2.35\text{ eV}`.

In ``QDEX``, single-particle excitation energies are corrected via an analytical and physically grounded hierarchy of **Scaled GW Quasiparticle Models** that incorporate quantum confinement, microscopic electrostatic screening, dielectric solvation, dynamic renormalization, and orbital relaxation.


1. Quasiparticle Theory and Hedin's GW Approximation
-----------------------------------------------------

In many-body perturbation theory, the true single-particle electron addition (electron affinity, EA) and removal (ionization potential, IP) energies are determined by the **quasiparticle (QP) equation**:

.. math::

   \left( \hat{T} + \hat{V}_{\mathrm{ext}}(\mathbf{r}) + \hat{V}_H(\mathbf{r}) \right) \psi_k(\mathbf{r}) + \int \Sigma(\mathbf{r}, \mathbf{r}'; \varepsilon_k^{\mathrm{QP}}) \psi_k(\mathbf{r}') \, d\mathbf{r}' = \varepsilon_k^{\mathrm{QP}} \psi_k(\mathbf{r})

where:

* :math:`\hat{V}_H(\mathbf{r})` is the classical Hartree potential.
* :math:`\Sigma(\mathbf{r}, \mathbf{r}'; \omega)` is the non-local, energy-dependent **electron self-energy** operator.


Hedin's Equations & The :math:`G_0W_0` Approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Hedin's pentagon of many-body equations, the self-energy is expanded to first order in the dynamically screened Coulomb interaction :math:`W`:

.. math::

   \Sigma(\mathbf{r}, \mathbf{r}'; \omega) = \frac{i}{2\pi} \int e^{i \omega' \eta} G_0(\mathbf{r}, \mathbf{r}'; \omega + \omega') W_0(\mathbf{r}, \mathbf{r}'; \omega') \, d\omega'

where:

* :math:`G_0` is the non-interacting single-particle Green's function.
* :math:`W_0 = \epsilon^{-1} v` is the screened Coulomb interaction, mediated by the dynamic microscopic dielectric function :math:`\epsilon(\mathbf{r}, \mathbf{r}'; \omega)`.

The quasiparticle energy shift :math:`\Delta \varepsilon_k^{\mathrm{QP}} = \varepsilon_k^{\mathrm{QP}} - \varepsilon_k^{\mathrm{DFT}}` is given to first order by:

.. math::

   \Delta \varepsilon_k^{\mathrm{QP}} = Z_k \, \operatorname{Re}\langle \psi_k | \Sigma(\varepsilon_k^{\mathrm{DFT}}) - v_{xc}^{\mathrm{DFT}} | \psi_k \rangle

where :math:`Z_k = \left( 1 - \left. \frac{\partial \operatorname{Re}\Sigma}{\partial \omega} \right|_{\varepsilon_k} \right)^{-1}` is the quasiparticle renormalization weight.

