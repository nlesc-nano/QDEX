Screening
=========

Part of :doc:`/interactions/index`.

.. figure:: /_static/figures/screening_profile.svg
   :width: 100%
   :alt: screening profile

   Resta-MNOK screening used by the direct kernel for CdSe (ε∞ = 6.2, d_NN = 2.60 Å). Short-range pairs are unscreened; long-range pairs are screened by 1/ε∞. The MNOK damping keeps the on-site value finite.


.. important::

   Resta and DIM screening builders omit external-medium screening from their BSE direct kernel. This is a model choice; their QP paths include an external reaction term. The two pieces have not been shown to cancel for arbitrary electron and hole densities.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.build_xs_kernel``
* CLI: ``--2e-integrals, --kernel, --eps-out``
* YAML: ``physics.2e-integrals, physics.kernel, physics.eps_out``

.. code-block:: python

   build_xs_kernel(shells, atom_symbols, coords, atom_ao_ranges, material_name=None, kernel_mode='bse', alpha=1.0, nthreads=1, C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None, C_occ_b_low=None, C_virt_b_low=None, eps_occ_b=None, eps_virt_b=None, eps_out=2.4, return_eps_info=False)


4. Dielectric Screening Kernels (``kernel``)
--------------------------------------------

The direct electron-hole attraction :math:`K^d` is mediated by the screened interaction :math:`W`. ``QDEX`` provides five distinct dielectric screening kernels:


1. Resta Screened Dielectric Kernel (``kernel: resta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In semiconductor nanoclusters, dielectric screening varies continuously from the bulk optical dielectric constant :math:`\epsilon_\infty` at large distances down to unscreened vacuum interaction (:math:`\epsilon = 1`) at short, intra-atomic distances.

``QDEX`` implements the **Resta model of electronic screening**:

.. math::

   W(r) = \frac{1}{\epsilon_\infty r} + \frac{1 - \epsilon_\infty^{-1}}{r} \exp\left( -\frac{r}{\lambda_s} \right)

where :math:`\lambda_s` is the **Thomas-Fermi screening length** of the valence electron gas:

.. math::

   \lambda_s = \sqrt{ \frac{\pi}{4 k_F} } = \left( \frac{\pi}{4 (3\pi^2 n_v)^{1/3}} \right)^{1/2}.

Damped over atomic centers with Ohno-Klopman hardness parameters, the discrete kernel is:

.. math::

   W_{AB}^{\mathrm{Resta}} = \frac{1}{\epsilon_\infty \sqrt{R_{AB}^2 + a_{AB}^2}} + \frac{1 - \epsilon_\infty^{-1}}{\sqrt{R_{AB}^2 + a_{AB}^2}} \exp\left( -\frac{\sqrt{R_{AB}^2 + a_{AB}^2}}{\lambda_s} \right).


2. Atomistic Discrete Dipole Interaction Kernel (``kernel: dim``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluates screening via the Discrete Dipole Interaction Model (DIM / Thole model). Each atom :math:`A` responds with an induced dipole :math:`\boldsymbol{\mu}_A = \alpha_A \mathbf{E}^{\mathrm{tot}}(\mathbf{R}_A)` to the electron-hole charge distribution:

.. math::

   (\mathbf{I} + \boldsymbol{\alpha}\mathbf{T}) \boldsymbol{\mu} = \boldsymbol{\alpha}\mathbf{E}_0

yielding an atom-specific screened potential :math:`W_{AB}^{\mathrm{DIM}} = S_{AB} \, \gamma_{AB}`.


3. Parameter-Free Microscopic ZDO-RPA Kernel (``kernel: rpa`` / ``xs-rpa``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Computes microscopic dielectric screening directly from the non-interacting transition polarizability of the quantum dot without any empirical parameters:

.. math::

   \Pi_{\mu \nu}^0 = 4 \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} \frac{C_{\mu i} C_{\mu a} C_{\nu i} C_{\nu a}}{\varepsilon_a - \varepsilon_i}

.. math::

   \boldsymbol{\epsilon} = \mathbf{I} + \boldsymbol{\Gamma} \boldsymbol{\Pi}^0, \quad \mathbf{W} = \boldsymbol{\epsilon}^{-1} \boldsymbol{\Gamma}.


4. Simplified BSE Kernel (``kernel: sbse``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Constructs the screened interaction according to Cho, Bintrim, and Berkelbach:

.. math::

   \mathbf{W} = (\mathbf{I} + \mathbf{J}_{\mathrm{solv}} \boldsymbol{\Pi}^0)^{-1} \mathbf{J}_{\mathrm{solv}}.


5. Uniform Dielectric Kernel (``kernel: bse``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scales the interaction uniformly by :math:`1/\epsilon_\infty` with an empirical scaling factor :math:`\alpha` (default :math:`\alpha = 1.0`).
