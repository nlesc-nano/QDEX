Fcwd
====

Part of :doc:`/recombination/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.compute_fcwd_rate``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   compute_fcwd_rate(E_gap_ev, V_el_ev, lambda_ev, sigma_ev)


6. Intermediate & Narrow Gap Decay: Franck-Condon Weighted Density of States (FCWD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In narrow-gap quantum dots (:math:`E_g < 1.0\text{ eV}`, e.g. infrared PbS, PbSe, HgTe, or InAs) or near crossing points with shallow defect states, the number of required accepting phonons is small (:math:`p < 15`). In this regime, the discrete multi-phonon expansion transitions into a continuous **Franck-Condon Weighted Density of States (FCWD)** evaluated from Fermi's Golden Rule:

.. math::

   k_{\mathrm{nr}}^{\mathrm{FCWD}} = \frac{2\pi}{\hbar} |V_{\mathrm{el}}|^2 \, \mathrm{FCWD}(E_g)

In the Marcus-Levich framework, the nuclear Franck-Condon factor takes a Gaussian line-shape:

.. math::

   \mathrm{FCWD}(E_g) = \frac{1}{\sqrt{2\pi \sigma^2}} \exp\left( -\frac{(E_g - \lambda)^2}{2\sigma^2} \right)

Data-Driven Extraction of FCWD Parameters
""""""""""""""""""""""""""""""""""""""""""

Every variable in this rate expression is evaluated directly from the NAMD trajectory:

1. **Gaussian Broadening** (:math:`\sigma`):
   :math:`\sigma` is **not an arbitrary broadening parameter**! It is the exact standard deviation of the energy gap fluctuations from the MD trajectory:

   .. math::

      \sigma = \sigma_E = \sqrt{\langle (E_g(t) - \langle E_g \rangle)^2 \rangle} = \sqrt{2 \lambda k_B T}

2. **Nuclear Reorganization Energy** (:math:`\lambda`):
   Computed from the trajectory variance: :math:`\lambda = \frac{\sigma_E^2}{2 k_B T}`.
3. **Electronic prefactor**:
   The wide-gap loss uses the Englman–Jortner law with the trajectory Huang–Rhys factor :math:`S`, the phonon energy taken from the peak of :math:`J(\omega)`, and a configured prefactor :math:`A_{\mathrm{nr}}` (default :math:`10^{13}\,\mathrm{s}^{-1}`). The root-mean-square intraband coupling :math:`\hbar\langle|d|\rangle` connects excited states to each other. It is not the exciton-to-ground matrix element, and it is not inserted as :math:`V_{\mathrm{el}}`. The Marcus–Levich Gaussian is appropriate only when the gap is a few :math:`\sigma`. It is not the loss used for a :math:`\mathrm{CsPbBr}_3` gap of about 3 eV.
