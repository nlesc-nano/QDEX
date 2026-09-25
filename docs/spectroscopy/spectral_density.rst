Spectral density
================

Part of :doc:`/spectroscopy/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.analysis``
* Callable: ``qdex.namd.analysis.compute_spectral_density``
* CLI: ``--sigma``
* YAML: ``namd.transient_absorption``

.. code-block:: python

   compute_spectral_density(signal, dt_fs)


7. Phonon Spectral Density J(ω): Mapping Electron-Phonon Coupling
-----------------------------------------------------------------


Mathematical Definition
~~~~~~~~~~~~~~~~~~~~~~~

The **Phonon Spectral Density** :math:`J(\omega)` is the Fourier transform of the energy gap autocorrelation function:

.. math::

   J(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} C(t) \, e^{i \omega t} \, dt

``QDEX`` stores the non-negative real part of a Hann-windowed Fourier transform of :math:`C(t)`, in :math:`\mathrm{cm}^{-1}`. That is the cosine transform of a real, even autocorrelation. It is not :math:`|\mathrm{FFT}|^2`. On a record of length :math:`T` the Rayleigh spacing is :math:`1/T` (about :math:`33\,\mathrm{cm}^{-1}` for a 1 ps trajectory at 2 fs), and that spacing is printed next to the peak.


Physical Meaning
~~~~~~~~~~~~~~~~

:math:`J(\omega)` provides the **frequency-resolved spectrum of nuclear vibrations that couple to the electronic transitions**. The area under a peak at frequency :math:`\omega` is proportional to the electron-phonon coupling strength (Huang-Rhys parameter :math:`S_\alpha`) for that vibrational mode:

.. math::

   J(\omega) = \pi \sum_\alpha \omega_\alpha^2 \, S_\alpha \, \delta(\omega - \omega_\alpha)


Identifying Active Phonon Modes During Cooling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By inspecting the peaks in :math:`J(\omega)`, researchers can directly identify which phonon modes facilitate carrier cooling and dephasing:

* **Low-Frequency Acoustic Modes** (:math:`< 50\text{ cm}^{-1}`): Acoustic phonons provide a continuous low-energy bath that mediates intra-band thermalization within dense manifolds.
* Optical Modes (e.g. Pb–Br / Pb–I Stretching, :math:`60 - 150\text{ cm}^{-1}`): Polar optical phonons create strong macroscopic electric fields (Fröhlich interaction), driving fast non-adiabatic transitions across intermediate energy gaps.
* **Organic Cation / Ligand Modes** (:math:`200 - 300\text{ cm}^{-1}`): In hybrid perovskites (:math:`\text{MAPbI}_3`), rotational and librational motions of methylammonium cations produce high-frequency peaks in :math:`J(\omega)` that help bridge larger energy spacings.

