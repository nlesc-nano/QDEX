Nacs tracking
=============

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.precompute``
* Callable: ``qdex.namd.precompute.align_phases_and_crossings``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   align_phases_and_crossings(S_mat, C_next, track_crossings=True, lock_above=0.5)


5. Trajectory Precomputation & Wavefunction Tracking
-----------------------------------------------------


Numerical Non-Adiabatic Couplings (NAC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Along the classical nuclear trajectory :math:`\mathbf{R}(t)`, the non-adiabatic coupling is evaluated numerically via finite differences:

.. math::

   d_{IJ}(t + \frac{\Delta t}{2}) = \langle \psi_I(t) | \frac{\partial}{\partial t} | \psi_J(t) \rangle \approx \frac{S_{IJ}(t, t+\Delta t) - S_{JI}^*(t, t+\Delta t)}{2 \Delta t}

where :math:`S_{IJ}(t, t+\Delta t) = \langle \psi_I(t) | \psi_J(t+\Delta t) \rangle` is the cross-frame state overlap. In ``QDEX``, the underlying atomic orbital cross-overlaps :math:`S_{\mu \nu}(t, t+\Delta t) = \int \chi_\mu(\mathbf{r}; \mathbf{R}(t)) \chi_\nu(\mathbf{r}; \mathbf{R}(t+\Delta t)) d\mathbf{r}` are evaluated analytically via Libint2.


Eliminating Gauge Phase Discontinuities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard electronic eigensolvers determine eigenvectors up to an arbitrary global phase factor :math:`e^{i \theta_I(t)}`. If left uncorrected, random sign flips between successive MD frames cause :math:`S_{II}(t, t+\Delta t) \approx -1`, producing spurious non-adiabatic couplings that are orders of magnitude too large.

``QDEX`` eliminates gauge discontinuities by applying a phase rotation:

.. math::

   |\psi_I(t+\Delta t)\rangle \leftarrow |\psi_I(t+\Delta t)\rangle \, e^{-i \theta_I}

where :math:`\theta_I = \operatorname{arg}(S_{II}(t, t+\Delta t))`. This guarantees that the diagonal overlap is strictly real and positive:

.. math::

   \operatorname{Re}(S_{II}(t, t+\Delta t)) \ge 0, \quad \operatorname{Im}(S_{II}(t, t+\Delta t)) = 0


Hungarian Matching for Trivial Crossings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A trivial crossing is a label swap between two orbitals that do not interact. Inside one nuclear step their diagonal overlap collapses and the character sits on an off-diagonal element. Following the old energy label then puts a numerical spike into :math:`d_{IJ}`.

``QDEX`` repairs only that case. The Hungarian assignment uses the cost matrix

.. math::

   C_{IJ} = 1 - |S_{IJ}(t, t+\Delta t)|^2

and a state is allowed to leave its own column only when :math:`|S_{II}| < 0.5`. An avoided crossing that still overlaps its own adiabatic label stays in the adiabatic basis that surface hopping propagates. The assignment is not a global diabatization.

