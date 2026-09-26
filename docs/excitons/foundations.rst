Foundations
===========

Part of :doc:`/excitons/index`.

.. figure:: /_static/figures/bse_kernels.svg
   :width: 100%
   :alt: bse kernels

   The two kernels of the TDA Bethe–Salpeter matrix. Exchange uses the bare interaction and transition charges :math:`q^{ia}_A`; the direct term uses the screened :math:`W` and the hole (:math:`q^{ij}_A`) and electron (:math:`q^{ab}_B`) densities.


.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.solver``
* Callable: ``qdex.solver.ExcitonSolver.solve``
* CLI: ``--excitation-mode, --include-direct-eh, --include-exchange``
* YAML: ``physics.excitation_mode, physics.include_direct_eh, physics.include_exchange``

.. code-block:: python

   solve(self, nroots=10, full_diag=False, tol=1e-05, excitation_mode='bse')


The description of neutral optical excitations in semiconductor nanostructures requires treating the two-particle correlated motion of an electron promoted to the conduction band and the hole left behind in the valence band.

QDEX provides four frameworks, from non-interacting transitions to the coupled Bethe–Salpeter
equation (BSE) in the Tamm–Dancoff approximation (TDA) (:doc:`frameworks`). The QP energies come from
:doc:`/quasiparticles/index`; the screened interaction of the direct term is the W of the same QP
model (:doc:`screened_kernel`).


1. The Two-Particle Excitation Problem
--------------------------------------

In a neutral optical excitation, an electron is removed from an occupied valence orbital :math:`i` and placed into an unoccupied conduction orbital :math:`a`, creating an electron-hole pair configuration :math:`|ia\rangle = a_a^\dagger a_i |0\rangle`.

The exact correlated excited-state wavefunction :math:`|\Psi_S\rangle` is expressed as a linear superposition of electron-hole configurations:

.. math::

   |\Psi_S\rangle = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, |ia\rangle

where :math:`X_{ia}^S` are the configuration interaction amplitudes and :math:`\Omega_S` is the corresponding optical excitation energy.


Tamm-Dancoff Approximation (TDA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In linear-response many-body perturbation theory, the full Bethe-Salpeter Equation contains both resonant excitations (:math:`\mathbf{A}`) and anti-resonant de-excitations (:math:`\mathbf{B}`):

.. math::

   \begin{pmatrix}
     \mathbf{A} & \mathbf{B} \\
     -\mathbf{B}^* & -\mathbf{A}^*
   \end{pmatrix}
   \begin{pmatrix}
     \mathbf{X} \\
     \mathbf{Y}
   \end{pmatrix}
   = \Omega
   \begin{pmatrix}
     \mathbf{X} \\
     \mathbf{Y}
   \end{pmatrix}.

Under the **Tamm-Dancoff Approximation (TDA)**, coupling to ground-state de-excitations is neglected (:math:`\mathbf{B} \approx \mathbf{0}`). This reduces the problem to a standard Hermitian eigenvalue problem:

.. math::

   \mathbf{A} \mathbf{X}_S = \Omega_S \mathbf{X}_S.

For the low-lying excitons of semiconductor nanocrystals the TDA is usually a small approximation (typically below 0.1 eV for band-edge states). It removes triplet instabilities, gives real excitation energies and a Hermitian problem, and halves the dimension of the eigenvalue problem. It is less reliable for high-energy states and for systems with small gaps relative to the kernel.


Singlet vs. Triplet Matrix Elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The elements of the resonant matrix :math:`A_{ia, jb}` depend on the spin multiplicity:

.. math::

   A_{ia, jb}^{\mathrm{singlet}} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} + 2 K_{ia, jb}^x - K_{ia, jb}^d

.. math::

   A_{ia, jb}^{\mathrm{triplet}} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} - K_{ia, jb}^d.

with the exchange and direct terms

.. math::

   K^x_{ia,jb} = (ia|v|jb) = \iint \psi_i^*(\mathbf r)\psi_a(\mathbf r)\,v(\mathbf r,\mathbf r')\,\psi_j(\mathbf r')\psi_b^*(\mathbf r'),

.. math::

   K^d_{ia,jb} = (ij|W|ab) = \iint \psi_i^*(\mathbf r)\psi_j(\mathbf r)\,W(\mathbf r,\mathbf r')\,\psi_a(\mathbf r')\psi_b^*(\mathbf r').

K\ :sup:`x` uses the bare Coulomb interaction and K\ :sup:`d` the static screened interaction of the
QP model. Their representation (MNOK or ZDO xs) is the one used for the QP correction
(:doc:`/quasiparticles/representation`).

The exchange term :math:`K^x` describes virtual annihilation of the electron–hole pair and its recreation elsewhere (the local-field term). In a closed-shell reference the spin-adapted singlet combination :math:`(|i\alpha\to a\alpha\rangle+|i\beta\to a\beta\rangle)/\sqrt2` couples to this process with weight 2, whereas the triplet combinations have zero net transition density and do not couple at all. This gives the factor :math:`2K^x` for singlets, :math:`0` for triplets, and the singlet–triplet (dark–bright, before SOC) splitting.


Relativistic 2-Component Spinor BSE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When spin-orbit coupling is enabled (``soc: true``), the single-particle spatial orbitals are expanded into 2-component Kramers spinors :math:`\psi_p(\mathbf{r}) = \begin{pmatrix} \phi_{p\alpha}(\mathbf{r}) \\ \phi_{p\beta}(\mathbf{r}) \end{pmatrix}`. The BSE Hamiltonian is transformed into the spinor electron-hole basis :math:`|IA\rangle = a_A^\dagger a_I |0\rangle`:

.. math::

   A_{IA, JB}^{\mathrm{spinor}} = \left( \varepsilon_A^{\mathrm{QP}} - \varepsilon_I^{\mathrm{QP}} \right) \delta_{IJ} \delta_{AB} + K_{IA, JB}^x - K_{IA, JB}^d

which carries the spin–orbit contribution to the exciton fine structure (bright–dark ordering) at the level of the single-particle SOC operator; the accuracy of the splitting is limited by the SOC window, the pseudopotential SOC terms and the bare-exchange representation.

