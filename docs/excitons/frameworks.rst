Frameworks
==========

Part of :doc:`/excitons/index`.

.. figure:: /_static/figures/excitation_frameworks.svg
   :width: 100%
   :alt: excitation frameworks

   The four ``excitation_mode`` choices, shown as the structure of the transition-space matrix.

.. rubric:: QDEX implementation

* Module: ``qdex.exciton_hamiltonian``, ``qdex.solver``
* CLI: ``--excitation-mode``, ``--include-direct-eh``, ``--include-exchange``
* YAML: ``physics.excitation_mode``, ``physics.include_direct_eh``, ``physics.include_exchange``

The four frameworks differ in which energies they use and which parts of K\ :sup:`x` and K\ :sup:`d`
they keep. In all of them K\ :sup:`x` is built with the bare v and K\ :sup:`d` with the W of the QP
model (:doc:`screened_kernel`).

.. list-table::
   :header-rows: 1
   :widths: 18 26 18 38

   * - ``excitation_mode``
     - Energy
     - K\ :sup:`x`, K\ :sup:`d`
     - What it describes
   * - ``independent_dft``
     - :math:`\varepsilon_a - \varepsilon_i`
     - none
     - KS transitions; no QP correction, no electron–hole interaction.
   * - ``independent_qp``
     - :math:`\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}`
     - none
     - Non-interacting quasiparticles. The QP gap contains the full ΔW self-image, and nothing
       compensates it: the energies are far too high and strongly solvent dependent.
   * - ``diagonal_bse``
     - :math:`\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} + 2K^x_{ia,ia} - K^d_{ia,ia}`
     - diagonal only
     - Each transition with its own exchange and its own electron–hole attraction, no mixing. The
       mutual image of each pair compensates its QP self-image, so most of the ΔW cancellation
       (:doc:`cancellation`) is already there.
   * - ``bse``
     - eigenvalues of :math:`A_{ia,jb}`
     - full
     - Coupled TDA BSE: transitions mix, the electron and hole correlate, and the oscillator
       strength redistributes.

**K**\ :sup:`x` **in each framework.** The diagonal element :math:`K^x_{ia,ia}` is the self-interaction
of the transition density. It shifts singlets up relative to triplets. The off-diagonal elements
couple transitions with large transition densities and push oscillator strength to the bright state.
Triplets (``--triplet``) have no K\ :sup:`x`. With SOC the spinor BSE has :math:`K^x - K^d` with no
factor 2.

**K**\ :sup:`d` **in each framework.** :math:`K^d_{ia,ia}` is the attraction of the electron density
:math:`|\psi_a|^2` to the hole density :math:`|\psi_i|^2` through W. With a shared W it contains the
mutual image. The off-diagonal :math:`K^d_{ia,jb}` mixes transitions and lowers S₁ further through
correlation.

Why diagonal BSE works in nanocrystals
--------------------------------------

The coupled BSE costs :math:`\mathcal O(N_{\mathrm{pairs}})` memory per root with Davidson and
:math:`\mathcal O(N_{\mathrm{pairs}}^3)` with full diagonalization. The diagonal BSE costs
:math:`\mathcal O(N_{\mathrm{pairs}})`:

.. math::

   \Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} + 2K^x_{ia,ia} - K^d_{ia,ia}.

* **Most of the binding is diagonal.** In a strongly confined dot the band-edge exciton is dominated
  by one or a few transitions. The share of the binding recovered is system dependent: 58 % for the
  supplied CdSe example.
* **Dynamics.** In non-adiabatic molecular dynamics the excited states are needed at every step.
  Diagonal BSE gives an approximate surface at low cost; check its error against the coupled BSE for
  the chosen active space.
