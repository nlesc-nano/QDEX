Excitations
===========

YAML section: ``excitations`` (spin–orbit coupling in ``soc``).

A neutral excitation promotes an electron from an occupied orbital i to an empty orbital a. QDEX
writes every excited state in the space of these transitions,

.. math::

   |\Psi_S\rangle = \sum_{i\in\mathrm{occ}}\sum_{a\in\mathrm{virt}} X^S_{ia}\,|ia\rangle ,

and obtains the energies Ω_S from a matrix in that space,

.. math::

   A_{ia,jb} = \Delta E_{ia}\,\delta_{ij}\delta_{ab} + K_{ia,jb},\qquad \mathbf A\,\mathbf X_S = \Omega_S\,\mathbf X_S .

The frameworks (``excitations.mode``) differ in two choices: which orbital energies enter
:math:`\Delta E_{ia}`, and how much of the electron–hole interaction K is kept.

.. figure:: /_static/figures/excitation_frameworks.svg
   :width: 100%
   :alt: excitation frameworks

   The frameworks as the structure of the transition-space matrix.

.. list-table::
   :header-rows: 1
   :widths: 20 24 20 36

   * - ``mode``
     - ΔE_ia
     - K
     - Describes
   * - ``independent_dft``
     - KS energies
     - none
     - KS transitions
   * - ``independent_qp``
     - QP energies
     - none
     - non-interacting quasiparticles (fundamental gap)
   * - ``diagonal_bse``
     - QP energies
     - K\ :sup:`x`, K\ :sup:`d`, diagonal only
     - each transition with its own electron–hole interaction
   * - ``bse``
     - QP energies
     - K\ :sup:`x`, K\ :sup:`d`, full
     - coupled excitons (TDA Bethe–Salpeter equation)
   * - ``sbse``, ``diagonal_sbse``
     - PBE + bulk GW
     - K\ :sup:`x`, K\ :sup:`d` with bulk W
     - the same solvers with ΔW removed from both the energies and the kernel
   * - ``stda``, ``diagonal_stda``
     - DFT energies
     - Grimme's γ\ :sup:`K`, a_x-dependent γ\ :sup:`J`
     - the simplified TDA of Grimme

Reading order:

1. :doc:`independent_dft` and :doc:`independent_qp`: no electron–hole interaction.
2. :doc:`kernel`: the two interaction terms K\ :sup:`x` and K\ :sup:`d`, and which W enters K\ :sup:`d`.
3. :doc:`diagonal_bse` and :doc:`bse`: the frameworks with K\ :sup:`x` and K\ :sup:`d`.
4. :doc:`cancellation`: why ΔW in the QP energies and in K\ :sup:`d` does not cancel completely.
5. :doc:`sbse`: dropping ΔW from both.
6. :doc:`stda`: Grimme's simplified TDA.
7. :doc:`configuration`.

.. toctree::
   :maxdepth: 1

   independent_dft
   independent_qp
   kernel
   diagonal_bse
   bse
   cancellation
   sbse
   stda
   configuration
