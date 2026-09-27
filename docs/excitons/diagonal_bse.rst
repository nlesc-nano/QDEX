Diagonal BSE (``diagonal_bse``)
===============================

Part of :doc:`/excitons/index`.

.. code-block:: yaml

   excitations:
     mode: diagonal_bse

**Energies.** Each transition keeps its own exchange and its own electron–hole attraction; the
off-diagonal elements of K are dropped:

.. math::

   \Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} + 2K^x_{ia,ia} - K^d_{ia,ia}
   \qquad(\text{triplet: } -K^d_{ia,ia}).

Intensities are those of the bare transition, :math:`f_{ia} \propto \Omega_{ia}|\boldsymbol\mu_{ia}|^2`.

**What it keeps.**

* The attraction :math:`K^d_{ia,ia}` of the electron density :math:`|\psi_a|^2` to the hole density
  :math:`|\psi_i|^2`. With a shared W it contains the mutual image, which compensates the self-image of
  the same transition in the QP energies. Most of the ΔW cancellation (:doc:`cancellation`) is
  therefore already present.
* The exchange :math:`K^x_{ia,ia}`, which shifts singlets above triplets.

**What it misses.**

* **Mixing.** Transitions do not combine, so the electron and hole cannot correlate their positions.
  The share of the coupled binding recovered is system dependent: 58 % for the supplied CdSe example.
* **Intensity redistribution.** Oscillator strength stays with each transition; it is not collected
  into the bright band-edge state.
* **Fine structure** from the coupling of near-degenerate transitions.

**Why use it.** It costs :math:`\mathcal O(N_{\mathrm{pairs}})`: only the diagonal of A is built and no
eigenproblem is solved. In non-adiabatic dynamics, where excited states are needed at every step of
thousands, it gives an approximate surface at low cost. Check its error against :doc:`bse` for the
chosen active space.
