Two-electron integrals
======================

YAML section: ``integrals``.

The QP correction and the excited states both need two-electron integrals over orbital-pair
densities: of ΔW for the self-energy, of W for the direct term K\ :sup:`d` and of v for the exchange
term K\ :sup:`x`. The ``integrals`` section chooses, once for all of them, how these integrals are
represented on the atomic basis.

.. code-block:: yaml

   integrals:
     representation: mnok      # mnok (atom-condensed) or xs (ZDO AO density pairs)
     charges: mulliken         # MNOK population partition: mulliken or lowdin (xs uses Löwdin)

.. toctree::
   :maxdepth: 1

   representation
