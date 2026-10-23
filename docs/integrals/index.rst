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
     # mnok_exponent: 2        # MNOK exponent beta
     # mnok_onsite: ip_ea      # on-site value IP - EA = 2 eta (benchmarked against xs)
     # one_centre: true        # xs only: same-atom AO products with exact one-centre integrals
     # one_centre_cutoff: 0    # also exact for atom pairs closer than this (A); 0: on-site only
     # one_centre_tail: dipole # charge-dipole and dipole-dipole of the one-centre densities beyond

.. toctree::
   :maxdepth: 1

   representation
