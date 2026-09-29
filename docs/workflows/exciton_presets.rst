Exciton presets
===============

Part of :doc:`/workflows/index`. Each preset shows only the ``excitations`` and ``soc`` sections; the QP
model and the kernel come from :doc:`qp_presets`.

Absorption spectrum
-------------------

.. code-block:: yaml

   excitations:
     mode: bse
     nhomos: 25
     nlumos: 25
     nroots: 40
     full_diag: false            # Davidson

Band-edge fine structure (SOC)
------------------------------

.. code-block:: yaml

   excitations:
     mode: bse
     nhomos: 25
     nlumos: 25
     nroots: 40

   soc:
     enabled: true
     window: 10.0

Singlet–triplet splitting
-------------------------

.. code-block:: yaml

   excitations:
     triplet: true               # no K^x; compare with the singlet run

Many frames (dynamics)
----------------------

.. code-block:: yaml

   excitations:
     mode: diagonal_bse          # each transition with its own K^x and K^d, no mixing
     nhomos: 30
     nlumos: 30
