Configuration
=============

Part of :doc:`/excitons/index`. The excited states are set in the ``excitations`` section; spin–orbit
coupling in ``soc``. All keys are listed in :doc:`/getting_started/configuration`.

.. code-block:: yaml

   excitations:
     mode: bse                 # independent_dft, independent_qp, diagonal_bse, bse (sbse, diagonal_sbse)
     # kernel: set by the QP model (qp for Resta/DIM); choose only for models without W
     nhomos: 25
     nlumos: 25
     nroots: 40
     full_diag: false          # Davidson
     triplet: false

   soc:
     enabled: true
     window: 10.0

The sBSE (bulk W, no finite-size QP correction; :doc:`screened_kernel`):

.. code-block:: yaml

   quasiparticles:
     model: none               # KS energies + bulk GW correction

   excitations:
     mode: sbse
     kernel: resta             # or dim
