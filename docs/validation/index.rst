Validation
==========

The QP and excitation models tested on CdSe dots against experiment and against an explicit GW
calculation.

* :doc:`cdse_experiment`: first-exciton energies at 1.2 and 2.0 nm, the experimental window, and how
  successive versions of the QP layer moved S₁ toward it.
* :doc:`evgw_cluster`: QP gaps and band-edge shifts of Cd₁₆Se₁₃Cl₆ against evGW\@PBE0.
* :doc:`qp_bse_sweep`: the script that reproduces these tables for any dot.
* :doc:`bulk_exciton_limit`: what is required for the bulk (Wannier–Mott) limit.

.. toctree::
   :maxdepth: 1

   cdse_experiment
   evgw_cluster
   qp_bse_sweep
   bulk_exciton_limit
