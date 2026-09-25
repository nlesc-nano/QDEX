Configuration
=============

Part of :doc:`/relativity/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.soc_utils``
* Callable: ``qdex.soc_utils.compute_spinor_subspace``
* CLI: ``--soc_flag, --gth_file``
* YAML: ``physics.soc, physics.soc_window_ev``

.. code-block:: python

   compute_spinor_subspace(atom_symbols, coords_ang, shells, C_AO, eps_Ha, S_AO, active_indices, gth_file, nthreads=1, soc_cache=None, assume_orthonormal=False, SC_AO=None, device='numpy', verbose=True)


7. CLI Flags & YAML Configuration Reference
-------------------------------------------


Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--soc_flag``
     - ``False``
     - Enable relativistic Spin-Orbit Coupling across the calculation.
   * - ``--gth_file <path>``
     - Built-in table
     - Path to custom GTH SOC pseudopotential parameter file (overriding internal database).
   * - ``--soc_window <float>``
     - ``10.0``
     - Energy window (in eV) centered at the Fermi level for selecting MOs in the SOC active space.
   * - ``--device <choice>``
     - ``auto``
     - Compute device: ``cpu``, ``cuda``, ``mps``, or ``numpy``.
   * - ``--namd-soc``
     - ``False``
     - Enable SOC specifically for Non-Adiabatic Molecular Dynamics trajectory precomputation.


YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   physics:
     soc: true
     soc_window_ev: 8.0
     gth_potentials: "GTH_SOC_POTENTIALS.txt"

   system:
     mo_file: "CsPbI3_MOs.mbse"
     xyz: "CsPbI3_QD.xyz"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-SR-GTH"
