Configuration
=============

Part of :doc:`/relativity/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.soc_utils``
* Callable: ``qdex.soc_utils.compute_spinor_subspace``
* CLI: ``--soc_flag, --gth_file``
* YAML: ``soc.enabled``, ``soc.window``

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
   * - ``--gth_functional <name>``
     - ``PBE``
     - Functional of the SOC constants read from the GTH file (entry ``GTH-<name>-q<n>``). Use the functional of the DFT calculation; the constants differ between functionals.
   * - ``--soc_window <float>``
     - largest ``|ewin|`` + 1 eV
     - SOC active space of the fuzzy bands: MOs within this many eV of mid-gap. The spinor diagonalization scales as the cube of the number of MOs; levels inside the plot window converge to a few meV with a 1 eV margin (Cs\ :sub:`324`\ Pb\ :sub:`216`\ Br\ :sub:`756`: at most 7 meV against a 2 eV margin for the 1000 spinors nearest the gap).
   * - ``--soc_bse_window <float>``
     - ``soc_window`` when given, else largest ``|ewin|`` + 2 eV
     - SOC is diagonalized for the MOs within this many eV of mid-gap and the BSE spinors are the window spinors living in the BSE MO window (des Cloizeaux projection), so their energies include the SOC coupling to orbitals outside the BSE window. ``0`` diagonalizes in the BSE window alone.
   * - ``--device <choice>``
     - ``auto``
     - Compute device: ``cpu``, ``cuda``, ``mps``, or ``numpy``.
   * - ``--namd-soc``
     - ``False``
     - Enable SOC specifically for Non-Adiabatic Molecular Dynamics trajectory precomputation.


YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   soc:
     enabled: true
     window: 6.0          # default: largest |ewin| + 1 eV
     bse_window: 7.0      # default: soc.window when set, else largest |ewin| + 2 eV
     gth_functional: PBE

   system:
     mo_file: "CsPbI3_MOs.mbse"
     xyz: "CsPbI3_QD.xyz"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-SR-GTH"
     gth_file: "GTH_SOC_POTENTIALS.txt"
