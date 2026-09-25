Configuration
=============

Part of :doc:`/electronic_structure/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.pdos_coop``
* Callable: ``qdex.pdos_coop.compute_pdos_and_coop``
* CLI: ``--charge_type, --run_fuzzy``
* YAML: ``physics.charge_type, fuzzy.run``

.. code-block:: python

   compute_pdos_and_coop(C, S, eps_eV, shells, pdos_atoms, coop_pairs, ewin, sigma=0.03, is_soc=False, prefix='sf', pops=None, population_bars=None, device='numpy')


7. CLI Flags & YAML Configuration Reference
-------------------------------------------

All ground-state electronic structure analyses can be triggered from the command-line interface or configured within a YAML file.


Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--run_fuzzy``
     - ``False``
     - Enable Fuzzy Band structure and PDOS/COOP analysis workflow.
   * - ``--cif <path>``
     - ``None``
     - Reference bulk crystal structure CIF file used to generate the Brillouin zone :math:`k`-path.
   * - ``--pdos_atoms <list>``
     - ``None``
     - List of atomic element symbols to project PDOS upon (e.g. ``--pdos_atoms Pb Br Cs``).
   * - ``--coop_pairs <list>``
     - ``None``
     - List of bonded pairs for COOP analysis (e.g. ``--coop_pairs Pb-Br Cs-Br``).
   * - ``--ewin <min> <max>``
     - ``-5.0 5.0``
     - Energy window (in eV) relative to the Fermi level for PDOS, COOP, and Fuzzy Band export.
   * - ``--pdos_sigma <float>``
     - ``0.10``
     - Gaussian smearing standard deviation :math:`\sigma_{\mathrm{pdos}}` (in eV) for continuous PDOS.
   * - ``--fuzzy_sigma <float>``
     - ``0.03``
     - Gaussian smearing standard deviation (in eV) along the energy axis for fuzzy band plots.
   * - ``--fold_to_bz``
     - ``False``
     - Fold plane-wave projections into the first Brillouin zone using reciprocal replicas.
   * - ``--g_shell <int>``
     - ``0``
     - Reciprocal lattice vector shell for BZ folding (0 = 1 replica, 1 = 27 replicas, 2 = 125 replicas).
   * - ``--cube``
     - ``False``
     - Export 3D volumetric Gaussian ``.cube`` files for frontier orbitals.
   * - ``--cube-spacing <float>``
     - ``0.5``
     - Grid resolution spacing in Angstroms for ``.cube`` files.
   * - ``--cube-nhomos <int>``
     - ``2``
     - Number of occupied frontier orbitals to export (e.g. HOMO, HOMO-1).
   * - ``--cube-nlumos <int>``
     - ``2``
     - Number of unoccupied frontier orbitals to export (e.g. LUMO, LUMO+1).
   * - ``--disable_cpp_cube``
     - ``False``
     - Disable the C++ Libint grid evaluator and fall back to pure Python grid generation.


YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   system:
     mo_file: "CsPbBr3_MOs.mbse"
     xyz: "CsPbBr3_QD.xyz"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-PBE-GTH"
     cif: "CsPbBr3_bulk.cif"

   fuzzy:
     run: true
     pdos_atoms: ["Pb", "Br", "Cs"]
     coop_pairs: ["Pb-Br"]
     ewin: [-4.0, 4.0]
     pdos_sigma: 0.08
     fuzzy_sigma: 0.03
     fold_to_bz: true
     g_shell: 1

   cube:
     export: true
     spacing_ang: 0.4
     nhomos: 3
     nlumos: 3
