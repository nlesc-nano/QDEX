Qp edges
========

Part of :doc:`/workflows/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system``, ``environment``, ``quasiparticles``, ``integrals``, ``excitations``, ``namd``

.. code-block:: python

   main()


Tutorial 3: Scaled GW Band Gap Correction & Absolute Edges
----------------------------------------------------------

Semi-local DFT underestimates the band gap. Here we activate the Scaled GW model with dielectric solvent screening (:math:`\epsilon_{\mathrm{out}} = 2.25` for toluene) and predict absolute IP and EA levels.


Configuration File (``tutorial3_gw.yaml``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   system:
     mo_file: "CsPbBr3_MOs.mbse"
     xyz: "CsPbBr3_QD.xyz"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-PBE-GTH"
     cif: "CsPbBr3_bulk.cif"

     material: "CSPBBR3"

   environment:
     eps_out: 2.25

   quasiparticles:
     model: sgw-resta
     energy_reference: vacuum    # absolute IP/EA

   analysis:
     run_fuzzy: true
     dashboard_energy_mode: both
     pdos_atoms: ["Pb", "Br"]


Execution
~~~~~~~~~

.. code-block:: bash

   qdex --config tutorial3_gw.yaml


Key Output
~~~~~~~~~~

The terminal prints the absolute band edges of the dot:

.. code-block:: text

   [Absolute Band Edges (IP & EA)]
     PBE HOMO / LUMO  :  e_H / e_L eV (CP2K eigenvalues; on the vacuum scale only for PERIODIC NONE ...)
     QP shifts        : HOMO d_H eV, LUMO d_L eV (orbital-resolved QP levels)
     Bulk edge split  : f_b = 0.500 of the bulk QP shift on the HOMO (symmetric split ...)
     QP HOMO (IP)     :  e_H + d_H eV   -> IP = -(e_H + d_H) eV
     QP LUMO (EA)     :  e_L + d_L eV   -> EA = -(e_L + d_L) eV
     Split range      : f_b = 0 -> 1 gives IP ... eV, EA ... eV (bulk shift ... eV; the gap does not depend on f_b)

The edges are the dot's own PBE eigenvalues plus the QP shift of each edge:

.. math::

   \mathrm{IP} = -\left(\varepsilon_\mathrm{HOMO}^\mathrm{PBE} + \Delta_\mathrm{HOMO}\right), \qquad
   \mathrm{EA} = -\left(\varepsilon_\mathrm{LUMO}^\mathrm{PBE} + \Delta_\mathrm{LUMO}\right).

A finite cluster computed with ``PERIODIC NONE`` and an isolated Poisson solver has its eigenvalues on
the vacuum scale (``MULTIPOLE`` is within about 0.04 eV of ``WAVELET`` and ``MT`` for the 1.2 nm CdSe
cluster in a 28 Å box), so dots of different sizes compare directly; keep a similar amount of vacuum
around each dot. A run with a periodic Poisson solver has an arbitrary eigenvalue zero, and its IP/EA are
meaningless (the 1.2 nm MO file in ``tests/CdSe/1.2nm`` lies rigidly 0.66 eV above the same cluster
computed with ``PERIODIC NONE`` and ``MULTIPOLE``).
With orbital-resolved QP levels (the Delta-W models, ``quasiparticles.levels: orbital``)
:math:`\Delta_\mathrm{HOMO}` and :math:`\Delta_\mathrm{LUMO}` are the model's own shifts of the two
orbitals. ``bulk`` adds the ``sgw-resta`` self-energy of the two edges to IP and EA only
(``quasiparticles.ip_ea: resta``, :doc:`/quasiparticles/models`), so that it gives the same edges as
``sgw-resta``; with rigid levels otherwise the gap correction is split between the edges with the model's
HOMO fraction. Periodic runs have no vacuum level, and no IP/EA is assigned.

Both ``fuzzy_dashboard_dft.html`` and ``fuzzy_dashboard_qp.html`` are created, allowing direct side-by-side comparison of DFT and quasiparticle band structures.


Splitting the bulk correction between the edges
-----------------------------------------------

The bulk QP correction :math:`\Delta_{\mathrm{bulk}}` (bulk QSGW gap − bulk PBE gap, with the vertex settings)
is a gap correction. Its HOMO share :math:`f_b` (HOMO :math:`-f_b\Delta_{\mathrm{bulk}}`, LUMO
:math:`+(1-f_b)\Delta_{\mathrm{bulk}}`) shifts every level by the same constant: it moves IP and EA together by
:math:`\Delta_{\mathrm{bulk}}` per unit of :math:`f_b` and leaves the gap, the levels relative to each other
and the BSE unchanged. ``quasiparticles.bulk_edge_split`` sets it:

* ``symmetric`` (default): :math:`f_b = 1/2` for every material, the rigid :math:`\pm\Delta/2` split of
  Biffi, Cho, Krahne and Berkelbach (arXiv:2210.01324);
* ``cluster``: per material, the value with which QDEX reproduces the evGW@PBE0 IP of the
  M\ :sub:`16`\ X\ :sub:`13`\ Cl\ :sub:`6`-type cluster anchors (``CLUSTER_EDGE_SPLIT`` in ``qdex/hardness.py``);
* ``bulk``: bulk band-edge shifts from GW with vertex corrections (Grüneis et al., *PRL* 112, 096401 (2014))
  and experiment (``BULK_EDGE_SPLIT``);
* a number.

The log prints IP and EA for :math:`f_b = 0` and :math:`f_b = 1` as the range of the split.

**Why the symmetric default.** :math:`f_b` is not known to better than about ±0.5 for nanocrystals:

* Bulk crystals: 75–98 % of the opening lowers the valence band of II-VI and III-V semiconductors (GWΓ and
  experiment vs PBE slab IPs; CdSe 0.94–0.95).
* Clusters of 1–2 nm (CP2K, DZVP-MOLOPT, ``PERIODIC NONE``/``MULTIPOLE``): after QDEX's finite-size term the
  evGW@PBE0 correction left over is almost entirely on the LUMO for the Zn, Cd, Ga and In compounds
  (:math:`f_b \approx 0`; CdSe 0.07 at 1.2 nm, 0.06 at 1.9 nm) and on both edges for Al pnictides, Pb
  chalcogenides and CsPbX\ :sub:`3` (0.3–0.9). The same split appears without QDEX: IP from evGW@PBE0 and from
  ΔSCF-PBE agree within 0.1 eV, the EA differ by about 0.8 eV.
* The split of cluster GW is fragile: G0W0@PBE, evGW@PBE and evGW@PBE0 spread by 0.5 eV in IP for the 1.2 nm
  CdSe cluster, and TZV2P lowers the LUMO by 0.25–0.30 eV against DZVP. In silicon nanocrystals, stochastic
  GW (Neuhauser et al., *PRL* 113, 076402 (2014)) and Tiago and Chelikowsky (*PRB* 73, 205334 (2006)) find the
  bulk-like part on the IP, while G0W0 with WEST (Govoni and Galli, *JCTC* 11, 2680 (2015)) does not.
* Measured IPs of CdSe nanocrystals (5.3–5.9 eV, Jasieniak et al., *ACS Nano* 5, 5888 (2011); Ehamparam et al.,
  *ACS Nano* 9, 8786 (2015)) lie about 1 eV below the bulk value even for unconfined dots: ligand and surface
  dipoles (0.3–2 eV) dominate the absolute levels, and experiment cannot fix :math:`f_b`.

Vertical IP and EA of the two CdSe clusters (vacuum, eV; QDEX with ``sgw-resta``, ``bulk_vertex: scaled``,
``bulk_vertex_factor: material``, ``bulk_residual: experimental``; ``bulk`` with ``ip_ea: resta`` gives the same):

.. list-table::
   :header-rows: 1

   * - Method
     - Cd\ :sub:`16`\ Se\ :sub:`13`\ Cl\ :sub:`6` IP / EA / gap
     - Cd\ :sub:`68`\ Se\ :sub:`55`\ Cl\ :sub:`26` IP / EA / gap
   * - PBE eigenvalues
     - 6.42 / 3.79 / 2.64
     - 5.76 / 4.28 / 1.48
   * - ΔSCF-PBE
     - 7.74 / 2.58 / 5.16
     - 6.60 / 3.51 / 3.09
   * - G0W0@PBE (DZVP)
     - 7.33 / 2.17 / 5.16
     - 6.30 / 3.07 / 3.22
   * - evGW0@PBE (TZV2P)
     - 7.60 / 2.27 / 5.33
     - —
   * - evGW0@PBE0 (ADMM)
     - 7.70 / 1.83 / 5.87
     - 6.63 / 2.75 / 3.88
   * - evGW@PBE0
     - 7.83 / 1.79 / 6.05
     - —
   * - QDEX, :math:`f_b = 1/2`
     - 8.41 / 1.86 / 6.55
     - 7.18 / 2.88 / 4.30
   * - QDEX, ``cluster`` (0.07)
     - 7.83 / 1.27 / 6.55
     - 6.64 / 2.34 / 4.30
   * - QDEX, ``bulk`` (0.95)
     - 9.02 / 2.47 / 6.55
     - 7.75 / 3.45 / 4.30

Over the 23 cluster anchors, the mean absolute errors of QDEX (same settings) against evGW@PBE0 are
0.46 / 0.29 eV (IP / EA) with :math:`f_b = 1/2`, 0.29 / 0.49 eV with :math:`f_b = 0` and 0.78 / 0.53 eV with
the bulk values; the gap error (+0.12 eV on average; +0.6 for the Zn, Cd and Hg chalcogenides, −0.6 to −0.9
for Al pnictides, Pb chalcogenides and CsPbX\ :sub:`3`) goes to whichever edge :math:`f_b` leaves it. In the
bulk limit :math:`f_b = 1/2` places both CdSe edges about 0.5 eV above the measured ones, within the same
uncertainty.
