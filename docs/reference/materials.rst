Materials
=========

Part of :doc:`/reference/index`.

.. important::

   The bulk quasiparticle reference of ``MATERIAL_DB`` is the literature 1.0 QSGW gap **with spin-orbit coupling** and our PBE+SOC gap at the same, experimental lattice: :math:`\Delta_\Sigma = E_g^{\mathrm{QSGW+SOC}}(a_{\mathrm{exp}}) - E_g^{\mathrm{PBE+SOC}}(a_{\mathrm{exp}})`. Index 7 is the spin-free PBE gap at :math:`a_{\mathrm{exp}}`, index 8 the spin-free equivalent :math:`E_g^{\mathrm{PBE}}(a_{\mathrm{exp}}) + \Delta_\Sigma`, so that index 8 − index 7 = :math:`\Delta_\Sigma` for spin-free and SOC dots alike. For CdSe: 0.644 eV PBE, 2.286 eV (:math:`\Delta_\Sigma = +1.642` eV, from QSGW+SO 2.16 eV and PBE+SOC 0.522 eV), :math:`\varepsilon_\infty = 6.2`. Sources, lattices and the geometry correction of PBE-relaxed dots: :doc:`/electronic_structure/bulk_bands`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--material``
* YAML: ``system.material``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)


11. Material Database Reference
-------------------------------

Every energy in the tables below is **spin-free (scalar-relativistic)**. The two bulk columns are the spin-free PBE gap at the experimental lattice (index 7) and the spin-free QSGW-equivalent gap (index 8); their difference :math:`\Delta_{\mathrm{bulk}}` is the self-energy correction :math:`\Delta_\Sigma`, obtained with spin-orbit coupling on both sides (next section). The cluster columns record the corresponding values for the finite anchor cluster calculated with CP2K (DZVP-RI / GTH / PBE + evGW), along with the frontier orbital fractions :math:`f_H` and :math:`f_L`.

Physical Rationale: 1.0 QSGW vs. 0.8 QSGW in Colloidal Quantum Dots (2–6 nm)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In macroscopic 3D bulk semiconductors, standard QSGW (under the random phase approximation, RPA) neglects attractive electron-hole ladder vertex corrections in the polarizability / screened interaction :math:`W`. In infinite bulk solids, these vertex corrections increase macroscopic screening (raising :math:`\varepsilon_\infty` by ~10–20%), which contracts the quasiparticle self-energy :math:`\Sigma = iGW` and shrinks the fundamental band gap. To compensate for this in extended bulk systems without solving the full 4-point vertex function, bulk electronic structure practitioners often apply an empirical scaling of 0.8 (e.g., "0.8 QSGW" or QSGW80: :math:`\Sigma_{\mathrm{eff}} = 0.8\,\Sigma_{\mathrm{QSGW}} + 0.2\,V_{\mathrm{LDA}}`).

However, for colloidal quantum dots and nanocrystals of typical experimental sizes (:math:`D \approx 2 - 6\ \mathrm{nm}`), **quantum confinement fundamentally quenches and suppresses these vertex corrections**:

1. **Discretization of the Particle-in-a-Box Spectrum:** The continuous bulk electron-hole continuum is replaced by widely spaced discrete subbands (:math:`\Delta\varepsilon \gg k_B T`), strongly detuning and suppressing virtual electron-hole polarization pairs.
2. **Truncation of the Excitonic Screening Cloud:** In bulk crystals, ladder corrections require a long-range spatial correlation volume spanning several unit cells (:math:`R \ge a_B`, where :math:`a_B` is the bulk exciton Bohr radius, :math:`3 - 10\ \mathrm{nm}`). In nanocrystals with :math:`R_{\mathrm{QD}} \le a_B`, the dielectric boundary physically truncates this correlation cloud.
3. **Dielectric Confinement and Classical Image Charges:** The strong dielectric mismatch between the semiconductor core (:math:`\varepsilon_{\mathrm{in}} \approx 6 - 10`) and the organic ligand/solvent shell (:math:`\varepsilon_{\mathrm{out}} \approx 2`) generates an intense classical self-polarization (image-charge) barrier. This surface-induced self-energy operates without bulk ladder reduction.

If a 0.8 QSGW bulk scissor were used in a 2–6 nm nanocrystal, the quasiparticle opening would be substantially underestimated (for example, yielding an exciton energy of ~2.2 eV for a 3.0 nm CdSe dot, whereas experimental sizing curves demonstrate 2.5–2.6 eV). Pure **1.0 QSGW** provides the physically appropriate high-energy, vertex-quenched asymptote that matches finite vacuum/cluster evGW anchors and accurately reproduces experimental sizing curves (e.g., Hens & Rodina, *Nano Lett.* 2022).

Spin-Orbit Coupling in the Bulk Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

QDEX solves the ground state on spin-free orbitals and adds spin-orbit coupling as spinors (GTH-SOC operator). The literature QSGW gaps include spin-orbit coupling, so the self-energy correction is taken with SOC on both sides and at the same lattice:

.. math::

   \Delta_\Sigma = E_g^{\mathrm{QSGW+SOC}}(a_{\mathrm{exp}}) - E_g^{\mathrm{PBE+SOC}}(a_{\mathrm{exp}}),\qquad
   E_g^{\mathrm{index\ 8}} = E_g^{\mathrm{PBE}}(a_{\mathrm{exp}}) + \Delta_\Sigma .

The same :math:`\Delta_\Sigma` then corrects spin-free and SOC dots, and the spin-orbit lowering of the gap is the PBE/GTH-SOC one in both the bulk reference and the dots: about :math:`\Delta_{\mathrm{so}}/3` in zinc blende (the :math:`\Gamma_8` valence edge), about :math:`\frac{2}{3}\Delta_{\mathrm{so}}^{\mathrm{CB}}` in the lead halide perovskites (the :math:`j = 1/2` conduction edge at R), and the L-point splitting in PbS and PbSe. This replaces the earlier conversion of literature gaps to spin-free ones with those formulas and experimental :math:`\Delta_{\mathrm{so}}`. Gaps are signed band-edge transitions (negative = inverted band order, e.g. the Hg compounds). Literature values and lattices: Deguchi et al., *Jpn. J. Appl. Phys.* 55, 051201 (2016) (III–V, Zn and Cd chalcogenides); Svane et al., *Phys. Rev. B* 81, 245120 (2010) (PbS, PbSe); Svane et al., *Phys. Rev. B* 84, 205205 (2011) (Hg chalcogenides); Huang and Lambrecht, *Phys. Rev. B* 93, 195211 (2016) (cubic CsPbX\ :sub:`3`, carried from their lattice to :math:`a_{\mathrm{exp}}` with their QSGW deformation potentials). Details and the full table: :doc:`/electronic_structure/bulk_bands`.

Perovskites
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0` (Å)
     - PBE bulk
     - QSGW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - Cs\ :sub:`3`\ Bi\ :sub:`2`\ Br\ :sub:`9`
     - 7.85
     - 3.330
     - 4.350
     - +1.020
     - 3.65
     - 7.35
     - +3.70
     - 0.46
     - 0.54
   * - CsPbCl\ :sub:`3`
     - 7.60
     - 1.771
     - 3.558
     - +1.787
     - 3.24
     - 6.90
     - +3.66
     - 0.50
     - 0.50
   * - CsPbBr\ :sub:`3`
     - 7.91
     - 1.278
     - 2.628
     - +1.350
     - 3.08
     - 6.83
     - +3.75
     - 0.43
     - 0.57
   * - CsPbI\ :sub:`3`
     - 8.37
     - 0.832
     - 2.076
     - +1.244
     - 2.69
     - 6.02
     - +3.33
     - 0.34
     - 0.66
   * - MAPbI\ :sub:`3`
     - 
     - 1.550
     - 2.730
     - +1.180
     - 
     - 
     - 
     - 
     - 
   * - FAPbI\ :sub:`3`
     - 
     - 1.450
     - 2.600
     - +1.150
     - 
     - 
     - 
     - 
     - 


II–VI Semiconductors
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0` (Å)
     - PBE bulk
     - QSGW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - ZnS
     - 4.82
     - 2.051
     - 4.128
     - +2.077
     - 3.45
     - 7.26
     - +3.81
     - 0.43
     - 0.57
   * - ZnSe
     - 4.94
     - 1.252
     - 3.221
     - +1.969
     - 3.29
     - 6.84
     - +3.55
     - 0.40
     - 0.60
   * - ZnTe
     - 5.19
     - 1.225
     - 2.940
     - +1.715
     - 2.98
     - 6.10
     - +3.12
     - 0.39
     - 0.61
   * - CdS
     - 5.15
     - 1.146
     - 2.863
     - +1.717
     - 2.72
     - 6.36
     - +3.65
     - 0.44
     - 0.56
   * - CdSe
     - 5.31
     - 0.644
     - 2.286
     - +1.642
     - 2.63
     - 6.03
     - +3.40
     - 0.41
     - 0.59
   * - CdTe
     - 5.44
     - 0.740
     - 2.261
     - +1.521
     - 2.82
     - 6.08
     - +3.26
     - 0.37
     - 0.63
   * - HgS
     - 5.08
     - -0.421
     - 0.585
     - +1.006
     - 2.23
     - 5.47
     - +3.25
     - 0.38
     - 0.62
   * - HgSe
     - 5.20
     - -0.872
     - -0.018
     - +0.854
     - 2.16
     - 5.16
     - +3.00
     - 0.35
     - 0.65
   * - HgTe
     - 5.41
     - -0.665
     - 0.373
     - +1.038
     - 2.18
     - 4.90
     - +2.72
     - 0.34
     - 0.66


III–V Semiconductors
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0` (Å)
     - PBE bulk
     - QSGW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - AlP
     - 5.18
     - 1.667
     - 2.732
     - +1.065
     - 1.94
     - 5.84
     - +3.90
     - 0.47
     - 0.53
   * - AlAs
     - 5.23
     - 1.522
     - 2.458
     - +0.936
     - 2.10
     - 5.73
     - +3.63
     - 0.45
     - 0.55
   * - AlSb
     - 5.56
     - 1.208
     - 1.804
     - +0.596
     - 1.93
     - 5.17
     - +3.24
     - 0.43
     - 0.57
   * - GaP
     - 5.32
     - 1.632
     - 2.487
     - +0.855
     - 1.98
     - 5.40
     - +3.41
     - 0.48
     - 0.52
   * - GaAs
     - 5.45
     - 0.496
     - 1.891
     - +1.395
     - 1.89
     - 5.03
     - +3.13
     - 0.45
     - 0.55
   * - GaSb
     - 5.72
     - 0.156
     - 1.304
     - +1.148
     - 1.48
     - 4.31
     - +2.83
     - 0.44
     - 0.56
   * - InP
     - 5.63
     - 0.680
     - 1.651
     - +0.971
     - 1.75
     - 4.99
     - +3.24
     - 0.45
     - 0.55
   * - InAs
     - 5.70
     - -0.247
     - 0.788
     - +1.035
     - 1.60
     - 4.61
     - +3.01
     - 0.43
     - 0.57
   * - InSb
     - 5.95
     - -0.153
     - 0.778
     - +0.931
     - 1.56
     - 4.30
     - +2.74
     - 0.42
     - 0.58


IV–VI Semiconductors
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0` (Å)
     - PBE bulk
     - QSGW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - PbS
     - 7.24
     - 0.253
     - 0.683
     - +0.430
     - 2.31
     - 5.76
     - +3.45
     - 0.43
     - 0.57
   * - PbSe
     - 7.39
     - 0.139
     - 0.582
     - +0.443
     - 2.27
     - 5.51
     - +3.24
     - 0.41
     - 0.59


