Materials
=========

Part of :doc:`/reference/index`.

.. important::

   Material entries in ``MATERIAL_DB`` are benchmark model inputs using **pure 1.0 QSGW scalar-relativistic (spin-free)** bulk values. For example, ``MATERIAL_DB["CDSE"]`` uses 0.62 eV PBE, 2.19 eV QSGW (:math:`\Delta_{\mathrm{bulk}} = +1.57\ \mathrm{eV}`), and :math:`\varepsilon_\infty = 6.2`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--material``
* YAML: ``system.material, physics.material``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)


11. Material Database Reference
-------------------------------

Every energy in the tables below is **spin-free (scalar-relativistic)**. The two bulk columns are a scalar PBE fundamental gap and a spin-free 1.0 QSGW gap. Their difference is :math:`\Delta_{\mathrm{bulk}} = E_g^{\mathrm{bulk, QSGW}} - E_g^{\mathrm{bulk, PBE}}`. The cluster columns record the corresponding values for the finite anchor cluster calculated with CP2K (DZVP-RI / GTH / PBE + evGW), along with the frontier orbital fractions :math:`f_H` and :math:`f_L`.

Physical Rationale: 1.0 QSGW vs. 0.8 QSGW in Colloidal Quantum Dots (2–6 nm)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In macroscopic 3D bulk semiconductors, standard QSGW (under the random phase approximation, RPA) neglects attractive electron-hole ladder vertex corrections in the polarizability / screened interaction :math:`W`. In infinite bulk solids, these vertex corrections increase macroscopic screening (raising :math:`\varepsilon_\infty` by ~10–20%), which contracts the quasiparticle self-energy :math:`\Sigma = iGW` and shrinks the fundamental band gap. To compensate for this in extended bulk systems without solving the full 4-point vertex function, bulk electronic structure practitioners often apply an empirical scaling of 0.8 (e.g., "0.8 QSGW" or QSGW80: :math:`\Sigma_{\mathrm{eff}} = 0.8\,\Sigma_{\mathrm{QSGW}} + 0.2\,V_{\mathrm{LDA}}`).

However, for colloidal quantum dots and nanocrystals of typical experimental sizes (:math:`D \approx 2 - 6\ \mathrm{nm}`), **quantum confinement fundamentally quenches and suppresses these vertex corrections**:

1. **Discretization of the Particle-in-a-Box Spectrum:** The continuous bulk electron-hole continuum is replaced by widely spaced discrete subbands (:math:`\Delta\varepsilon \gg k_B T`), strongly detuning and suppressing virtual electron-hole polarization pairs.
2. **Truncation of the Excitonic Screening Cloud:** In bulk crystals, ladder corrections require a long-range spatial correlation volume spanning several unit cells (:math:`R \ge a_B`, where :math:`a_B` is the bulk exciton Bohr radius, :math:`3 - 10\ \mathrm{nm}`). In nanocrystals with :math:`R_{\mathrm{QD}} \le a_B`, the dielectric boundary physically truncates this correlation cloud.
3. **Dielectric Confinement and Classical Image Charges:** The strong dielectric mismatch between the semiconductor core (:math:`\varepsilon_{\mathrm{in}} \approx 6 - 10`) and the organic ligand/solvent shell (:math:`\varepsilon_{\mathrm{out}} \approx 2`) generates an intense classical self-polarization (image-charge) barrier. This surface-induced self-energy operates without bulk ladder reduction.

If a 0.8 QSGW bulk scissor were used in a 2–6 nm nanocrystal, the quasiparticle opening would be substantially underestimated (for example, yielding an exciton energy of ~2.2 eV for a 3.0 nm CdSe dot, whereas experimental sizing curves demonstrate 2.5–2.6 eV). Pure **1.0 QSGW** provides the physically appropriate high-energy, vertex-quenched asymptote that matches finite vacuum/cluster evGW anchors and accurately reproduces experimental sizing curves (e.g., Hens & Rodina, *Nano Lett.* 2022).

Relativistic Spin-Orbit Coupling (SOC) Correction Conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because QDEX solves the Kohn-Sham ground state and Bethe-Salpeter equation on spin-free spatial orbitals (with optional subsequent SOC spinor transformations), all reference bulk gaps in ``MATERIAL_DB`` are strictly **scalar-relativistic (spin-free)**:

* **Zinc-Blende II–VI and III–V Semiconductors:** Spin-orbit coupling splits the triply degenerate :math:`\Gamma_{15v}` valence band into a 4-fold :math:`\Gamma_{8v}` (heavy/light hole) and a 2-fold :math:`\Gamma_{7v}` (split-off) state. The valence band maximum is pushed upward by :math:`+\frac{1}{3}\Delta_{\mathrm{so}}`, while the :math:`s`-like conduction band minimum (:math:`\Gamma_{1c}`) undergoes negligible first-order shift. Consequently:

  .. math::

     E_g^{\mathrm{spin-free}} = E_g^{\mathrm{with\ SOC}} + \frac{1}{3}\Delta_{\mathrm{so}}

* **Lead Halide Perovskites (Cubic Phase):** The band edges at :math:`R` feature an :math:`s`-like valence band maximum (Pb :math:`6s` - halide :math:`np` antibonding, negligible SOC shift) and a :math:`p`-like conduction band minimum (Pb :math:`6p`). SOC splits the Pb :math:`6p` manifold into :math:`j=1/2` and :math:`j=3/2`, lowering the :math:`j=1/2` conduction edge by :math:`\approx \frac{2}{3}\Delta_{\mathrm{so}}(\mathrm{Pb}\ 6p) \approx 0.95 - 1.05\ \mathrm{eV}`. Consequently:

  .. math::

     E_g^{\mathrm{spin-free}} = E_g^{\mathrm{with\ SOC}} + \frac{2}{3}\Delta_{\mathrm{so}}^{\mathrm{CB}}


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
     - 3.33
     - 4.35
     - +1.02
     - 3.65
     - 7.35
     - +3.70
     - 0.46
     - 0.54
   * - CsPbCl\ :sub:`3`
     - 7.60
     - 2.45
     - 4.20
     - +1.75
     - 3.24
     - 6.90
     - +3.66
     - 0.50
     - 0.50
   * - CsPbBr\ :sub:`3`
     - 7.91
     - 1.85
     - 3.40
     - +1.55
     - 3.08
     - 6.83
     - +3.75
     - 0.43
     - 0.57
   * - CsPbI\ :sub:`3`
     - 8.37
     - 1.55
     - 2.80
     - +1.25
     - 2.69
     - 6.02
     - +3.33
     - 0.34
     - 0.66
   * - MAPbI\ :sub:`3`
     - 
     - 1.55
     - 2.73
     - +1.18
     - 
     - 
     - 
     - 
     - 
   * - FAPbI\ :sub:`3`
     - 
     - 1.45
     - 2.60
     - +1.15
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
     - 2.08
     - 4.17
     - +2.09
     - 3.45
     - 7.26
     - +3.81
     - 0.43
     - 0.57
   * - ZnSe
     - 4.94
     - 1.27
     - 3.26
     - +1.99
     - 3.29
     - 6.84
     - +3.55
     - 0.40
     - 0.60
   * - ZnTe
     - 5.19
     - 1.45
     - 2.96
     - +1.51
     - 2.98
     - 6.10
     - +3.12
     - 0.39
     - 0.61
   * - CdS
     - 5.15
     - 1.12
     - 2.78
     - +1.66
     - 2.72
     - 6.36
     - +3.65
     - 0.44
     - 0.56
   * - CdSe
     - 5.31
     - 0.62
     - 2.19
     - +1.57
     - 2.63
     - 6.03
     - +3.40
     - 0.41
     - 0.59
   * - CdTe
     - 5.44
     - 0.86
     - 2.12
     - +1.26
     - 2.82
     - 6.08
     - +3.26
     - 0.37
     - 0.63
   * - HgS
     - 5.08
     - −0.15
     - 0.31
     - +0.46
     - 2.23
     - 5.47
     - +3.25
     - 0.38
     - 0.62
   * - HgSe
     - 5.20
     - −0.40
     - 0.05
     - +0.45
     - 2.16
     - 5.16
     - +3.00
     - 0.35
     - 0.65
   * - HgTe
     - 5.41
     - −0.30
     - 0.23
     - +0.53
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
     - 1.59
     - 2.74
     - +1.15
     - 1.94
     - 5.84
     - +3.90
     - 0.47
     - 0.53
   * - AlAs
     - 5.23
     - 1.39
     - 2.46
     - +1.07
     - 2.10
     - 5.73
     - +3.63
     - 0.45
     - 0.55
   * - AlSb
     - 5.56
     - 1.18
     - 1.80
     - +0.62
     - 1.93
     - 5.17
     - +3.24
     - 0.43
     - 0.57
   * - GaP
     - 5.32
     - 1.61
     - 2.49
     - +0.88
     - 1.98
     - 5.40
     - +3.41
     - 0.48
     - 0.52
   * - GaAs
     - 5.45
     - 0.49
     - 1.89
     - +1.40
     - 1.89
     - 5.03
     - +3.13
     - 0.45
     - 0.55
   * - GaSb
     - 5.72
     - 0.11
     - 1.20
     - +1.09
     - 1.48
     - 4.31
     - +2.83
     - 0.44
     - 0.56
   * - InP
     - 5.63
     - 0.46
     - 1.65
     - +1.19
     - 1.75
     - 4.99
     - +3.24
     - 0.45
     - 0.55
   * - InAs
     - 5.70
     - −0.42
     - 0.80
     - +1.22
     - 1.60
     - 4.61
     - +3.01
     - 0.43
     - 0.57
   * - InSb
     - 5.95
     - −0.61
     - 0.77
     - +1.38
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
     - 0.15
     - 0.73
     - +0.58
     - 2.31
     - 5.76
     - +3.45
     - 0.43
     - 0.57
   * - PbSe
     - 7.39
     - 0.05
     - 0.65
     - +0.60
     - 2.27
     - 5.51
     - +3.24
     - 0.41
     - 0.59


