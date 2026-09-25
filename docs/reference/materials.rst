Materials
=========

Part of :doc:`/reference/index`.

.. important::

   Material entries are model inputs, with scalar-relativistic PBE and GW gap conventions where given. Current ``MATERIAL_DB["CDSE"]`` uses 0.64 eV PBE, 1.91 eV GW, and epsilon_inf=6.2; alternative CdSe numbers must be separate, cited data sets.

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

Every energy in the tables below is spin-free. The two bulk columns are a scalar PBE fundamental gap and a spin-free GW gap. Their difference is :math:`\Delta_{\mathrm{bulk}}`. The cluster columns are the same difference on the finite anchor cluster, plus the separate HOMO and LUMO shifts that fix :math:`f_{\mathrm{homo}}` and :math:`f_{\mathrm{lumo}}`.


Perovskites
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 8 8 8 8 8 8 8 6 6

   * - Material
     - :math:`R_0` (Å)
     - PBE bulk
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - Cs\ :sub:`3`\ Bi\ :sub:`2`\ Br\ :sub:`9`
     - 7.85
     - 3.33
     - 3.65
     - +0.32
     - 3.65
     - 7.35
     - +3.70
     - 0.46
     - 0.54
   * - CsPbCl\ :sub:`3`
     - 7.60
     - 2.45
     - 3.65
     - +1.20
     - 3.24
     - 6.90
     - +3.66
     - 0.50
     - 0.50
   * - CsPbBr\ :sub:`3`
     - 7.91
     - 1.85
     - 2.85
     - +1.00
     - 3.08
     - 6.83
     - +3.75
     - 0.43
     - 0.57
   * - CsPbI\ :sub:`3`
     - 8.37
     - 1.55
     - 2.45
     - +0.90
     - 2.69
     - 6.02
     - +3.33
     - 0.34
     - 0.66
   * - MAPbI\ :sub:`3`
     - 
     - 1.55
     - 2.50
     - +0.95
     - 
     - 
     - 
     - 
     - 
   * - FAPbI\ :sub:`3`
     - 
     - 1.45
     - 2.35
     - +0.90
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
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - ZnS
     - 4.82
     - 2.11
     - 3.73
     - +1.62
     - 3.45
     - 7.26
     - +3.81
     - 0.43
     - 0.57
   * - ZnSe
     - 4.94
     - 1.25
     - 2.76
     - +1.51
     - 3.29
     - 6.84
     - +3.55
     - 0.40
     - 0.60
   * - ZnTe
     - 5.19
     - 1.18
     - 2.35
     - +1.17
     - 2.98
     - 6.10
     - +3.12
     - 0.39
     - 0.61
   * - CdS
     - 5.15
     - 1.14
     - 2.55
     - +1.41
     - 2.72
     - 6.36
     - +3.65
     - 0.44
     - 0.56
   * - CdSe
     - 5.26
     - 0.64
     - 1.91
     - +1.27
     - 2.63
     - 6.03
     - +3.39
     - 0.41
     - 0.59
   * - CdTe
     - 5.44
     - 0.61
     - 1.62
     - +1.01
     - 2.82
     - 6.08
     - +3.26
     - 0.37
     - 0.63
   * - HgS
     - 5.08
     - 0.00
     - 0.50
     - +0.50
     - 2.23
     - 5.47
     - +3.25
     - 0.38
     - 0.62
   * - HgSe
     - 5.20
     - −0.42
     - 0.12
     - +0.54
     - 2.16
     - 5.16
     - +3.00
     - 0.35
     - 0.65
   * - HgTe
     - 5.41
     - −0.61
     - 0.10
     - +0.71
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
     - GW bulk
     - :math:`\Delta_{\mathrm{bulk}}`
     - PBE clus
     - GW clus
     - :math:`\Delta(R_0)`
     - :math:`f_H`
     - :math:`f_L`
   * - AlP
     - 5.18
     - 1.59
     - 2.44
     - +0.85
     - 1.94
     - 5.84
     - +3.90
     - 0.47
     - 0.53
   * - AlAs
     - 5.23
     - 1.39
     - 2.15
     - +0.76
     - 2.10
     - 5.73
     - +3.63
     - 0.45
     - 0.55
   * - AlSb
     - 5.56
     - 1.18
     - 1.83
     - +0.65
     - 1.93
     - 5.17
     - +3.24
     - 0.43
     - 0.57
   * - GaP
     - 5.32
     - 1.61
     - 2.36
     - +0.75
     - 1.98
     - 5.40
     - +3.41
     - 0.48
     - 0.52
   * - GaAs
     - 5.45
     - 0.49
     - 1.42
     - +0.93
     - 1.89
     - 5.03
     - +3.13
     - 0.45
     - 0.55
   * - GaSb
     - 5.72
     - 0.11
     - 0.77
     - +0.66
     - 1.48
     - 4.31
     - +2.83
     - 0.44
     - 0.56
   * - InP
     - 5.63
     - 0.46
     - 1.44
     - +0.98
     - 1.75
     - 4.99
     - +3.24
     - 0.45
     - 0.55
   * - InAs
     - 5.70
     - −0.42
     - 0.35
     - +0.77
     - 1.60
     - 4.61
     - +3.01
     - 0.43
     - 0.57
   * - InSb
     - 5.95
     - −0.61
     - 0.23
     - +0.84
     - 1.56
     - 4.30
     - +2.74
     - 0.42
     - 0.58

