Cluster size
============

Part of :doc:`/reference/index`.

Every run prints the size of the dot right after reading the geometry, whether or not the QP or
excitation model uses it, and writes it to ``cluster_size.json`` (and into ``qp_provenance.json``).
The NAMD precompute prints it for the first frame.

.. code-block:: text

   --- Cluster size ---
     Atoms              : 149 (Cd 68, Cl 26, Se 55)
     Core               : Cd68Se55
     SAXS-visible atoms : 149 (Cd, Cl, Se)
     Diameter (SAXS)    : 1.93 nm   (Debye I(q), sphere fit; reference size)
     Diameter (Guinier) : 1.93 nm   (R_g = 7.48 A)
     Diameter (volume)  : 1.87 nm   (core formula units x bulk volume)
     Diameter (hull)    : 1.84 nm   (core convex hull + 1.25 A)
     Principal extents  : 1.81 x 1.75 x 1.57 nm (anisotropy 1.15)
     Used by the model  : none (the size enters only through the orbitals)

**Reference size: the SAXS diameter.** Experimental sizing curves (Aubert, Hens et al. 2022) use the
diameter from small-angle X-ray scattering. SAXS sees the electron-density contrast with the solvent:
the inorganic part scatters, organic ligands (about the solvent's density) do not. QDEX computes what a
SAXS measurement of the cluster would give:

.. math::

   I(q) = \sum_{ij} Z_i Z_j\,\frac{\sin(q r_{ij})}{q r_{ij}},\qquad
   I(q) \approx I_0\,\Big[\frac{3(\sin qR - qR\cos qR)}{(qR)^3}\Big]^2 ,\qquad d_{\mathrm{SAXS}} = 2R,

summed over the inorganic atoms and fitted over the first lobe of the sphere form factor. The
Guinier value, :math:`d = 2\sqrt{5/3}\,R_g` with the electron-weighted radius of gyration, is printed
beside it; for compact dots the two agree.

* **Visible atoms.** By default all elements except H, C, N, O, P, B, Si and F; inorganic surface
  atoms such as halides count. ``system.inorganic_elements`` (``--inorganic-elements``) sets the list
  explicitly.
* **Other definitions.** The formula-unit diameter, :math:`(6 N_{\mathrm{fu}} V_{\mathrm{fu}}/\pi)^{1/3}`,
  counts only the core; the hull diameter is the volume-equivalent convex hull of the core atom centres
  plus 1.25 Å.
* **Which size the model uses.** The QP models that need a radius (``brus``, ``gw``, the sphere
  reaction field of the Resta and DIM models) use the hull radius; the block says so. ``bulk``,
  ``none`` and the sTDA use no radius.
* **Cost.** One pass over the pair distances: 0.1 s for 150 atoms, 4 s for 10⁴ atoms.

``benchmarks/compare_models.py --exp-ref`` evaluates the sizing curves at the SAXS diameter.
