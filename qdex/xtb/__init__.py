"""
g-xTB interface for QDEX.

* ``calculator``: ASE calculator that runs the xtb binary with ``--gxtb``.
* ``md``: ASE molecular dynamics that writes one QDEX frame per production step.
* ``molden``: converts an xtb molden file into a per-atom basis file and ``MOs.mbse``.

g-xTB uses a charge-dependent (q-vSZP) basis, so the contraction coefficients change
from atom to atom and from frame to frame. Every frame therefore carries its own basis.
"""
