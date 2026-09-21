import numpy as np
import time
import libint_cpp

def compute_dipole_ao(shells, origin=None, nthreads=1):
    """
    Computes AO dipole matrices (mu_x, mu_y, mu_z) using libint_cpp.
    """
    if origin is None:
        origin = [0.0, 0.0, 0.0]
    
    # libint_cpp.dipole returns a (3, nbf, nbf) numpy array
    mu_ao = libint_cpp.dipole(shells, origin, nthreads)
    return mu_ao[0], mu_ao[1], mu_ao[2]


def compute_cross_overlap_ao(shells1, shells2, nthreads=1):
    """
    Computes cross-overlap between two different geometries in the same basis set.
    Returns (nbf1, nbf2) numpy array.
    """
    return libint_cpp.cross_overlap_geometries(shells1, shells2, nthreads)


