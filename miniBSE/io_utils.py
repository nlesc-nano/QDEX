import os
import re
import time
import math
import numpy as np
import collections
from scipy.sparse import issparse, csr_matrix
import libint_cpp


from miniBSE.constants import BOHR_PER_ANG

BOHR_PER_ANGSTROM = BOHR_PER_ANG


# ============================================================
# XYZ READER
# ============================================================

def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()

    nat = int(lines[0].strip())
    syms = []
    coords = []

    for line in lines[2:2 + nat]:
        p = line.split()
        syms.append(p[0])
        coords.append([float(x) for x in p[1:4]])

    return syms, np.asarray(coords)


# ============================================================
# BASIS PARSER (CP2K MOLOPT)
# ============================================================
import collections
import numpy as np

def parse_basis(fname, wanted, required_elements=None):
    basis = collections.defaultdict(list)
    required_elements = None if required_elements is None else set(required_elements)
    
    with open(fname) as f:
        # Filter out comments and completely empty lines upfront
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]
        
    it = iter(lines)
    
    for line in it:
        parts = line.split()
        if len(parts) < 2:
            continue
            
        elem = parts[0]
        bnames = parts[1:]
        
        # Prefer an exact CP2K basis name.  A unique prefix remains supported
        # for compatibility with element-specific suffixes such as ``-q13``.
        relevant_element = required_elements is None or elem in required_elements
        matching_names = (
            [wanted] if wanted in bnames else [b for b in bnames if b.startswith(wanted)]
        ) if relevant_element else []
        if len(matching_names) > 1:
            raise ValueError(
                f"Ambiguous basis prefix '{wanted}' for element {elem}: "
                + ", ".join(matching_names)
            )
        match_found = bool(matching_names)
        
        if not match_found:
            # Safely skip this block
            try:
                nset = int(next(it).split()[0])
                for _ in range(nset):
                    hdr = next(it).split()
                    nexp = int(hdr[3])
                    for _ in range(nexp):
                        next(it)
            except StopIteration:
                break
            continue
            
        # We found a matching element + basis name.
        if elem in basis:
            raise ValueError(
                f"Ambiguous basis selection for element {elem}: more than one block matches '{wanted}'. "
                "Use the complete CP2K basis name."
            )

        print(f"  [Basis] {elem}: selected {matching_names[0]}")

        # Extract the matched basis
        try:
            nset = int(next(it).split()[0])
            for _ in range(nset):
                hdr = next(it).split()
                lmin = int(hdr[1])
                nexp = int(hdr[3])
                counts = list(map(int, hdr[4:]))
                
                exps_list = []
                coef_rows = []
                for _ in range(nexp):
                    row = next(it).split()
                    exps_list.append(float(row[0]))
                    coef_rows.append([float(c) for c in row[1:]])
                    
                exps = np.array(exps_list)
                coef_cols = np.array(coef_rows).T
                
                idx = 0
                for j, n_shells in enumerate(counts):
                    l = lmin + j
                    for _ in range(n_shells):
                        coefs = coef_cols[idx].copy()
                        basis[elem].append((l, exps.copy(), coefs))
                        idx += 1
        except StopIteration:
            break

    return basis



# ============================================================
# BUILD SHELL DICTS FOR LIBINT
# ============================================================

def build_shell_dicts(syms, coords_ang, basis_dict):

    shells = []

    for atom_idx, (sym, xyz_ang) in enumerate(zip(syms, coords_ang)):

        if sym not in basis_dict:
            raise KeyError(f"No basis for element {sym}")

        xyz_bohr = np.asarray(xyz_ang) * BOHR_PER_ANGSTROM

        for l, exps, coefs in basis_dict[sym]:
            shells.append(dict(
                sym=sym,
                atom_idx=atom_idx,
                l=int(l),
                exps=np.asarray(exps, dtype=float),
                coefs=np.asarray(coefs, dtype=float),
                center=xyz_bohr,
                pure=True
            ))

    return shells


# ============================================================
# AO COUNT
# ============================================================

def count_ao_from_shells(shells):
    return sum(2 * int(sh["l"]) + 1 for sh in shells)


# ============================================================
# AO RANGES PER ATOM
# ============================================================

def build_atom_ao_ranges(shells):

    atom_ranges = {}
    ao_counter = 0

    for sh in shells:
        atom_idx = sh["atom_idx"]
        nbf = 2 * int(sh["l"]) + 1

        if atom_idx not in atom_ranges:
            atom_ranges[atom_idx] = [ao_counter, ao_counter + nbf]
        else:
            atom_ranges[atom_idx][1] += nbf

        ao_counter += nbf

    return [tuple(atom_ranges[i]) for i in sorted(atom_ranges.keys())]


# ============================================================
# MO UTILITIES
# ============================================================

_NUM_RE = re.compile(r"""
    [\+\-]?
    (?:
        \d+\.\d* |
        \.\d+ |
        \d+
    )
    (?:[EeDd][\+\-]?\d+)?
""", re.VERBOSE)


def _extract_numbers(s):
    toks = _NUM_RE.findall(s)
    return [float(t.replace("D", "E").replace("d", "E")) for t in toks]


import struct

def read_mos_mbse(path, n_ao_total, verbose=False):
    """
    Reads molecular orbitals from a binary .mbse (MBSEMO) file.

    Format specification (MBSEMO2):
      - 8 bytes magic: b'MBSEMO2\\x00' (or b'MBSEMO1\\x00')
      - 4 bytes uint32: version (e.g. 2)
      - 4 bytes uint32: endianness marker (0x01020304)
      - 8 bytes uint64: header_size (88)
      - 8 bytes uint64: n_ao
      - 8 bytes uint64: n_mo
      - 8 bytes uint64: n_occ
      - 8 bytes uint64: n_elec
      - 8 bytes uint64: flags
      - 8 bytes uint64: offset_C
      - 8 bytes uint64: offset_eps
      - 8 bytes uint64: offset_occ
      - Data at offset_C: n_ao * n_mo float64s in Fortran (column-major) order: C[:, i] for MO i
      - Data at offset_eps: n_mo float64s
      - Data at offset_occ: n_mo float64s
    """
    t0 = time.perf_counter()
    with open(path, "rb") as f:
        header = f.read(88)
        if len(header) < 88:
            raise ValueError(f"Corrupt or truncated MBSE MO file '{path}': header < 88 bytes.")

        magic = header[:8]
        if not magic.startswith(b"MBSEMO"):
            raise ValueError(f"Invalid MBSE MO magic header in '{path}': expected b'MBSEMO...', got {magic}")

        version, endian = struct.unpack("<II", header[8:16])
        header_size, n_ao, n_mo, n_occ, n_elec, flags, off_C, off_eps, off_occ = struct.unpack("<9Q", header[16:88])

        if n_ao != n_ao_total:
            raise ValueError(
                f"AO count mismatch in {path}: file contains {n_ao} AOs, "
                f"but current basis set requires {n_ao_total} AOs."
            )

        f.seek(off_C)
        raw_C = np.fromfile(f, dtype=np.float64, count=n_ao * n_mo)
        if raw_C.size != n_ao * n_mo:
            raise ValueError(f"Truncated C matrix in {path}: expected {n_ao * n_mo} elements, got {raw_C.size}.")
        C = raw_C.reshape((n_ao, n_mo), order="F")

        f.seek(off_eps)
        eps = np.fromfile(f, dtype=np.float64, count=n_mo)
        if eps.size != n_mo:
            raise ValueError(f"Truncated eps array in {path}: expected {n_mo} elements, got {eps.size}.")

        f.seek(off_occ)
        occ = np.fromfile(f, dtype=np.float64, count=n_mo)
        if occ.size != n_mo:
            raise ValueError(f"Truncated occ array in {path}: expected {n_mo} elements, got {occ.size}.")

    if verbose:
        dt = time.perf_counter() - t0
        size_gib = os.path.getsize(path) / (1024.0 ** 3)
        rate = size_gib / dt if dt > 0.0 else float("inf")
        print(
            f"[MOs:MBSE] Loaded {size_gib:.2f} GiB in {dt:.4f} s "
            f"({rate:.2f} GiB/s) | C shape {C.shape} (n_occ={n_occ}, n_mo={n_mo})"
        )

    return C, eps, occ


def write_mos_mbse(path, C, eps, occ, n_occ=None, n_elec=None, flags=1, version=2):
    """
    Writes molecular orbitals to a binary .mbse (MBSEMO2) file.
    """
    C = np.asfortranarray(C, dtype=np.float64)
    eps = np.asarray(eps, dtype=np.float64)
    occ = np.asarray(occ, dtype=np.float64)
    n_ao, n_mo = C.shape

    if eps.shape != (n_mo,):
        raise ValueError(f"eps shape {eps.shape} does not match n_mo {n_mo}")
    if occ.shape != (n_mo,):
        raise ValueError(f"occ shape {occ.shape} does not match n_mo {n_mo}")

    if n_occ is None:
        n_occ = int(np.sum(occ > 0.5))
    if n_elec is None:
        n_elec = int(np.round(np.sum(occ)))

    header_size = 88
    off_C = header_size
    off_eps = off_C + n_ao * n_mo * 8
    off_occ = off_eps + n_mo * 8

    magic = f"MBSEMO{version}\x00".encode("ascii")[:8]
    endian = 0x01020304

    with open(path, "wb") as f:
        f.write(magic)
        f.write(struct.pack("<II", version, endian))
        f.write(struct.pack("<9Q", header_size, n_ao, n_mo, n_occ, n_elec, flags, off_C, off_eps, off_occ))
        f.write(C.tobytes(order="F"))
        f.write(eps.tobytes())
        f.write(occ.tobytes())


def read_mos_auto(path, n_ao_total, verbose=False, cache=False):

    ext = os.path.splitext(path)[-1].lower()

    if ext == ".npz":
        d = np.load(path, allow_pickle=False)
        if "C" in d.files:
            C = d["C"]
        else:
            C = csr_matrix((d["data"], d["indices"], d["indptr"]),
                           shape=d["shape"])
        eps, occ = d["eps"], d["occ"]

        if verbose:
            print(f"[MOs] Loaded NPZ: {C.shape}")

        return C, eps, occ

    if ext == ".mbse":
        return read_mos_mbse(path, n_ao_total, verbose=verbose)

    # Check magic header for extension-agnostic detection
    if os.path.exists(path) and os.path.getsize(path) >= 88:
        try:
            with open(path, "rb") as f:
                head = f.read(6)
            if head == b"MBSEMO":
                return read_mos_mbse(path, n_ao_total, verbose=verbose)
        except Exception:
            pass

    cache_path = f"{path}.minibse.npz"
    source_stat = os.stat(path)
    if cache and os.path.exists(cache_path):
        with np.load(cache_path, allow_pickle=False) as d:
            valid = (
                int(d["source_size"]) == source_stat.st_size
                and int(d["source_mtime_ns"]) == source_stat.st_mtime_ns
                and int(d["n_ao_total"]) == n_ao_total
            )
            if valid:
                t0 = time.perf_counter()
                C, eps, occ = d["C"], d["eps"], d["occ"]
                if verbose:
                    print(
                        f"[MOs:cache] Loaded {cache_path} in {time.perf_counter() - t0:.4f} s "
                        f"| C shape {C.shape}"
                    )
                return C, eps, occ
        if verbose:
            print(f"[MOs:cache] Ignoring stale cache {cache_path}")

    C, eps, occ = read_mos_txt_cc(path, n_ao_total, verbose=verbose)
    if cache:
        t0 = time.perf_counter()
        np.savez(
            cache_path, C=C, eps=eps, occ=occ,
            source_size=np.int64(source_stat.st_size),
            source_mtime_ns=np.int64(source_stat.st_mtime_ns),
            n_ao_total=np.int64(n_ao_total),
        )
        if verbose:
            print(
                f"[MOs:cache] Wrote reusable binary cache {cache_path} "
                f"in {time.perf_counter() - t0:.2f} s"
            )
    return C, eps, occ

def read_mos_uks(path_alpha, path_beta, n_ao, verbose=False, cache=False):
    """
    Reads alpha and beta MO files from a CP2K UKS calculation.
    Returns:
        C_alpha, eps_alpha, occ_alpha  — alpha spin channel
        C_beta,  eps_beta,  occ_beta   — beta spin channel
    """
    C_alpha, eps_alpha, occ_alpha = read_mos_auto(path_alpha, n_ao, verbose=verbose, cache=cache)
    C_beta,  eps_beta,  occ_beta  = read_mos_auto(path_beta,  n_ao, verbose=verbose, cache=cache)
    return C_alpha, eps_alpha, occ_alpha, C_beta, eps_beta, occ_beta


def read_mos_txt_cc(path, n_ao_total, verbose=False):
    t0 = time.perf_counter()
    
    # 1 line of Python. C++ handles everything else.
    C, eps, occ = libint_cpp.parse_cp2k_mos(path, n_ao_total)
    
    if verbose:
        dt = time.perf_counter() - t0
        size_gib = os.path.getsize(path) / (1024.0 ** 3)
        rate = size_gib / dt if dt > 0.0 else float("inf")
        print(
            f"[MOs:C++] Parsed {size_gib:.2f} GiB in {dt:.4f} s "
            f"({rate:.2f} GiB/s) | C shape {C.shape}"
        )
        
    return C, eps, occ

import collections

def parse_gth_soc_potentials(path, elements_to_parse):
    """Parses GTH potentials for SOC parameters."""
    ecp_dict = collections.defaultdict(lambda: {'so': []})
    needed_elements = set(elements_to_parse.keys())

    def _collect_coeffs(line_iter, n, init):
        coeffs = list(init)
        while len(coeffs) < n:
            line = next(line_iter).strip()
            if line and not line.startswith('#'):
                coeffs.extend([float(x) for x in line.split()])
        return coeffs

    try:
        with open(path, "r") as f:
            line_iter = iter(f.readlines())
    except FileNotFoundError:
        print(f"Warning: Potential file not found at {path}. SOC will be zero.")
        return ecp_dict

    for line in line_iter:
        if not needed_elements: break
        parts = line.strip().split()
        if not parts or parts[0] not in needed_elements: continue
        
        sym, q = parts[0], elements_to_parse.get(parts[0])
        if q is None or not any(f"q{q}" in p for p in parts): continue
        
        try:
            next(line_iter); next(line_iter) # Skip header
            n_soc_sets = int(next(line_iter).strip().split()[0])
            for l in range(n_soc_sets):
                proj_line = next(line_iter)
                while not proj_line.strip() or proj_line.strip().startswith('#'):
                    proj_line = next(line_iter)
                proj_parts = proj_line.split()
                r, nprj = float(proj_parts[0]), int(proj_parts[1])
                n_coeffs = nprj * (nprj + 1) // 2
                h = _collect_coeffs(line_iter, n_coeffs, proj_parts[2:])
                k = _collect_coeffs(line_iter, n_coeffs, []) if l > 0 else []
                ecp_dict[sym]['so'].append({'l': l, 'r': r, 'nprj': nprj, 'h_coeffs': h, 'k_coeffs': k})
            needed_elements.remove(sym)
        except Exception as e:
            continue
            
    return ecp_dict

def get_vxc_ao_matrix(txt_path, n_ao):
    """
    Retrieves the Vxc AO matrix. 
    Uses a raw .bin cache to instantly load massive arrays on subsequent runs.
    """
    import os
    import numpy as np
    import libint_cpp

    bin_path = txt_path + ".raw.bin"

    # 1. Try instant binary load via C++
    if os.path.exists(bin_path):
        print(f"  [Vxc] Found raw binary cache. Loading instantly via C++...")
        return libint_cpp.load_raw_binary(bin_path, n_ao)

    # 2. Fallback to C++ text parser
    print(f"  [Vxc] Cache not found. Parsing text block-matrix with C++ (First time only)...")
    V_ao = libint_cpp.parse_cp2k_block_matrix(txt_path, n_ao)

    # 3. Save as raw binary for future runs
    print(f"  [Vxc] Caching raw binary matrix for high-speed future access...")
    with open(bin_path, "wb") as f:
        f.write(V_ao.tobytes()) # Zero-overhead raw dump
        
    return V_ao
