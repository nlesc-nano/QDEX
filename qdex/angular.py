"""Angular coverage Omega of a state of a quantum dot, and the S / P / D / Facet / Localized classes.

With normalized atom populations p_A and directions u_A = (R_A - c)/|R_A - c| from the centre c of the inorganic
core, the angular density rho(u) = sum_A p_A delta(u - u_A) has multipoles a_lm = sum_A p_A Y*_lm(u_A) and angular
power C_l = sum_m |a_lm|^2 (rotation invariant), and

    Omega_L = C_0 / sum_{even l <= L} C_l,

the fraction of directions the even (inversion-symmetric) part of the density covers at angular resolution L, an
angular participation ratio. Only even l enter: the density of an envelope state of angular momentum l has even
multipoles only (up to 2l), while the odd ones come from parity breaking (the dot's internal field pushing a state
to one side, l = 1; the tetrahedral shape, l = 3) and would blur the S / P / D symmetry. Their weight,
sum_{odd l} C_l / C_0, is large for a patch on one facet (3-6 in CdSe dots) and small for envelope states (<1).

Omega is 1 for an isotropic density (1S envelope), 5/9 for |Y_1m|^2 (a 1P envelope, no density at the centre), and
1 / sum_{even l <= L} (2l+1) for a single atom; an atom closer than half a bond to the centre has no direction and
counts as isotropic. L grows with the dot, L = max(4, round(pi R / (2 d))), so that the resolution stays about two
bonds d on the surface.

QDEX reports Omega relative to the ideal 1S envelope evaluated on the dot's own atoms, min(Omega_L / Omega_1S, 1), so
that 1 means "as isotropic as a 1S envelope" on that lattice (the raw 1S value is 0.8-0.98: the lattice is discrete
and the dot faceted).

The class boundaries come from reference densities evaluated on the dot's own atoms with the same L: hard-wall
envelopes 1S, 1P, 1D (|j_l(x_l1 r / R_env) Y_l0|^2, R_env = R + d/2), a cap of one facet (the outer layer within
the half-angle of 1/8 of the sphere, averaged over directions) and a single atom:

    S-like     Omega >= (S + P)/2
    P-like     (P + D)/2 <= Omega < (S + P)/2
    D-like     (D + facet)/2 <= Omega < (P + D)/2
    Facet      sqrt(facet * atom) <= Omega < (D + facet)/2
    Localized  Omega < sqrt(facet * atom)
"""
import numpy as np

CLASS_NAMES = ("S-like", "P-like", "D-like", "Facet", "Localized")
BAND_CLASSES = ("S-like", "P-like", "D-like")
PASSIVANTS = {"H", "C", "N", "O", "F", "Cl", "Br", "I", "P", "S"}


def _sph_harm(l, m, theta, phi):
    try:
        from scipy.special import sph_harm_y
        return sph_harm_y(l, m, theta, phi)
    except ImportError:                       # scipy < 1.15
        from scipy.special import sph_harm
        return sph_harm(m, l, phi, theta)


def harmonics(dirs, L):
    """Y_lm (row per (l, m), l = 0..L) at the unit vectors dirs (n, 3); also the l of every row."""
    dirs = np.asarray(dirs, dtype=float)
    theta = np.arccos(np.clip(dirs[:, 2], -1.0, 1.0))
    phi = np.arctan2(dirs[:, 1], dirs[:, 0])
    rows, ls = [], []
    for l in range(L + 1):
        for m in range(-l, l + 1):
            rows.append(_sph_harm(l, m, theta, phi))
            ls.append(l)
    return np.array(rows), np.array(ls)


def angular_power(weights, Y, ls):
    """C_l (l = 0..L, rows) of the angular densities whose atom weights are the columns of ``weights``."""
    w = np.atleast_2d(np.asarray(weights, dtype=float).T).T
    a = np.conj(Y) @ w                                     # (n_lm, n_state)
    power = np.abs(a) ** 2
    return np.array([power[ls == l].sum(axis=0) for l in range(int(ls.max()) + 1)])


def omega_of(weights, Y, ls, parity="even"):
    """Omega_L of the angular densities whose atom weights are the columns of ``weights`` (n_atom, n_state):
    C_0 over the sum of C_l of the even l (parity="even", default) or of all l (parity="all")."""
    C = angular_power(weights, Y, ls)
    keep = (np.arange(len(C)) % 2 == 0) if parity == "even" else np.ones(len(C), dtype=bool)
    return C[0] / np.maximum(C[keep].sum(axis=0), 1e-300)


def odd_power(weights, Y, ls):
    """Sum of C_l over odd l relative to C_0: the one-sidedness of the density (polarization, a patch)."""
    C = angular_power(weights, Y, ls)
    return C[1::2].sum(axis=0) / np.maximum(C[0], 1e-300)


def core_mask(syms, core_elements=None):
    syms = [str(s).strip().capitalize() for s in syms]
    if core_elements:
        core = {str(e).strip().capitalize() for e in core_elements}
    else:
        core = {s for s in syms if s not in PASSIVANTS} or set(syms)
    return np.array([s in core for s in syms], dtype=bool)


class AngularCoverage:
    """Geometry of one dot: centre, radius, L, harmonics and the class boundaries of Omega."""

    def __init__(self, coords_ang, syms, core_elements=None, L=None, n_cap_dirs=24):
        X = np.asarray(coords_ang, dtype=float)
        self.core = core_mask(syms, core_elements)
        self.centre = X[self.core].mean(axis=0)
        v = X - self.centre
        r = np.linalg.norm(v, axis=1)
        self.R = float(r[self.core].max())
        Xc = X[self.core]
        dd = np.linalg.norm(Xc[:, None, :] - Xc[None, :, :], axis=2) if len(Xc) <= 4000 else None
        if dd is not None and len(Xc) > 1:
            np.fill_diagonal(dd, np.inf)
            self.bond = float(np.median(dd.min(axis=1)))
        else:
            self.bond = 2.6
        self.L = int(L) if L else max(4, int(round(np.pi * self.R / (2.0 * self.bond))))
        self.dirs = v / np.maximum(r[:, None], 1e-9)
        self.Y, self.ls = harmonics(self.dirs, self.L)
        # an atom at the centre has no direction: its weight counts as isotropic (l = 0 only)
        self.Y[np.ix_(self.ls > 0, r < 0.5 * self.bond)] = 0.0
        # Omega is reported relative to the ideal 1S envelope on this dot's atoms, so that the scale 0..1 is used
        # in full (a discrete, faceted dot has raw Omega_1S of 0.8-0.98) and 1 means "as isotropic as a 1S envelope"
        self.norm = 1.0
        _, refs_raw = self._references(r, n_cap_dirs)
        self.norm = float(refs_raw["S"])
        self.thresholds, self.references = self._references(r, n_cap_dirs)

    def omega_raw(self, atom_weights, parity="even"):
        """Omega_L of each state as defined above; atom_weights (n_atom, n_state), negative parts clipped."""
        w = np.clip(np.asarray(atom_weights, dtype=float), 0.0, None)
        w = w / np.maximum(w.sum(axis=0, keepdims=True), 1e-300)
        return omega_of(w, self.Y, self.ls, parity=parity)

    def omega(self, atom_weights, parity="even"):
        """Omega relative to the ideal 1S envelope of this dot, min(Omega_L / Omega_1S, 1)."""
        return np.minimum(self.omega_raw(atom_weights, parity=parity) / self.norm, 1.0)

    def odd(self, atom_weights):
        w = np.clip(np.asarray(atom_weights, dtype=float), 0.0, None)
        w = w / np.maximum(w.sum(axis=0, keepdims=True), 1e-300)
        return odd_power(w, self.Y, self.ls)

    def _references(self, r, n_cap_dirs):
        from scipy.special import spherical_jn
        from scipy.optimize import brentq
        core = self.core
        R_env = self.R + 0.5 * self.bond
        refs = {}
        for name, l in (("S", 0), ("P", 1), ("D", 2)):
            grid = np.arange(0.5, 12.0, 0.01)
            f = spherical_jn(l, grid)
            i = int(np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0])
            x = brentq(lambda t: spherical_jn(l, t), grid[i], grid[i + 1])
            amp = spherical_jn(l, x * np.minimum(r, R_env) / R_env)
            if l:
                amp = amp * np.real(_sph_harm(l, 0, np.arccos(np.clip(self.dirs[:, 2], -1, 1)),
                                              np.arctan2(self.dirs[:, 1], self.dirs[:, 0])))
            w = np.where(core, amp ** 2, 0.0)
            refs[name] = float(self.omega(w[:, None])[0])
        # one facet: outer layer within a cone covering 1/8 of the sphere (1 - cos(theta) = 1/4), averaged over axes
        outer = core & (r > self.R - self.bond)
        k = np.arange(n_cap_dirs) + 0.5
        pol = np.arccos(1 - 2 * k / n_cap_dirs)
        az = np.pi * (1 + 5 ** 0.5) * k
        axes = np.stack([np.cos(az) * np.sin(pol), np.sin(az) * np.sin(pol), np.cos(pol)], axis=1)
        caps = (self.dirs @ axes.T >= 0.75) & outer[:, None]
        caps = caps[:, caps.sum(axis=0) > 0].astype(float)
        refs["facet"] = float(np.mean(self.omega(caps))) if caps.size else 0.1
        refs["atom"] = 1.0 / sum(2 * l + 1 for l in range(0, self.L + 1, 2)) / self.norm
        thr = {
            "S": 0.5 * (refs["S"] + refs["P"]),
            "P": 0.5 * (refs["P"] + refs["D"]),
            "D": 0.5 * (refs["D"] + refs["facet"]),
            "Facet": float(np.sqrt(refs["facet"] * refs["atom"])),
        }
        return thr, refs

    def classify(self, omega):
        """Class name of each Omega (CLASS_NAMES)."""
        t = self.thresholds
        out = []
        for o in np.atleast_1d(omega):
            if o >= t["S"]:
                out.append("S-like")
            elif o >= t["P"]:
                out.append("P-like")
            elif o >= t["D"]:
                out.append("D-like")
            elif o >= t["Facet"]:
                out.append("Facet")
            else:
                out.append("Localized")
        return out

    def summary(self):
        t, r = self.thresholds, self.references
        note = ("; the dot is too small for S, P and D envelopes to differ in Omega"
                if r["S"] - r["D"] < 0.15 else "")
        return (f"L = {self.L} (R = {self.R:.2f} A, bond {self.bond:.2f} A), relative to the ideal 1S value {self.norm:.2f}; "
                f"references S {r['S']:.2f}, P {r['P']:.2f}, "
                f"D {r['D']:.2f}, facet {r['facet']:.3f}, atom {r['atom']:.4f}; boundaries S >= {t['S']:.2f}, "
                f"P >= {t['P']:.2f}, D >= {t['D']:.2f}, facet >= {t['Facet']:.3f}{note}")


def atom_populations(pops, shells, n_atoms):
    """Atom populations (n_atoms, n_states) summed from Mulliken AO populations (n_ao, n_states)."""
    ao_atom = np.concatenate([[int(sh.get("atom_idx", 0))] * (2 * int(sh["l"]) + 1) for sh in shells])
    P = np.zeros((n_atoms, np.asarray(pops).shape[1]))
    np.add.at(P, ao_atom, np.asarray(pops, dtype=float))
    return P


_CACHE = {}


def coverage_for(coords_ang, syms, core_elements=None):
    """AngularCoverage of a geometry, cached (the same dot is classified several times per run)."""
    X = np.asarray(coords_ang, dtype=float)
    key = (X.shape, float(X.sum()), float(np.abs(X).sum()), tuple(core_elements or ()))
    if key not in _CACHE:
        _CACHE.clear()
        _CACHE[key] = AngularCoverage(X, syms, core_elements)
    return _CACHE[key]


def states_omega(pops, shells, coords_ang, syms, core_elements=None, chunk=2000):
    """Omega, class and IPR (sum_mu P_mu^2 of the normalized AO populations) of the states (columns of pops)."""
    P = np.asarray(pops)
    cov = coverage_for(coords_ang, syms, core_elements)
    n = P.shape[1]
    om, ipr = np.empty(n), np.empty(n)
    for i0 in range(0, n, chunk):                    # chunks: atom populations of every MO of a large dot
        Pc = np.asarray(P[:, i0:i0 + chunk], dtype=float)
        om[i0:i0 + chunk] = cov.omega(atom_populations(Pc, shells, len(syms)))
        tot = Pc.sum(axis=0)
        tot = np.where(np.abs(tot) < 1e-12, 1.0, tot)
        ipr[i0:i0 + chunk] = np.sum((Pc / tot) ** 2, axis=0)
    return om, cov.classify(om), ipr, cov
