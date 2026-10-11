"""Birch-Murnaghan fit of the HLE17 energy-volume scan (scan/s*/cp2k_job.out); prints a0 (A) on the last line."""
import glob, re, numpy as np
from scipy.optimize import curve_fit

HA_EV, A3_GPA = 27.211386245988, 160.21766208
def bm(V, E0, V0, B0, B1):
    x = (V0 / V) ** (2 / 3)
    return E0 + 9 * V0 * B0 / 16 * ((x - 1) ** 3 * B1 + (x - 1) ** 2 * (6 - 4 * x))
rows = []
for d in sorted(glob.glob("scan/s[0-9]*") + glob.glob("scan/exp")):
    out = open(f"{d}/cp2k_job.out").read()
    E = float(re.findall(r"ENERGY\| Total FORCE_EVAL.*?(-?\d+\.\d+)\s*$", out, re.M)[-1]) * HA_EV
    inp = open(f"{d}/cp2k_job.in").read()
    A = [list(map(float, re.search(rf"^\s+{c}\s+(\S+)\s+(\S+)\s+(\S+)", inp, re.M).groups())) for c in "ABC"]
    V = abs(np.linalg.det(np.array(A)))
    rows.append((V, E, (4 * V) ** (1 / 3)))
    print(f"{d}: a {(4 * V) ** (1 / 3):.4f} A, E {E:.6f} eV")
V, E, _ = map(np.array, zip(*rows))
p, _ = curve_fit(bm, V, E, p0=(E.min(), V[np.argmin(E)], 0.3, 4.5))
rms = np.sqrt(np.mean((bm(V, *p) - E) ** 2)) * 1000
print(f"BM fit: B0 {p[2] * A3_GPA:.1f} GPa, B1 {p[3]:.2f}, rms {rms:.3f} meV")
print(f"{(4 * p[1]) ** (1 / 3):.4f}")
