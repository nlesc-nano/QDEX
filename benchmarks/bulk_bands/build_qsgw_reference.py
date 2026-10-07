"""Bulk self-energy correction Delta_Sigma from the literature QSGW+SOC gaps and our PBE+SOC gaps at a_exp.

    python build_qsgw_reference.py [PBE_LATTICE_SOC_JSON]

Inputs: QSGW_SOC_LITERATURE and MATERIAL_DB_BULK_PBE (qdex/hardness.py), exp_lattice_soc_gaps.json
(signed spin-free and PBE+SOC gaps at a_exp, summarize_soc_gaps.py on the make_exp_lattice_inputs.py
runs) and, optionally, the same summary at the PBE lattice (qdex/data/bulk_bands).

    E_QSGW+SOC(a_exp) = E_QSGW+SOC(a_lit) + D * 3 ln(a_exp / a_lit)      (D: dEg/dlnV of the source, else PBE)
    Delta_Sigma       = E_QSGW+SOC(a_exp) - E_PBE+SOC(a_exp)               (same lattice, same SOC treatment)
    MATERIAL_DB[8]    = E_PBE,SF(a_exp) + Delta_Sigma                       (spin-free QSGW-equivalent gap)

so that Delta_bulk = index 8 - index 7 = Delta_Sigma: one correction for the spin-free and the SOC dots,
and the spin-orbit lowering of the gap is the PBE one, as in the dots. Prints the MATERIAL_DB_BULK_PBE
entries and the new index-8 values.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from make_lattice_scan import NAMES  # noqa: E402
from qdex.hardness import MATERIAL_DB, MATERIAL_DB_BULK_PBE, QSGW_SOC_LITERATURE  # noqa: E402

if __name__ == "__main__":
    exp = json.load(open(os.path.join(HERE, "exp_lattice_soc_gaps.json")))
    pbe = json.load(open(sys.argv[1])) if len(sys.argv) > 1 else {}
    print(f"{'mat':8s} {'QSGW+SO':>8s} {'a_lit':>6s} {'a_exp':>6s} {'->a_exp':>8s} | {'PBE_SF':>7s} {'PBE_SOC':>8s} "
          f"{'dSOC':>6s} | {'D_Sigma':>8s} {'idx8 new':>8s} {'idx8 old':>8s} {'D_old':>6s}")
    lines = []
    for key, name in {**NAMES, "CDSE_WZ": "CdSe_wz", "CDS_WZ": "CdS_wz"}.items():
        b = MATERIAL_DB_BULK_PBE[key]
        e = exp[name]
        a_exp, a_pbe = b["a_exp"], b["a_pbe"]
        d_pbe = (b["gap_pbe_lattice"] - b["gap_exp_lattice"]) / (3.0 * math.log(a_pbe / a_exp))
        if key in QSGW_SOC_LITERATURE:
            g_lit, a_lit, _, src, d_src = QSGW_SOC_LITERATURE[key]
            g_soc = g_lit + (d_src if d_src is not None else d_pbe) * 3.0 * math.log(a_exp / a_lit)
            d_sigma = g_soc - e["gap_soc"]
        else:                                    # no literature value: keep the stored index 8
            g_lit, a_lit, src = float("nan"), float("nan"), "stored"
            d_sigma = MATERIAL_DB[key][8] - e["gap_sf"]
            g_soc = e["gap_soc"] + d_sigma
        idx8 = e["gap_sf"] + d_sigma
        p = pbe.get(name, {})
        print(f"{key:8s} {g_lit:8.3f} {a_lit:6.3f} {a_exp:6.3f} {g_soc:8.3f} | {e['gap_sf']:7.3f} {e['gap_soc']:8.3f} "
              f"{e['delta_soc']:6.3f} | {d_sigma:8.3f} {idx8:8.3f} {MATERIAL_DB[key][8]:8.2f} "
              f"{MATERIAL_DB[key][8] - MATERIAL_DB[key][7]:6.3f}  {src}")
        lines.append(f'    "{key}": dict(a_exp={a_exp:.3f}, gap_exp_lattice={b["gap_exp_lattice"]:.3f}, '
                     f'gap_soc_exp_lattice={e["gap_soc"]:.3f}, a_pbe={a_pbe:.3f}, gap_pbe_lattice={b["gap_pbe_lattice"]:.3f}, '
                     f'gap_soc_pbe_lattice={p["gap_soc"]:.3f}, ' if p else
                     f'    "{key}": dict(a_exp={a_exp:.3f}, gap_exp_lattice={b["gap_exp_lattice"]:.3f}, '
                     f'gap_soc_exp_lattice={e["gap_soc"]:.3f}, a_pbe={a_pbe:.3f}, gap_pbe_lattice={b["gap_pbe_lattice"]:.3f}, ')
        lines[-1] += f'qsgw_soc_exp_lattice={g_soc:.3f}, delta_sigma={d_sigma:.3f}, B0={b["B0"]:.1f}),'
    print("\n".join(lines))
