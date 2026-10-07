"""Signed spin-free and PBE+SOC gaps from the qdex.bulk_soc metadata (<name>.json).

    python summarize_soc_gaps.py JSON_DIR OUT_JSON

The gap is the band-edge transition of the material, signed so that an inverted band order is
negative: zinc blende with a direct gap, E(Gamma6) - E(Gamma8) (E(Gamma1) - E(Gamma15) spin-free) from
the s-p band order; rock salt (L) and cubic perovskites (R), the edge band order; indirect zinc
blende (AlP, AlAs, AlSb, GaP), the fundamental gap. delta_soc = gap_sf - gap_soc.
"""
import glob
import json
import os
import sys


def signed_gaps(j):
    sf, soc = j["sf"]["gap"], j["soc"]["gap"]
    direct_gamma = all(abs(x) < 1e-6 for x in j["sf"]["cbm_k"])
    if j.get("gamma_band_order") and (direct_gamma or j["gamma_band_order"].get("inverted")):
        sf, soc = j["gamma_band_order"]["sf"], j["gamma_band_order"]["soc"]
    elif j.get("edge_band_order") and j["edge_band_order"].get("reliable", True):
        sf, soc = j["edge_band_order"]["sf"], j["edge_band_order"]["soc"]
    trip = j.get("gamma_vb_triplet") or j.get("cbm_triplet") or {}
    return dict(gap_sf=round(sf, 4), gap_soc=round(soc, 4), delta_soc=round(sf - soc, 4),
                delta_so=round(trip["delta_so"], 4) if "delta_so" in trip else None,
                inverted_soc=bool(soc < 0), a=j.get("lattice", {}).get("a"))


if __name__ == "__main__":
    src, out = sys.argv[1], sys.argv[2]
    table = {}
    for f in sorted(glob.glob(os.path.join(src, "*.json"))):
        j = json.load(open(f))
        if "sf" in j and "soc" in j:
            table[j.get("name", os.path.basename(f)[:-5])] = signed_gaps(j)
    json.dump(table, open(out, "w"), indent=1)
    for k, v in table.items():
        print(f"{k:14s} SF {v['gap_sf']:7.3f}  SOC {v['gap_soc']:7.3f}  delta_soc {v['delta_soc']:6.3f}")
