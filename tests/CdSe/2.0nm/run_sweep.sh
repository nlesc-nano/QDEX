#!/bin/bash
# CdSe 2.0 nm model sweep (Delta-W QP models use the shared-W kernel "qp").
# For the organised comparison with checks use benchmarks/compare_models.py.
# Originally reported in
# audit/AUDIT_2026-09-25_QP_EXCITED_STATES.md (spin-free, 4 threads).
#   cd tests/CdSe && ./run_sweep.sh [outdir]
set -u
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here/sweep}"
mkdir -p "$out"
[ -f "$here/MOs_cleaned_20ang.txt" ] || gunzip -k "$here/MOs_cleaned_20ang.txt.gz"

run() {
  name=$1; shift
  d="$out/$name"; mkdir -p "$d"
  for f in MOs_cleaned_20ang.txt geom.xyz BASIS_MOLOPT_UZH GTH_SOC_POTENTIALS.txt; do ln -sf "$here/$f" "$d/$f"; done
  sed 's/soc_flag: true/soc_flag: false/; s/nthreads: 12/nthreads: 4/' "$here/config.yaml" > "$d/config.yaml"
  [ -n "${DAV:-}" ] && sed -i 's/full_diag: true/full_diag: false/' "$d/config.yaml"
  (cd "$d" && qdex --config config.yaml "$@" > run.out 2>&1; echo "EXIT $?" >> run.out)
  qp=$(grep "Final QP gap" "$d/run.out" | tail -1 | awk '{print $(NF-1)}')
  s1=$(awk -F, 'NR==2{print $3}' "$d/exciton_results.csv" 2>/dev/null)
  printf "%-22s QP=%-7s S1=%s\n" "$name" "$qp" "$s1"
}

run gw_eps1
run gw_eps2.4 --eps-out 2.4
run brus --qp_gap brus
run pbe --qp_gap pbe
run sgw-resta --qp_gap sgw-resta --kernel qp
run sgw-resta_eps2.4 --qp_gap sgw-resta --eps-out 2.4 --kernel qp
run sgw-resta-pure --qp_gap sgw-resta-pure --kernel qp
run sgw-dim --qp_gap sgw-dim --kernel qp
run evgw-resta --qp_gap evgw-resta --kernel qp
run qsgw-resta --qp_gap qsgw-resta --kernel qp
run qsgw-dim --qp_gap qsgw-dim --kernel qp
run sgw_sbse --qp_gap sgw --kernel qp
run diag_bse --excitation-mode diagonal_bse
run indep_qp --excitation-mode independent_qp
run as50 --nhomos 50 --nlumos 50
run k_dim --qp_gap gw --kernel dim
run k_sbse --qp_gap gw --kernel sbse
run k_mnok --qp_gap gw --kernel bse
run k_xsresta --qp_gap gw --kernel xs-resta
run lowdin --charge_type lowdin
run triplet --triplet
DAV=1 run as100_davidson --nhomos 100 --nlumos 100 --nroots 10
