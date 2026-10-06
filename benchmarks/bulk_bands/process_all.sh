#!/bin/bash
# Build <name>.bs.gz, <name>_soc.bs.gz and <name>.json in qdex/data/bulk_bands from the CP2K runs.
#   process_all.sh RUNS_DIR CIF_ROOT BASIS_MOLOPT_UZH GTH_SOC_POTENTIALS [names...]
# RUNS_DIR/<name>/ holds the CP2K output of benchmarks/bulk_bands/inputs/<name> (*.csr, *.bs,
# *_trexio.h5, cp2k_job.out).
set -u
RUNS=$1; CIFS=$2; BASIS=$3; GTH=$4; shift 4
REPO=$(cd "$(dirname "$0")/../.." && pwd)
names=${@:-$(ls "$RUNS")}
for n in $names; do
  cif=$(ls "$CIFS"/*/bulk_cifs/$n.cif | head -1)
  (cd "$REPO" && python -m qdex.bulk_soc "$RUNS/$n" --basis "$BASIS" --gth "$GTH" --name $n --cif "$cif" \
      -d "$REPO/qdex/data/bulk_bands" --gzip 2>&1 | grep "Bulk SOC\]\|Error\|Traceback" | sed "s/^/$n: /")
done
