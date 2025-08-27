#!/usr/bin/env bash
# Launch all windows listed in a windows.list file.
# Each line: /abs/or/rel/path/to/lambda_xx <gpu_id>
# Requires CUDA MPS (start once per node): nvidia-cuda-mps-control -d
set -euo pipefail
winlist="${1:?usage: launch_rbfe_windows.sh path/to/windows.list}"
if ! pgrep -f "nvidia-cuda-mps-control" >/dev/null 2>&1; then
  nvidia-cuda-mps-control -d || true
fi
# Concurrency: run all lines; adjust by head -n <N> to limit
while read -r wdir gid; do
  (
    cd "$wdir"
    export CUDA_VISIBLE_DEVICES="$gid"
    # Equil
    gmx mdrun -deffnm eq -ntmpi 1 -ntomp 1 -pin on -nb gpu -pme gpu -bonded gpu -update gpu
    # Production
    gmx mdrun -deffnm prod -ntmpi 1 -ntomp 1 -pin on -nb gpu -pme gpu -bonded gpu -update gpu
  ) &
done < "$winlist"
wait
echo "All windows in $winlist finished."