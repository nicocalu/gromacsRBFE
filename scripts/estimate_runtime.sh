#!/usr/bin/env bash
# Micro-benchmark a single window (200 ps) and project walltime per mutation.
# Usage: bash scripts/estimate_runtime.sh mutations/A33Y complex lambda_00 6
set -euo pipefail
mutdir="${1:?mutations/<MUT>}"
leg="${2:?complex|solvent}"
lam="${3:?lambda_xx}"
win_per_gpu="${4:-6}"
wdir="$mutdir/$leg/$lam"
test -f "$wdir/prod.tpr" || { echo "Missing $wdir/prod.tpr (run grompp_all.sh first)"; exit 1; }

tmpdir="$wdir/.bench"
mkdir -p "$tmpdir"
cp "$wdir/prod.tpr" "$tmpdir/bench.tpr"
pushd "$tmpdir" >/dev/null
# 200 ps = 100000 steps @2 fs
gmx mdrun -s bench.tpr -nsteps 100000 -noconfout -g bench.log -ntmpi 1 -ntomp 1 -pin on -nb gpu -pme gpu -bonded gpu -update gpu >/dev/null 2>&1 || true
popd >/dev/null

ns_per_day=$(grep -E "Performance:" "$tmpdir/bench.log" | tail -n1 | awk '{print $(NF-2)}' || echo "0")
if [[ "$ns_per_day" == "0" || -z "$ns_per_day" ]]; then
  echo "Could not read performance from bench.log"
  exit 1
fi

# Assume roughly linear sharing with MPS across windows/GPU
per_window_rate=$(python3 - <<PY
r=float("$ns_per_day")/float("$win_per_gpu")
print(f"{r:.2f}")
PY
)

# Detect profile by reading nsteps from mdps
eq_steps=$(grep -m1 "^nsteps" "$mutdir/$leg/$lam/eq.mdp" | awk '{print $3}')
prod_steps=$(grep -m1 "^nsteps" "$mutdir/$leg/$lam/prod.mdp" | awk '{print $3}')
eq_ns=$(python3 - <<PY
print(f"{int("$eq_steps")*0.002/1000.0:.3f}")
PY
)
prod_ns=$(python3 - <<PY
print(f"{int("$prod_steps")*0.002/1000.0:.3f}")
PY
)

time_per_window_hr=$(python3 - <<PY
rate=float("$per_window_rate")
eq=float("$eq_ns"); prod=float("$prod_ns")
tot=eq+prod
hours = 24.0*tot/rate
print(f"{hours:.2f}")
PY
)

# Count lambdas
nlam=$(ls "$mutdir/$leg"/lambda_* | wc -l | awk '{print $1}')
echo "Single-window baseline (1/GPU): ${ns_per_day} ns/day"
echo "Assumed windows/GPU: $win_per_gpu -> per-window rate: ${per_window_rate} ns/day"
echo "Per-window length: eq ${eq_ns} ns + prod ${prod_ns} ns = $(python3 - <<PY
print(f"{float("$eq_ns")+float("$prod_ns"):.3f}")
PY
) ns"
echo "Projected time per window at this sharing: ${time_per_window_hr} h"

tot_windows=$(python3 - <<PY
print(int("$nlam")*2)
PY
)
# With all windows launched concurrently, walltime ~ time_per_window_hr (plus overhead)
echo "Windows per mutation: ${tot_windows} (both legs). If run concurrently, expected walltime ~ ${time_per_window_hr} h per mutation (± overhead)."