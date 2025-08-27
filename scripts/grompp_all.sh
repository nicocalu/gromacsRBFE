#!/usr/bin/env bash
# Compile TPRs for all lambdas for a given mutation using existing hybrid topologies from pmx.
# Usage: bash scripts/grompp_all.sh mutations/A33Y
set -euo pipefail
mutdir="${1:?usage: grompp_all.sh mutations/<MUT>}"
for leg in complex solvent; do
  legdir="$mutdir/$leg"
  test -d "$legdir" || { echo "Missing $legdir"; exit 1; }
  # Expect hybrid/conf.gro and hybrid/topol.top
  test -f "$legdir/hybrid/conf.gro" || { echo "Missing $legdir/hybrid/conf.gro"; exit 1; }
  test -f "$legdir/hybrid/topol.top" || { echo "Missing $legdir/hybrid/topol.top"; exit 1; }
  for lamdir in "$legdir"/lambda_*; do
    pushd "$lamdir" >/dev/null
    # Equil
    gmx grompp -f eq.mdp -c ../hybrid/conf.gro -p ../hybrid/topol.top -o eq.tpr -maxwarn 1
    # Prod
    gmx grompp -f prod.mdp -c eq.tpr -t eq.cpt -p ../hybrid/topol.top -o prod.tpr -maxwarn 1
    popd >/dev/null
  done
done
echo "TPRs ready under $mutdir/{complex,solvent}/lambda_xx/"