#!/usr/bin/env bash
# Master RBFE pipeline script
# Usage: ./run_rbfe_pipeline.sh complex.pdb subs.txt
set -euo pipefail

# Input validation
pdb_file="${1:?Missing PDB file}"
subs_file="${2:?Missing mutations file (format: A33Y, one per line)}"
[ -f "$pdb_file" ] || { echo "PDB file not found: $pdb_file"; exit 1; }
[ -f "$subs_file" ] || { echo "Mutations file not found: $subs_file"; exit 1; }
command -v gmx >/dev/null 2>&1 || { echo "GROMACS (gmx) not found in PATH"; exit 1; }

# Directories setup
mkdir -p wt_complex/include wt_solvent/include mdp mutations logs
cp scripts/mdp/*.mdp mdp/

# Extract chains from PDB
echo "Extracting chains from $pdb_file..."
gmx editconf -f "$pdb_file" -o wt_complex/complex.pdb 2>/dev/null
cat > selection.txt << EOF
chain A
chain B
q
EOF
gmx make_ndx -f wt_complex/complex.pdb -o index.ndx < selection.txt &> logs/make_ndx.log || { echo "Failed to create index"; exit 1; }

# Extract nanobody (chain B) for solvent leg
echo 1 | gmx trjconv -s wt_complex/complex.pdb -f wt_complex/complex.pdb -o wt_solvent/nanobody.pdb -n index.ndx &> logs/extract_nb.log || { echo "Failed to extract nanobody"; exit 1; }

# Create resid map (needed for pmx)
echo "Creating residue mapping..."
python3 scripts/make_resid_map.py --pdb "$pdb_file" --out resid_map.tsv

# Prepare WT complex system
echo "Preparing WT complex system..."
cd wt_complex
echo -e "1\n1" | gmx pdb2gmx -f complex.pdb -o conf.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh &> ../logs/pdb2gmx_complex.log || { echo "Failed at pdb2gmx for complex"; exit 1; }
gmx editconf -f conf.gro -o box.gro -c -d 1.0 -bt dodecahedron &> ../logs/editconf_complex.log || { echo "Failed at editconf for complex"; exit 1; }
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top &> ../logs/solvate_complex.log || { echo "Failed at solvate for complex"; exit 1; }
gmx grompp -f ../mdp/em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 &> ../logs/grompp_ions_complex.log || { echo "Failed at grompp for ions"; exit 1; }
echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -neutral &> ../logs/genion_complex.log || { echo "Failed at genion for complex"; exit 1; }

# Minimization and equilibration
gmx grompp -f ../mdp/em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1 &> ../logs/grompp_em_complex.log || { echo "Failed at grompp for EM"; exit 1; }
gmx mdrun -v -deffnm em &> ../logs/mdrun_em_complex.log || { echo "Failed at mdrun for EM"; exit 1; }
gmx grompp -f ../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 &> ../logs/grompp_nvt_complex.log || { echo "Failed at grompp for NVT"; exit 1; }
gmx mdrun -v -deffnm nvt &> ../logs/mdrun_nvt_complex.log || { echo "Failed at mdrun for NVT"; exit 1; }
gmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 &> ../logs/grompp_npt_complex.log || { echo "Failed at grompp for NPT"; exit 1; }
gmx mdrun -v -deffnm npt &> ../logs/mdrun_npt_complex.log || { echo "Failed at mdrun for NPT"; exit 1; }

# Final production structure
cp npt.gro conf.gro
cd ..

# Prepare WT solvent system (nanobody only)
echo "Preparing WT solvent system..."
cd wt_solvent
echo -e "1\n1" | gmx pdb2gmx -f nanobody.pdb -o conf.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh &> ../logs/pdb2gmx_solvent.log || { echo "Failed at pdb2gmx for solvent"; exit 1; }
gmx editconf -f conf.gro -o box.gro -c -d 1.0 -bt dodecahedron &> ../logs/editconf_solvent.log || { echo "Failed at editconf for solvent"; exit 1; }
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top &> ../logs/solvate_solvent.log || { echo "Failed at solvate for solvent"; exit 1; }
gmx grompp -f ../mdp/em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1 &> ../logs/grompp_ions_solvent.log || { echo "Failed at grompp for ions"; exit 1; }
echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -neutral &> ../logs/genion_solvent.log || { echo "Failed at genion for solvent"; exit 1; }

# Minimization and equilibration
gmx grompp -f ../mdp/em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1 &> ../logs/grompp_em_solvent.log || { echo "Failed at grompp for EM"; exit 1; }
gmx mdrun -v -deffnm em &> ../logs/mdrun_em_solvent.log || { echo "Failed at mdrun for EM"; exit 1; }
gmx grompp -f ../mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 &> ../logs/grompp_nvt_solvent.log || { echo "Failed at grompp for NVT"; exit 1; }
gmx mdrun -v -deffnm nvt &> ../logs/mdrun_nvt_solvent.log || { echo "Failed at mdrun for NVT"; exit 1; }
gmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 &> ../logs/grompp_npt_solvent.log || { echo "Failed at grompp for NPT"; exit 1; }
gmx mdrun -v -deffnm npt &> ../logs/mdrun_npt_solvent.log || { echo "Failed at mdrun for NPT"; exit 1; }

# Final production structure
cp npt.gro conf.gro
cd ..

# Initialize project structure with lambda windows
echo "Initializing project structure..."
python3 scripts/init_project.py --mutations "$subs_file" --profile quick

# Generate pmx hybrid topologies
echo "Generating hybrid topologies with pmx..."
python3 scripts/pmx_make_hybrids.py --mutations "$subs_file" \
  --wt-complex wt_complex --wt-solvent wt_solvent \
  --resid-map resid_map.tsv --run

# Prepare all lambda windows
echo "Preparing all lambda windows..."
for mut_dir in mutations/*/; do
  mut=$(basename "$mut_dir")
  echo "Preparing $mut..."
  bash scripts/grompp_all.sh "mutations/$mut"
done

# Estimate runtime for first mutation
first_mut=$(head -n1 "$subs_file")
if [ -n "$first_mut" ]; then
  echo "Estimating runtime for $first_mut..."
  bash scripts/estimate_runtime.sh "mutations/$first_mut" complex lambda_00 6
fi

# Create master run script
cat > run_all_mutations.sh << 'EOF'
#!/usr/bin/env bash
set -euo pipefail
# Run all mutations sequentially
for mut_dir in mutations/*/; do
  mut=$(basename "$mut_dir")
  echo "Running $mut..."
  bash scripts/launch_rbfe_windows.sh "$mut_dir/windows.list"
  python3 scripts/analyze_ddg.py --root mutations
done
EOF
chmod +x run_all_mutations.sh

echo "============================================================="
echo "RBFE setup complete! To run all mutations:"
echo "./run_all_mutations.sh"
echo ""
echo "To run a specific mutation:"
echo "bash scripts/launch_rbfe_windows.sh mutations/MUTCODE/windows.list"
echo ""
echo "After completion, analyze results with:"
echo "python3 scripts/analyze_ddg.py --root mutations"
echo "============================================================="