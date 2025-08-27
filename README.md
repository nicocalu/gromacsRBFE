RBFE scaffold for Nanobody–Antigen (AMBER ff19SB + TIP3P, LJ-PME)

What you need prepared once (WT references)
- wt_complex/: equilibrated WT bound complex (GROMACS files: conf.gro, topol.top, include/itp as needed).
- wt_solvent/: equilibrated WT nanobody (or antigen, whichever leg you’re transforming in solvent) in a box with ions to neutralize.
- Protonation fixed (His tautomers set), same across legs.

Install
- GROMACS 2025.2 with CUDA offload (PME/bonded/update).
- Python: pip install -r requirements.txt
- Enable CUDA MPS on the node once: nvidia-cuda-mps-control -d

Workflow (first tests, 5–10 mutations)
1) Create mutation list
   - mutations.txt with lines like:
     A33Y
     A52F
     B76Y

2) Initialize scaffold
   - Quick screen (16 λ: 0.25 ns eq + 1.0 ns prod):
     python3 scripts/init_project.py --mutations mutations.txt --profile quick
   - Standard screen (24 λ: 0.5 ns eq + 2.0 ns prod):
     python3 scripts/init_project.py --mutations mutations.txt --profile standard

   This creates:
   - mutations/<MUT>/{complex,solvent}/lambda_xx/
   - mdp/ with templates, and per-window prod/eq mdp copies.
   - windows lists for each mutation in mutations/<MUT>/windows.list (8 GPUs, round-robin).

3) Generate hybrid topologies with pmx (one-time per mutation and leg)
   - Use your equilibrated WT reference (wt_complex, wt_solvent) and run pmx to create the hybrid:
     Example (pseudocode; adapt to your pmx setup):
       pmx mutate --structure wt_complex/conf.gro --topology wt_complex/topol.top \
         --mutation "A:33Y" --out mutations/A33Y/complex/hybrid
       pmx mutate --structure wt_solvent/conf.gro --topology wt_solvent/topol.top \
         --mutation "A:33Y" --out mutations/A33Y/solvent/hybrid

   - Ensure each leg directory has:
       hybrid/conf.gro
       hybrid/topol.top
       hybrid/include/* (if itps are used)

   - Then compile TPRs for all λ for that mutation:
       bash scripts/grompp_all.sh mutations/A33Y
     (This picks eq.mdp/prod.mdp in each lambda dir and uses init-lambda-state accordingly.)

4) Benchmark and plan walltime (optional but recommended)
   - bash scripts/estimate_runtime.sh mutations/A33Y complex lambda_00 6
     This runs a 200 ps micro-benchmark in that window, estimates ns/day at 6 windows/GPU, and prints a walltime estimate for the current profile.

5) Run the mutation
   - Launch windows for one mutation (fills 48 windows by default = 6/GPU):
     bash scripts/launch_rbfe_windows.sh mutations/A33Y/windows.list
   - Repeat for other mutations in sequence (or run two simultaneously if CPU cores permit).

6) Analyze ΔΔG
   - python3 scripts/analyze_ddg.py --root mutations
     Outputs a succinct table: mutation, ΔG_complex, ΔG_solvent, ΔΔG (kJ/mol and kcal/mol).

Notes
- Keep all RBFE starting structures from a single equilibrated WT reference (do NOT use per-mutation docked starts for RBFE).
- For charge-changing mutations, start with neutral-only in first tests. We can add alchemical counter-ions later.
- If any λ pair shows poor dhdl overlap (especially near end states), add an extra vdW λ near 0.95–1.0 and rerun those windows.

GPU/CPU
- Default launcher uses: -ntmpi 1 -ntomp 1 -pin on -nb gpu -pme gpu -bonded gpu -update gpu
- With 8× A100 and 64 CPU cores, 48 windows concurrent (6/GPU) is a good starting point.
