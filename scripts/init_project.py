#!/usr/bin/env python3
import argparse, os, shutil, re, sys, textwrap

QUICK = {
  "name": "quick",
  "eq_steps": 125000,  # 0.25 ns @ 2 fs
  "prod_steps": 500000,  # 1.0 ns
  "coul": [0.00, 0.25, 0.50, 0.75, 1.00] + [1.0]*11,
  "vdw":  [0,0,0,0,0, 0.10,0.20,0.30,0.45,0.60,0.75,0.88,0.94,0.97,0.99,1.00],
}
STANDARD = {
  "name": "standard",
  "eq_steps": 250000,  # 0.5 ns
  "prod_steps": 1000000,  # 2.0 ns
  "coul": [0.00,0.20,0.40,0.60,0.80,1.00] + [1.0]*18,
  "vdw":  [0,0,0,0,0,0, 0.05,0.10,0.15,0.22,0.30,0.40,0.52,0.65,0.78,0.88,0.94,0.97,1.00,1.00,1.00,1.00,1.00,1.00],
}

MDP_BASE = {
"common": textwrap.dedent("""\
integrator               = md
dt                       = 0.002
nstxout-compressed       = 0
nstenergy                = 2500
nstlog                   = 2500
nstcalcenergy            = 50
cutoff-scheme            = Verlet
rcoulomb                 = 1.0
rvdw                     = 1.0
coulombtype              = PME
pme-grid-spacing         = 0.12
vdwtype                  = PME
DispCorr                 = no
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau-p                    = 5.0
ref-p                    = 1.0
compressibility          = 4.5e-5
constraints              = h-bonds
constraint-algorithm     = lincs
lincs-order              = 4
lincs-iter               = 1
free-energy              = yes
sc-algorithm             = gmx
sc-power                 = 1
sc-alpha                 = 0.5
sc-sigma                 = 0.3
calc-lambda-neighbors    = 3
separate-dhdl-file       = yes
dhdl-print-energy        = yes
dhdl-derivatives         = yes
dhdl-output-interval     = 500
"""),
"eq_header": "",
"prod_header": textwrap.dedent("""\
; lambdas set below
""")
}

EQ_HEADER = "init-lambda-state        = {LAMBDA_INDEX}\n"
PROD_TAIL_TMPL = """coul-lambdas             = {COUL}
vdw-lambdas              = {VDW}
bonded-lambdas           = {VDW}
mass-lambdas             = {VDW}
"""

EM_MDP = """integrator = steep
nsteps = 50000
emtol = 1000
cutoff-scheme = Verlet
rcoulomb = 1.0
rvdw = 1.0
coulombtype = PME
pme-grid-spacing = 0.12
vdwtype = PME
DispCorr = no
constraints = none
nstlog = 500
nstenergy = 500
"""

NVT_MDP = """integrator = md
dt = 0.002
nsteps = 50000
tcoupl = V-rescale
tc-grps = System
tau-t = 1.0
ref-t = 300
pcoupl = no
cutoff-scheme = Verlet
rcoulomb = 1.0
rvdw = 1.0
coulombtype = PME
pme-grid-spacing = 0.12
vdwtype = PME
DispCorr = no
constraints = h-bonds
constraint-algorithm = lincs
lincs-order = 4
lincs-iter = 1
"""

NPT_MDP = """integrator = md
dt = 0.002
nsteps = 250000
tcoupl = V-rescale
tc-grps = System
tau-t = 1.0
ref-t = 300
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 5.0
ref-p = 1.0
compressibility = 4.5e-5
cutoff-scheme = Verlet
rcoulomb = 1.0
rvdw = 1.0
coulombtype = PME
pme-grid-spacing = 0.12
vdwtype = PME
DispCorr = no
constraints = h-bonds
constraint-algorithm = lincs
lincs-order = 4
lincs-iter = 1
"""

def write_file(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(content)

def parse_mutations(mfile):
    muts=[]
    rx=re.compile(r"^([A-Z])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")
    with open(mfile) as f:
        for line in f:
            s=line.strip()
            if not s or s.startswith("#"): continue
            m=rx.match(s)
            if not m:
                print(f"Skipping malformed line: {s}", file=sys.stderr)
                continue
            chain, resi, newaa = m.group(1), int(m.group(2)), m.group(3)
            muts.append((s, chain, resi, newaa))
    return muts

def init_mdp_templates(profile):
    # base mdp dir
    write_file("mdp/em.mdp", EM_MDP)
    write_file("mdp/nvt.mdp", NVT_MDP)
    write_file("mdp/npt.mdp", NPT_MDP)
    # rbfe eq/prod base without lambda arrays
    for kind in ("eq","prod"):
        base = MDP_BASE["common"]
        if kind=="eq":
            base += EQ_HEADER.format(LAMBDA_INDEX="{LIDX}")
            steps = profile["eq_steps"]
            write_file("mdp/rbfe_eq.base.mdp", base + f"nsteps = {steps}\n")
        else:
            steps = profile["prod_steps"]
            write_file("mdp/rbfe_prod.base.mdp", base + f"nsteps = {steps}\n" + MDP_BASE["prod_header"])

def make_lambda_dirs(mutdir, lambdas):
    L = len(lambdas["vdw"])
    for i in range(L):
        lname=f"lambda_{i:02d}"
        for leg in ("complex","solvent"):
            d=os.path.join(mutdir,leg,lname)
            os.makedirs(d, exist_ok=True)
            # materialize eq/prod mdp with proper init-lambda-state and arrays
            with open("mdp/rbfe_eq.base.mdp") as f: eq_base=f.read()
            eq=eq_base.replace("{LIDX}", str(i))
            write_file(os.path.join(d,"eq.mdp"), eq)
            with open("mdp/rbfe_prod.base.mdp") as f: prod_base=f.read()
            prod = prod_base + PROD_TAIL_TMPL.format(
                COUL=" ".join([("{:.2f}".format(x) if isinstance(x,float) else str(x)) for x in lambdas["coul"]]),
                VDW=" ".join([("{:.2f}".format(x) if isinstance(x,float) else str(x)) for x in lambdas["vdw"]]),
            )
            write_file(os.path.join(d,"prod.mdp"), prod)

def write_windows_list(mutdir, gpus=8, per_gpu=6):
    # Assign windows round-robin across GPUs; both legs in same list
    lambdas = sorted([d for d in os.listdir(os.path.join(mutdir,"complex")) if d.startswith("lambda_")])
    entries=[]
    gpu_ids=list(range(gpus))
    idx=0
    # interleave complex and solvent windows for load balance
    for i, ln in enumerate(lambdas):
        for leg in ("complex","solvent"):
            entries.append((os.path.join(mutdir,leg,ln), gpu_ids[idx % len(gpu_ids)]))
            idx+=1
    # limit concurrency to gpus*per_gpu lines when you run; rest can be queued or run sequentially
    with open(os.path.join(mutdir,"windows.list"),"w") as f:
        for d,g in entries:
            f.write(f"{d} {g}\n")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mutations", required=True, help="file with lines like A33Y, B57F")
    ap.add_argument("--profile", choices=["quick","standard"], default="quick")
    args=ap.parse_args()
    profile = QUICK if args.profile=="quick" else STANDARD
    lambdas={"coul": profile["coul"], "vdw": profile["vdw"]}
    init_mdp_templates(profile)
    muts = parse_mutations(args.mutations)
    if not muts:
        print("No valid mutations parsed.", file=sys.stderr)
        sys.exit(1)
    for mcode, chain, resi, newaa in muts:
        mutdir = os.path.join("mutations", mcode)
        os.makedirs(mutdir, exist_ok=True)
        # drop a note for pmx
        write_file(os.path.join(mutdir,"PMX_MUTATION.txt"), f"{chain}:{resi}{newaa}\n")
        make_lambda_dirs(mutdir, lambdas)
        write_windows_list(mutdir, gpus=8, per_gpu=6)
    print(f"Initialized {len(muts)} mutation(s) with profile '{profile['name']}'.")
    print("Next: generate hybrids with pmx per mutation, then:")
    print("  bash scripts/grompp_all.sh mutations/<MUT>")
    print("  bash scripts/launch_rbfe_windows.sh mutations/<MUT>/windows.list")

if __name__=="__main__":
    main()