#!/usr/bin/env python3
import argparse, os, re, subprocess, sys, shutil, textwrap

MUT_RE = re.compile(r"^([A-Z])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")  # e.g., A33Y

def parse_mutations(path):
    muts=[]
    with open(path) as f:
        for line in f:
            s=line.strip()
            if not s or s.startswith("#"): continue
            m=MUT_RE.match(s)
            if not m:
                print(f"Skipping malformed mutation code: {s}", file=sys.stderr)
                continue
            chain, resid_chain, newaa = m.group(1), int(m.group(2)), m.group(3)
            muts.append((s, chain, resid_chain, newaa))
    return muts

def load_resid_map(path):
    # TSV: chain resid_chain resname resid_gmx
    mp={}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            parts=line.strip().split()
            if len(parts) < 4: continue
            chain, resid_chain, _, resid_gmx = parts[0], int(parts[1]), parts[2], int(parts[3])
            mp[(chain, resid_chain)] = resid_gmx
    return mp

def find_pmx_mutate_cmd():
    # Priority: env override
    env = os.getenv("PMX_MUTATE_CMD", "").strip()
    if env:
        return env.split()
    # Try common options
    candidates = [
        ["pmx", "mutate"],
        ["pmx-mutate"],
        [sys.executable, "-m", "pmx.scripts.mutate"],
    ]
    for cmd in candidates:
        if shutil.which(cmd[0]) or cmd[0].endswith("python") or cmd[0].endswith("python3"):
            return cmd
    return candidates[-1]

def build_mut_string(chain, resid_chain, resid_gmx, newaa, prefer_chain=True):
    # pmx typically accepts "A:33Y" if chains are present in input; if using .gro, chain info is absent,
    # so pass numeric resid only (e.g., "33Y"). We try chain-prefixed string by default.
    if prefer_chain:
        return f"{chain}:{resid_chain}{newaa}"
    else:
        return f"{resid_gmx}{newaa}"

def write_and_maybe_run(cmds, script_path, run=False, cwd=None):
    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    with open(script_path, "w") as f:
        f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        for c in cmds:
            f.write(c + "\n")
    os.chmod(script_path, 0o755)
    print(f"Wrote {script_path}")
    if run:
        print(f"Running {script_path} ...")
        subprocess.check_call([script_path], cwd=cwd)
        print("Done.")

def main():
    ap = argparse.ArgumentParser(description="Generate pmx hybrid topologies for RBFE from WT references.")
    ap.add_argument("--mutations", required=True, help="file with codes like A33Y (one per line)")
    ap.add_argument("--wt-complex", required=True, help="dir with conf.gro and topol.top for bound complex (WT)")
    ap.add_argument("--wt-solvent", required=True, help="dir with conf.gro and topol.top for solvent leg (WT)")
    ap.add_argument("--resid-map", default="", help="optional TSV from make_resid_map.py to map chain A/B resid -> GROMACS resid")
    ap.add_argument("--run", action="store_true", help="execute pmx mutate commands after writing scripts")
    args = ap.parse_args()

    for req in (os.path.join(args.wt_complex,"conf.gro"), os.path.join(args.wt_complex,"topol.top")):
        if not os.path.isfile(req):
            print(f"Missing WT complex file: {req}", file=sys.stderr); sys.exit(1)
    for req in (os.path.join(args.wt_solvent,"conf.gro"), os.path.join(args.wt_solvent,"topol.top")):
        if not os.path.isfile(req):
            print(f"Missing WT solvent file: {req}", file=sys.stderr); sys.exit(1)

    resid_map = load_resid_map(args.resid_map) if args.resid_map else {}
    prefer_chain = True if args.resid_map else False  # if we have mapping, we know chain IDs are correct; without it, prefer numeric
    pmx_cmd = find_pmx_mutate_cmd()
    print(f"Using pmx mutate command: {' '.join(pmx_cmd)}")

    muts = parse_mutations(args.mutations)
    if not muts:
        print("No valid mutations.", file=sys.stderr); sys.exit(1)

    for code, chain, resid_chain, newaa in muts:
        # Map to GROMACS resid if provided; else assume resid_chain equals resid in .gro
        resid_gmx = resid_map.get((chain, resid_chain), resid_chain)
        mut_string = build_mut_string(chain, resid_chain, resid_gmx, newaa, prefer_chain=prefer_chain)

        for leg, wtdir in (("complex", args.wt_complex), ("solvent", args.wt_solvent)):
            outdir = os.path.join("mutations", code, leg, "hybrid")
            os.makedirs(outdir, exist_ok=True)

            cmd = [
                *pmx_cmd,
                "--structure", os.path.join(wtdir, "conf.gro"),
                "--topology", os.path.join(wtdir, "topol.top"),
                "--mutation", mut_string,
                "--out", outdir
            ]
            # Some pmx installs expect different flag names; include alternates commented for quick edit.
            # e.g., "--ff", "amberff19SB", "--keepnames"
            shell_line = " ".join([subprocess.list2cmdline([c]) if " " in c else c for c in cmd])

            # Also copy include folder if present (so hybrid topol can #include it)
            inc_src = os.path.join(wtdir, "include")
            inc_dst = os.path.join(outdir, "include")
            copy_inc = ""
            if os.path.isdir(inc_src) and not os.path.isdir(inc_dst):
                copy_inc = f'cp -r "{inc_src}" "{inc_dst}"'

            script_path = os.path.join("mutations", code, f"pmx_mutate_{code}_{leg}.sh")
            lines = [shell_line]
            if copy_inc:
                lines.append(copy_inc)
            write_and_maybe_run(lines, script_path, run=args.run)

        # Convenience wrapper to run both legs
        wrapper = os.path.join("mutations", code, f"pmx_mutate_{code}_both.sh")
        combo = textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -euo pipefail
        bash "$(dirname "$0")/pmx_mutate_{code}_complex.sh"
        bash "$(dirname "$0")/pmx_mutate_{code}_solvent.sh"
        """)
        with open(wrapper, "w") as f:
            f.write(combo)
        os.chmod(wrapper, 0o755)
        print(f"Wrote {wrapper}")

    print("PMX mutation scripts generated. Review and run with --run or execute the per-mutation scripts.")

if __name__ == "__main__":
    main()