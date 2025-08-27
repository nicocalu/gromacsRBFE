#!/usr/bin/env python3
import argparse, os, glob, math
import pandas as pd
from alchemlyb.parsing.gmx import load_dhdl
from alchemlyb.preprocessing import extract_u_nk, decorrelate_u_nk
from alchemlyb.estimators import MBAR

R_kJ = 8.314462618e-3  # kJ/mol/K

def gather_leg(leg_root):
    # Expect prod.dhdl.xvg in each lambda dir
    files = sorted(glob.glob(os.path.join(leg_root, "lambda_*", "prod.dhdl.xvg")))
    if not files:
        return None
    dfs = [load_dhdl(f) for f in files]
    # concat preserving metadata
    data = pd.concat(dfs, axis=0, sort=False)
    return data

def leg_deltaG_kJ(leg_root, T=300.0):
    data = gather_leg(leg_root)
    if data is None or data.empty:
        return None, None
    u_nk = extract_u_nk(data, T=T)
    u_nk = decorrelate_u_nk(u_nk)
    est = MBAR().fit(u_nk)
    # free energy difference from first to last state (kT)
    df_kT = est.delta_f_.iloc[0, -1]
    ddf_kT = est.d_delta_f_.iloc[0, -1]
    kT = R_kJ * T
    return df_kT * kT, ddf_kT * kT

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root", default="mutations", help="root directory containing mutation subdirs")
    args=ap.parse_args()
    rows=[]
    for mut in sorted(os.listdir(args.root)):
        mdir=os.path.join(args.root, mut)
        if not os.path.isdir(mdir): continue
        leg_c=os.path.join(mdir,"complex")
        leg_s=os.path.join(mdir,"solvent")
        if not (os.path.isdir(leg_c) and os.path.isdir(leg_s)): continue
        Gc, dGc = leg_deltaG_kJ(leg_c)
        Gs, dGs = leg_deltaG_kJ(leg_s)
        if Gc is None or Gs is None:
            continue
        ddG = Gc - Gs
        dddG = math.sqrt((dGc or 0)**2 + (dGs or 0)**2)
        rows.append({
            "mutation": mut,
            "DeltaG_complex_kJmol": Gc,
            "SE_complex_kJmol": dGc,
            "DeltaG_solvent_kJmol": Gs,
            "SE_solvent_kJmol": dGs,
            "DeltaDeltaG_kJmol": ddG,
            "SE_DeltaDeltaG_kJmol": dddG,
            "DeltaDeltaG_kcalmol": ddG/4.184,
            "SE_DeltaDeltaG_kcalmol": dddG/4.184
        })
    if not rows:
        print("No results found.")
        return
    df=pd.DataFrame(rows).sort_values("DeltaDeltaG_kJmol")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(df.to_string(index=False, float_format=lambda x: f"{x: .3f}"))

if __name__=="__main__":
    main()