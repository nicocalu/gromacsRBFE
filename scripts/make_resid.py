#!/usr/bin/env python3
import argparse, re

def main():
    ap = argparse.ArgumentParser(description="Build chain+resid -> resid (int) map from WT PDB.")
    ap.add_argument("--pdb", required=True, help="WT complex PDB with chain IDs and correct residue numbering")
    ap.add_argument("--out", required=True, help="Output TSV: chain resid resname -> resid_int")
    args = ap.parse_args()

    # Parse PDB ATOM/HETATM lines
    # Columns: chainID at col 22, resSeq at cols 23-26 (1-based indexing)
    mapping = {}  # key: (chain, resseq) -> (resname, resid_int)
    with open(args.pdb) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 27:
                continue
            chain = line[21].strip() or " "
            resseq = line[22:26].strip()
            resname = line[17:20].strip()
            try:
                resid_int = int(resseq)
            except ValueError:
                # Skip insertion codes; you can pre-clean PDB if needed
                continue
            key = (chain, resid_int)
            if key not in mapping:
                mapping[key] = (resname, resid_int)

    with open(args.out, "w") as out:
        out.write("#chain\tresid_chain\tresname\tresid_gmx\n")
        for (chain, resid_chain), (resname, resid_gmx) in sorted(mapping.items(), key=lambda x: (x[0][0], x[0][1])):
            out.write(f"{chain}\t{resid_chain}\t{resname}\t{resid_gmx}\n")

    print(f"Wrote mapping for {len(mapping)} residues to {args.out}")

if __name__ == "__main__":
    main()