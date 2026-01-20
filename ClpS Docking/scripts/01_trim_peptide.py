#!/usr/bin/env python3
import argparse

def main():
    p = argparse.ArgumentParser(description="Trim peptide chain residues to a range; keep all other chains.")
    p.add_argument("--in", dest="inp", required=True, help="Input PDB")
    p.add_argument("--out", required=True, help="Output PDB")
    p.add_argument("--pep-chain", default="B", help="Peptide chain ID (default: B)")
    p.add_argument("--start", type=int, default=1, help="First residue number to keep on peptide chain")
    p.add_argument("--end", type=int, default=5, help="Last residue number to keep on peptide chain")
    args = p.parse_args()

    with open(args.inp) as fin, open(args.out, "w") as fout:
        for line in fin:
            if not line.startswith("ATOM"):
                fout.write(line)
                continue

            chain = line[21]
            resi = int(line[22:26])

            # Keep everything not on peptide chain; keep only range on peptide chain
            if chain != args.pep_chain or (args.start <= resi <= args.end):
                fout.write(line)

if __name__ == "__main__":
    main()
