#!/usr/bin/env python3
import argparse
import os

def main():
    p = argparse.ArgumentParser(description="Generate one resfile per peptide sequence.")
    p.add_argument("--peptides", required=True, help="Text file: one peptide per line")
    p.add_argument("--outdir", required=True, help="Output directory for .resfile files")
    p.add_argument("--chain", default="B", help="Chain letter to apply mutations to (default: B)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.peptides) as f:
        for pep in f:
            pep = pep.strip()
            if not pep:
                continue

            out_path = os.path.join(args.outdir, f"{pep}.resfile")
            with open(out_path, "w") as out:
                out.write("NATRO\nstart\n")
                for i, aa in enumerate(pep, start=1):
                    out.write(f"{i} {args.chain} PIKAA {aa}\n")

            print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
