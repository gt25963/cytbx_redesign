#!/usr/bin/env python3
"""
RQ2 Chai-1 input prep: one protein chain + retained HEM + swapped cofactor.
Usage: python prep_chai_fasta_rq2.py <input.fasta> <output.fasta> <cofactor>
  <cofactor> is FMN or U10 (selects the cofactor SMILES).
Chain order written: protein (A), HEM (B), cofactor (C) -> indices 0,1,2.
"""
import re, sys

HEM_SMILES = "CC1=C(CCC(O)=O)C2=[N]3C1=Cc1c(C)c(C=C)c4C=C5C(C)=C(C=C)C6=[N]5[Fe]3(n14)n1c(=C6)c(C)c(CCC(O)=O)c1=C2"
COFACTOR_SMILES = {
    "FMN": "Cc1cc2c(cc1C)N(C3=NC(=O)NC(=O)C3=N2)C[C@@H]([C@@H]([C@@H](COP(=O)(O)O)O)O)O",
    "U10": r"CC1=C(C(=O)C(=C(C1=O)OC)OC)C\C=C(/C)\CC\C=C(/C)\CC\C=C(/C)\CC\C=C(/C)\CC\C=C(/C)\CC\C=C(/C)\CC\C=C(/C)\CC\C=C(/C)\CCC=C(C)C",
}

def reformat(input_fasta, output_fasta, cofactor):
    if cofactor not in COFACTOR_SMILES:
        sys.exit(f"Unknown cofactor {cofactor}; expected FMN or U10")
    cof_smiles = COFACTOR_SMILES[cofactor]
    with open(input_fasta) as fin, open(output_fasta, "w") as fout:
        lines = [l.strip() for l in fin]
        i = 0
        while i < len(lines):
            line = lines[i]
            if line in ("--", ""):
                i += 1; continue
            if line.startswith(">"):
                m = re.search(r">([^,]+).*id=(\d+)", line)
                if m and i + 1 < len(lines):
                    pname, pid = m.group(1).strip(), m.group(2).strip()
                    seq = lines[i + 1].split(":")[0]   # single chain only
                    fout.write(f">protein|name={pname}_id{pid}_chain1\n{seq}\n")
                    fout.write(f">ligand|name=heme_b_retained\n{HEM_SMILES}\n")
                    fout.write(f">ligand|name={cofactor.lower()}_swapped\n{cof_smiles}\n")
                    i += 2; continue
                print(f"Skipping invalid header: {line}")
            i += 1
    print(f"RQ2 Chai input ({cofactor}) saved to {output_fasta}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: python prep_chai_fasta_rq2.py <input.fasta> <output.fasta> <FMN|U10>")
    reformat(sys.argv[1], sys.argv[2], sys.argv[3])
