#!/usr/bin/env python3

# Q8 Chai input prep: one protein chain + BOTH retained hemes + surface Q8.
# Q8 does NOT replace a haem (unlike FMN); both hemes are kept.
# Chain order: protein (A), HEM_B (B), HEM_C (C), Q8 (D) -> protein<->Q8 = [0][3].

import re, sys
HEM_SMILES = "CC1=C(CCC(O)=O)C2=[N]3C1=Cc1c(C)c(C=C)c4C=C5C(C)=C(C=C)C6=[N]5[Fe]3(n14)n1c(=C6)c(C)c(CCC(O)=O)c1=C2"
Q8_SMILES = r"COC1=C(OC)C(=O)C(=C(C)C1=O)C\C=C(C)\CC\C=C(C)\CC\C=C(C)\CC\C=C(C)\CC\C=C(C)\CC\C=C(C)\CC\C=C(C)\CCC=C(C)C"
def reformat(input_fasta, output_fasta):
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
                    seq = lines[i + 1].split(":")[0] ## single chain only 
                    fout.write(f">protein|name={pname}_id{pid}_chain1\n{seq}\n")
                    fout.write(f">ligand|name=heme_b_retained\n{HEM_SMILES}\n")
                    fout.write(f">ligand|name=heme_c_retained\n{HEM_SMILES}\n")
                    fout.write(f">ligand|name=q8_surface\n{Q8_SMILES}\n")
                    i += 2; continue
                print(f"Skipping invalid header: {line}")
            i += 1
    print(f"Q8 Chai input saved to {output_fasta}")
if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python prep_chai_fasta_q8.py <input.fasta> <output.fasta>")
    reformat(sys.argv[1], sys.argv[2])
