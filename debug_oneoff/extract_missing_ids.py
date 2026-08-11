#!/usr/bin/env python3
"""
Extract specific sequence-ID blocks from a Chai-1 multi-entry FASTA.

Each sequence ID forms a "block" of paired header/sequence lines:
  >protein|name=top_scoring.cif_id<N>_chain1
  <seq>
  >ligand|name=heme_b_<k>
  <smiles>
  >ligand|name=heme_b_<k+1>
  <smiles>
  >protein|name=top_scoring.cif_id<N>_chain2
  <seq>
  >ligand|name=heme_b_<k+2>
  <smiles>
  >ligand|name=heme_b_<k+3>
  <smiles>

Usage:
    python extract_missing_ids.py <input_fasta> <id1,id2,id3,...> <output_fasta>
"""
#recovery script: pulls out only candidates listed by id from a larger chai input fasta, for re-running/checking specific designs without resubmitting the whole batch
import sys
import re

def parse_fasta_records(path):
    # read fasta into pairs (header, seqeunce), one per record
    records = []
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "\n".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, "\n".join(seq_lines)))
    return records


def group_into_blocks(records):
    # group records into per-candidate blocks, starting a new block whenever a chain1 protein header appears (marks start of new desing)
    blocks = []
    current_block = []
    for header, seq in records:
        is_chain1_start = bool(re.search(r"protein\|name=.*_chain1$", header))
        if is_chain1_start and current_block:
            blocks.append(current_block)
            current_block = []
        current_block.append((header, seq))
    if current_block:
        blocks.append(current_block)
    return blocks


def block_id(block):
    # pull numeric desing id out of a block's chain1 header
    header = block[0][0]
    m = re.search(r"_id(\d+)_chain1$", header)
    if not m:
        return None
    return int(m.group(1))


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    input_fasta, id_list_str, output_fasta = sys.argv[1:4]
    target_ids = set(int(x) for x in id_list_str.split(","))

    records = parse_fasta_records(input_fasta)
    blocks = group_into_blocks(records)

    print(f"Parsed {len(records)} FASTA records into {len(blocks)} blocks")

    #keep ounly the blocks matching a requested target id
    found_ids = set()
    output_records = []
    for block in blocks:
        bid = block_id(block)
        if bid in target_ids:
            output_records.extend(block)
            found_ids.add(bid)

    #flag any requestedd ids that weren't found at all, so nothing goes missing silenetly 
    missing_from_input = target_ids - found_ids
    if missing_from_input:
        print(f"WARNING: these target IDs were not found in the input fasta at all: {sorted(missing_from_input)}")

    with open(output_fasta, "w") as f:
        for header, seq in output_records:
            f.write(header + "\n")
            f.write(seq + "\n")

    print(f"Extracted {len(found_ids)} blocks ({len(output_records)} fasta records) to {output_fasta}")
    print(f"IDs extracted: {sorted(found_ids)}")


if __name__ == "__main__":
    main()
