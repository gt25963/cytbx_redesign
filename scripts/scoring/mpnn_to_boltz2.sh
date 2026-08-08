#!/usr/bin/env bash
#
# Convert LigandMPNN FASTA output into Boltz2 YAML input files.
#
# Usage:
#   ./fasta_to_boltz2.sh input.fasta [-l "HEM:2,NAG:1"] [-o boltz2_inputs] [-c B]
#
# Arguments:
#   input.fasta    FASTA file with LigandMPNN sequences
#
# Options:
#   -l LIGANDS     Ligands (comma-separated, optional counts with colon)
#                  e.g. "HEM", "HEM:2", "HEM:2,NAG:1"
#                  Default: HEM
#   -o OUTDIR      Output directory (default: boltz2_inputs)
#   -c CHAIN_START Chain ID letter to start ligands from (default: B)
#
# Example:
#   ./fasta_to_boltz2.sh ligandmpnn_output.fasta -l "HEM:2,NAG:1" -c L
#

set -euo pipefail

ligand_ccd="HEM"
outdir="boltz2_inputs"
chain_start="B"

if [ $# -lt 1 ]; then
  echo "Usage: $0 input.fasta [-l LIGANDS] [-o OUTDIR] [-c CHAIN_START]" >&2
  exit 1
fi

fasta="$1"
shift

while getopts ":l:o:c:" opt; do
  case "$opt" in
    l) ligand_ccd="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    c) chain_start="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done

mkdir -p "$outdir"

awk -v outdir="$outdir" -v ligand_ccd="$ligand_ccd" -v chain_start="$chain_start" '
  BEGIN {
    skip_block = 1
    in_block = 0
    yaml_count = 0
    total_chains = 0
    total_ligands = 0

    # Parse ligand specification
    nitems = split(ligand_ccd, items, ",")
    for (i = 1; i <= nitems; i++) {
        split(items[i], parts, ":")
        lig_type[i] = parts[1]
        lig_count[i] = (parts[2] == "" ? 1 : parts[2])
    }

    # Starting ligand chain letter
    if (length(chain_start) != 1 || chain_start !~ /^[A-Z]$/) {
      print "Error: invalid chain start" > "/dev/stderr"
      exit 1
    }
    start_ascii = index("ABCDEFGHIJKLMNOPQRSTUVWXYZ", toupper(chain_start)) + 64
  }

  function write_yaml() {
    # Short filename: protein name + id number
    if (match(name, /^([^,]+),.*id=([0-9]+)/, arr)) {
        shortname = arr[1] "_id" arr[2]
    } else {
        shortname = "sequence_" ++yaml_count
    }
    file = outdir "/" shortname ".yaml"

    print "sequences:" > file

    # Split protein chains
    if (seq ~ /:/) {
        nchain = split(seq, chains, ":")
    } else {
        nchain = 1
        chains[1] = seq
    }

    for (i=1; i <= nchain; i++) {
        cid = sprintf("[%c]", 64 + i)
        print "  - protein:" >> file
        print "      id: " cid >> file
      print "      sequence: " chains[i] >> file
      print "      msa: empty" >> file
        total_chains++
    }

    # Group ligands
    curr_ascii = start_ascii
    for (i=1; i <= nitems; i++) {
        ids = ""
        for (j=1; j <= lig_count[i]; j++) {
            lid = sprintf("%c", curr_ascii++)
            if (ids == "") ids = lid
            else ids = ids ", " lid
            total_ligands++
        }
        print "  - ligand:" >> file
        print "      id: [" ids "]" >> file
        print "      ccd: " lig_type[i] >> file
    }

    yaml_count++
  }

  /^>/ {
    if (in_block && skip_block == 0 && seq != "") write_yaml()

    name = $0
    sub(/^>/, "", name)
    in_block = 1

    # Skip first reference block
    if (skip_block > 0) {
        skip_block--
        in_block = 0
    }
    seq = ""
    next
  }

  /^[A-Z:]+$/ { if (in_block) seq = seq $0 }

  END {
    if (in_block && skip_block == 0 && seq != "") write_yaml()

    print "? Summary:"
    print "  YAML files generated: " yaml_count
    print "  Total protein chains: " total_chains
    print "  Total ligand entries: " total_ligands
  }
' "$fasta"

echo "Done! YAML files written to: $outdir"