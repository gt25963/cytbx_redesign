#!/usr/bin/env bash
set -euo pipefail
ligand_ccd="HEM"
outdir="boltz2_inputs"
chain_start=""
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
awk -v outdir="$outdir" -v ligand_ccd="$ligand_ccd" -v chain_start_req="$chain_start" '
  BEGIN {
    skip_block = 1
    in_block = 0
    yaml_count = 0
    total_chains = 0
    total_ligands = 0
    nitems = split(ligand_ccd, items, ",")
    for (i = 1; i <= nitems; i++) {
        split(items[i], parts, ":")
        lig_type[i] = parts[1]
        lig_count[i] = (parts[2] == "" ? 1 : parts[2])
    }
  }
  function write_yaml() {
    if (match(name, /^([^,]+),.*id=([0-9]+)/, arr)) {
        shortname = arr[1] "_id" arr[2]
    } else {
        shortname = "sequence_" ++yaml_count
    }
    file = outdir "/" shortname ".yaml"
    print "sequences:" > file
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
    auto_start_ascii = 64 + nchain + 1
    if (chain_start_req != "") {
        req_ascii = index("ABCDEFGHIJKLMNOPQRSTUVWXYZ", toupper(chain_start_req)) + 64
        if (req_ascii <= 64 + nchain) {
            print "WARNING: requested ligand chain start \"" chain_start_req \
                  "\" collides with protein chain letters (A-" sprintf("%c", 64+nchain) \
                  "); auto-bumping to " sprintf("%c", auto_start_ascii) > "/dev/stderr"
            curr_ascii = auto_start_ascii
        } else {
            curr_ascii = req_ascii
        }
    } else {
        curr_ascii = auto_start_ascii
    }
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
