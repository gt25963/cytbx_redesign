#!/usr/bin/env bash
set -euo pipefail
#Convert a multi-entry FASTA (LigandMPNN output) into per-design Boltz-2 YAML input files, one per sequence, with protein chains plus attached ligand(s)
ligand_ccd="HEM"
outdir="boltz2_inputs"
chain_start=""
if [ $# -lt 1 ]; then
  echo "Usage: $0 input.fasta [-l LIGANDS] [-o OUTDIR] [-c CHAIN_START]" >&2
  exit 1
fi
fasta="$1"
shift
#Parse optional flags: -l ligand CCD code(s)/counts, -o output dir, -c chain letter to start ligands at
while getopts ":l:o:c:" opt; do
  case "$opt" in
    l) ligand_ccd="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    c) chain_start="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done
mkdir -p "$outdir"
#Main conversion logic in awk: reads the fasta, groups each design's chain(s) and ligand(s) into one YAML block per candidate
awk -v outdir="$outdir" -v ligand_ccd="$ligand_ccd" -v chain_start_req="$chain_start" '
  BEGIN {
    skip_block = 1
    in_block = 0
    yaml_count = 0
    total_chains = 0
    total_ligands = 0
    #Parse the ligand spec (e.g. "HEM:2,FMN:1") into type/count pairs
    nitems = split(ligand_ccd, items, ",")
    for (i = 1; i <= nitems; i++) {
        split(items[i], parts, ":")
        lig_type[i] = parts[1]
        lig_count[i] = (parts[2] == "" ? 1 : parts[2])
    }
  }
  function write_yaml() {
    #Name the output file after the design id parsed from the fasta header
    if (match(name, /^([^,]+),.*id=([0-9]+)/, arr)) {
        shortname = arr[1] "_id" arr[2]
    } else {
        shortname = "sequence_" ++yaml_count
    }
    file = outdir "/" shortname ".yaml"
    print "sequences:" > file
    #Handle multi-chain sequences (colon-separated) for homo-oligomeric designs
    if (seq ~ /:/) {
        nchain = split(seq, chains, ":")
    } else {
        nchain = 1
        chains[1] = seq
    }
    #Write one protein block per chain, lettered A, B, C...; no MSA (matches the no-MSA AF3 setting used throughout, for speed on a de novo protein)
    for (i=1; i <= nchain; i++) {
        cid = sprintf("[%c]", 64 + i)
        print "  - protein:" >> file
        print "      id: " cid >> file
        print "      sequence: " chains[i] >> file
        print "      msa: empty" >> file
        total_chains++
    }
    #Ligand chain letters start right after the protein chains, unless a specific start letter was requested and doesn't collide with them
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
    #Write one ligand block per ligand type, assigning as many chain letters as that ligand's requested count (e.g. two HEM copies get two letters)
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
  #Detect a new fasta header, write out the previous block's YAML before starting the next
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
  #Accumulate sequence lines belonging to the current block
  /^[A-Z:]+$/ { if (in_block) seq = seq $0 }
  END {
    #Write the final block and print a summary
    if (in_block && skip_block == 0 && seq != "") write_yaml()
    print "? Summary:"
    print "  YAML files generated: " yaml_count
    print "  Total protein chains: " total_chains
    print "  Total ligand entries: " total_ligands
  }
' "$fasta"
echo "YAML files written to: $outdir"
