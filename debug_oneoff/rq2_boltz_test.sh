#!/bin/bash
#SBATCH --account=brics.b5ae
#SBATCH --partition=workq
#SBATCH --time=01:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=rq2_boltz_test
#SBATCH --output=logs/rq2_boltz_test_%j.out
#Exit immediately on any error, undefined variable, or failed pipeline step
set -euo pipefail

#One-off test: pull a single designed sequence from the discontinued U10 pocket track and check how it scores in Boltz-2 before committing to a full-batch run
BASE=/scratch/b5ae/mvg2713124.b5ae/cytbx_pipeline
cd "$BASE"
POSE=poseB
DESIGN=1
FASTA=$BASE/rq2/design/U10_pocket/cycle1_relaxed_${POSE}/LigandMPNN/outputs/seqs/holo_hemC_relaxed.fa
WORKDIR=$BASE/rq2/design/U10_pocket/cycle1_relaxed_${POSE}/boltz_test
mkdir -p "$WORKDIR"

#Pull out just the one target design's sequence from the multi-sequence fasta
SEQ=$(awk -v id="$DESIGN" '
  /^>/ {keep=0; if ($0 ~ ("id="id",") || $0 ~ ("design_"id"$") || $0 ~ ("_"id"_")) keep=1; next}
  keep==1 {print; exit}
' "$FASTA")
echo "Sequence length: ${#SEQ}"

#Build a minimal Boltz-2 YAML input: protein chain, retained HEM, plus U10
YAML="$WORKDIR/u10_${POSE}_design${DESIGN}.yaml"
cat > "$YAML" << YEOF
version: 1
sequences:
  - protein:
      id: A
      sequence: ${SEQ}
  - ligand:
      id: B
      ccd: HEM
  - ligand:
      id: C
      ccd: U10
YEOF
echo "Wrote $YAML"

#Run Boltz-2 on this single test candidate
singularity run --nv /projects/b5ae/containers/boltz_venv.sif \
  boltz predict "$YAML" --out_dir "$WORKDIR/out" --use_msa_server false --override
