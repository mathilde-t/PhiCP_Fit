#!/bin/bash
'''
Author : Mathilde Witt
Description : This script renames the pickle files obtaiens from 
              phi_cp_merges_hists.sh which are presorted by the user 
              and puts them in  the same strcture but with shortened 
              names needed for the  PhiCPfit_BarPlot_DESYTau_cf_0p3 
              script aka plotting of PhiCP histogramms.
'''
# ---------------------------
# Options
# ---------------------------
# --dry-run:      Simulate actions without copying files
# First argument: Target subdirectory aka version name from cf 
#                 (default: bugfree_samples_full_15_nov)
DRYRUN=false
TARGET_SUBDIR="${1:-bugfree_samples_full_15_nov}"

# Check if --dry-run flag is set
if [[ "$1" == "--dry-run" ]]; then
    DRYRUN=true
    TARGET_SUBDIR="${2:-bugfree_samples_full_15_nov}"
fi

echo "Dry run: $DRYRUN"
echo "Target subdirectory: $TARGET_SUBDIR"

# ---------------------------
# Top-level Configurations
# ---------------------------
# List of era configurations to process
CONFS=(
    run3_2022_preEE_mutau
    run3_2022_postEE_mutau
    run3_2023_preBPix_mutau
    run3_2023_postBPix_mutau
)

# List of processes to match in datasets
PROCS=(
    h_ggf_htt_sm_prod_sm_filtered
    h_ggf_htt_mm_prod_sm_filtered
    h_ggf_htt_cpo_prod_sm_filtered
    h_vbf_htt_sm_filtered
    h_vbf_htt_mm_filtered
    h_vbf_htt_cpo_filtered
)

# Base output directory for copied files
OUTBASE="./" 

# ---------------------------
# Helper Function: Copy and Rename Files
# ---------------------------
# Args:
#   $1: Source file path
#   $2: Era identifier
#   $3: Process type (ggf/vbf)
#   $4: Channel type (cpo/mm/sm)
copy_to_target() {
    local SRC="$1"
    local ERA="$2"
    local PROC="$3"
    local CHANNEL="$4"

    # Map process type to top-level folder (ggF/VBF)
    if [[ "$PROC" == ggf* ]]; then
        TOP="ggF"
    else
        TOP="VBF"
    fi

    # Map channel type to subfolder (cpo/mm/sm)
    case "$CHANNEL" in
        cpo) SUBF="cpo" ;;
        mm)  SUBF="mm"  ;;
        sm)  SUBF="sm"  ;;
        *)   echo "Unknown channel $CHANNEL"; return ;;
    esac

    # Create target directory if it doesn't exist
    TARGETDIR="$OUTBASE/$TOP/$SUBF"
    mkdir -p "$TARGETDIR"

    # Extract and rename file: ERA + original tail
    BASENAME=$(basename "$SRC")
    TAIL=$(echo "$BASENAME" | sed 's/.*hist__/hist__/; s/hist__hist__/hist__/')
    NEWNAME="${ERA}_${TAIL}"

    # Copy file (or simulate if dry-run)
    if $DRYRUN; then
        echo "[DRY RUN] $SRC  ->  $TARGETDIR/$NEWNAME"
    else
        cp "$SRC" "$TARGETDIR/$NEWNAME"
        echo "Copied: $TARGETDIR/$NEWNAME"
    fi
}

# ---------------------------
# Main Processing Loop
# ---------------------------
for CONF in "${CONFS[@]}"; do
    echo "Processing $CONF"
    for DATASET in $CONF/*/; do
        DATASET=$(basename "$DATASET")

        # Check if dataset matches any process
        MATCH=false
        for PROC in "${PROCS[@]}"; do
            if [[ "$DATASET" == *"$PROC"* ]]; then
                MATCH=true
                break
            fi
        done
        [[ "$MATCH" == false ]] && continue

        # Construct target path for pickle files
        TARGET="$CONF/$DATASET/nominal/calib__main/sel__main/red__cf_default/prod__main/hist__httcp_hist_producer/$TARGET_SUBDIR"
        [[ ! -d "$TARGET" ]] && continue

        # Extract era, process, and channel from dataset name
        ERA=$(echo "$CONF" | tr '[:upper:]' '[:lower:]' \
            | sed 's/run3_2022_postEE/postee/; s/run3_2022_preEE/preee/; s/run3_2023_postBPix/postbpix/; s/run3_2023_preBPix/prebpix/')
        PROC=$(echo "$DATASET" | awk -F'_' '{print $2"_"$3}' | tr '[:upper:]' '[:lower:]')
        CHANNEL=$(echo "$DATASET" | awk -F'_' '{print $4}' | tr '[:upper:]' '[:lower:]')

        # Copy and rename all pickle files in target directory
        for FILE in "$TARGET"/*.pickle; do
            [ -e "$FILE" ] || continue
            copy_to_target "$FILE" "$ERA" "$PROC" "$CHANNEL"
        done
    done
done

echo "Done!"