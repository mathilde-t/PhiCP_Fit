#!/bin/bash

###########################################
# Usage instructions
#
# Run normally (copies files) :
#     ./phi_cp_merges_hists.sh
#     ./phi_cp_merges_hists.sh <target_subdirectory>
#
# Run in dry-run mode (no files are actually copied) :
#     ./phi_cp_merges_hists.sh --dry-run
#     ./phi_cp_merges_hists.sh --dry-run <target_subdirectory>
#
# Notes:
# - The <target_subdirectory> refers to the folder inside each dataset
#   from which pickle files are collected (default: bugfree_samples_full_15_nov)
# - Dry-run mode prints all actions without performing them
# - The output directory is fixed at /eos/user/m/mwitt/phi_cp_merges_hists
###########################################


###########################################
# Options
###########################################
DRYRUN=false

# Check whether --dry-run was passed
if [[ "$1" == "--dry-run" ]]; then
    DRYRUN=true
    # The actual argument (TARGET_SUBDIR) shifts by one position
    TARGET_SUBDIR="${2:-bugfree_samples_full_15_nov}"
else
    TARGET_SUBDIR="${1:-bugfree_samples_full_15_nov}"
fi

# Small helper function for copy commands
copy() {
    local SRC="$1"
    local DST="$2"

    if $DRYRUN; then
        echo "[DRY RUN] would copy: $SRC  -->  $DST"
    else
        cp "$SRC" "$DST"
        echo "        copied: $(basename "$DST")"
    fi
}

###########################################
# Start
###########################################
OUTDIR="/eos/user/m/mwitt/phi_cp_merges_hists"
mkdir -p "$OUTDIR"

echo " --> Dry run: $DRYRUN"
echo " --> Using target subdirectory: $TARGET_SUBDIR"
echo ""

CONFS=(
    run3_2022_preEE_mutau
    run3_2022_postEE_mutau
    run3_2023_preBPix_mutau
    run3_2023_postBPix_mutau
)
PROCS=(
    h_ggf_htt_sm_prod_sm_filtered
    h_ggf_htt_mm_prod_sm_filtered
    h_ggf_htt_cpo_prod_sm_filtered
    h_vbf_htt_sm_filtered
    h_vbf_htt_mm_filtered
    h_vbf_htt_cpo_filtered
)

echo "Starting copy operation…"
echo ""

for CONF in "${CONFS[@]}"; do
    echo "===> Processing $CONF"

    for DATASET in $CONF/*/; do
        DATASET=$(basename "$DATASET")

        MATCH_FOUND=false
        for PROC in "${PROCS[@]}"; do
            if [[ "$DATASET" == *"$PROC"* ]]; then
                MATCH_FOUND=true
                break
            fi
        done

        # Skip datasets that do not match any process
        if [[ "$MATCH_FOUND" == false ]]; then
            echo "     (skip: $DATASET does not correspond to any PROC)"
            continue
        fi

        TARGET="$CONF/$DATASET/nominal/calib__main/sel__main/red__cf_default/prod__main/hist__httcp_hist_producer/$TARGET_SUBDIR"

        echo "     Checking: $TARGET"
        if [[ ! -d "$TARGET" ]]; then
            echo "     (skip: no $TARGET_SUBDIR found in $DATASET)"
            continue
        fi

        echo "     -> Relevant process identified, searching for pickle files…"

        ERA=$(echo "$CONF" | tr '[:upper:]' '[:lower:]' \
            | sed 's/run3_2022_postEE/postee/; s/run3_2022_preEE/preee/; s/run3_2023_postBPix/postbpix/; s/run3_2023_preBPix/prebpix/')

        PROCESS=$(echo "$DATASET" | awk -F'_' '{print $2"_"$3}' | tr '[:upper:]' '[:lower:]')
        CHANNEL=$(echo "$DATASET" | awk -F'_' '{print $4}' | tr '[:upper:]' '[:lower:]')

        # Loop over all pickle files
        for FILE in "$TARGET"/*.pickle; do
            [ -e "$FILE" ] || continue
            BASENAME=$(basename "$FILE")

            # Keep only files containing 'phi_cp_mu_a1'
            if [[ "$BASENAME" != *"phi_cp_mu_a1"* ]]; then
                continue
            fi

            # Construct the new filename
            NEWNAME="${ERA}_${PROCESS}_${CHANNEL}_hist__${BASENAME}"

            copy "$FILE" "$OUTDIR/$NEWNAME"
        done
    done
    echo ""
done

echo "Finished!"
