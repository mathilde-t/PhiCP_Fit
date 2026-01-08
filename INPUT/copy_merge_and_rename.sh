#!/bin/bash

# ---------------------------
# Options
# ---------------------------
DRYRUN=false
TARGET_SUBDIR="${1:-bugfree_samples_full_15_nov}"

# Check dry-run
if [[ "$1" == "--dry-run" ]]; then
    DRYRUN=true
    TARGET_SUBDIR="${2:-bugfree_samples_full_15_nov}"
fi

echo "Dry run: $DRYRUN"
echo "Target subdirectory: $TARGET_SUBDIR"

# ---------------------------
# Top-level config
# ---------------------------
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

# Local output base folder
OUTBASE="./"  # change if needed

# Helper function
copy_to_target() {
    local SRC="$1"
    local ERA="$2"
    local PROC="$3"
    local CHANNEL="$4"

    # Map PROC -> top folder
    if [[ "$PROC" == ggf* ]]; then
        TOP="ggF"
    else
        TOP="VBF"
    fi

    # CHANNEL -> second folder
    case "$CHANNEL" in
        cpo) SUBF="cpo" ;;
        mm)  SUBF="mm"  ;;
        sm)  SUBF="sm"  ;;
        *)   echo "Unknown channel $CHANNEL"; return ;;
    esac

    # Prepare target folder
    TARGETDIR="$OUTBASE/$TOP/$SUBF"
    mkdir -p "$TARGETDIR"

    # Construct new name (ERA + tail)
    BASENAME=$(basename "$SRC")
    TAIL=$(echo "$BASENAME" | sed 's/.*hist__/hist__/; s/hist__hist__/hist__/')
    NEWNAME="${ERA}_${TAIL}"

    # Copy
    if $DRYRUN; then
        echo "[DRY RUN] $SRC  ->  $TARGETDIR/$NEWNAME"
    else
        cp "$SRC" "$TARGETDIR/$NEWNAME"
        echo "Copied: $TARGETDIR/$NEWNAME"
    fi
}

# ---------------------------
# Main loop
# ---------------------------
for CONF in "${CONFS[@]}"; do
    echo "Processing $CONF"
    for DATASET in $CONF/*/; do
        DATASET=$(basename "$DATASET")

        MATCH=false
        for PROC in "${PROCS[@]}"; do
            if [[ "$DATASET" == *"$PROC"* ]]; then
                MATCH=true
                break
            fi
        done
        [[ "$MATCH" == false ]] && continue

        TARGET="$CONF/$DATASET/nominal/calib__main/sel__main/red__cf_default/prod__main/hist__httcp_hist_producer/$TARGET_SUBDIR"
        [[ ! -d "$TARGET" ]] && continue

        # Determine ERA
        ERA=$(echo "$CONF" | tr '[:upper:]' '[:lower:]' \
            | sed 's/run3_2022_postEE/postee/; s/run3_2022_preEE/preee/; s/run3_2023_postBPix/postbpix/; s/run3_2023_preBPix/prebpix/')

        # PROCESS (ggf/vbf) and CHANNEL (cpo/mm/sm)
        PROC=$(echo "$DATASET" | awk -F'_' '{print $2"_"$3}' | tr '[:upper:]' '[:lower:]')
        CHANNEL=$(echo "$DATASET" | awk -F'_' '{print $4}' | tr '[:upper:]' '[:lower:]')

        # Copy all pickle files
        for FILE in "$TARGET"/*.pickle; do
            [ -e "$FILE" ] || continue
            copy_to_target "$FILE" "$ERA" "$PROC" "$CHANNEL"
        done
    done
done

echo "Done!"

