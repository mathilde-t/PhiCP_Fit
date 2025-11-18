#!/bin/bash

# Prüfen, ob ein Ordner angegeben wurde
if [ -z "$1" ]; then
    echo "Fehler: Bitte gib einen Ordner an, z.B.:"
    echo "   ./rename_pickles.sh 04_DesyTau_phi_reorder_main_struc"
    exit 1
fi

# Ordner setzen
BASE="$1"

# Prüfen, ob der Ordner existiert
if [ ! -d "$BASE" ]; then
    echo "Fehler: Der Ordner $BASE existiert nicht!"
    exit 1
fi

echo "🔍 Suche Pickle-Dateien in: $BASE"
echo

find "$BASE" -type f -name "*.pickle" | while read -r f; do
    dir=$(dirname "$f")
    file=$(basename "$f")

    # Era extrahieren (preee / postee / prebpix / postbpix)
    if [[ "$file" =~ preee|postee|prebpix|postbpix ]]; then
        era=$(echo "$file" | grep -oE 'preee|postee|prebpix|postbpix')
    else
        echo "! Kein Era in $file – überspringe"
        continue
    fi

    # Alles nach letztem "hist__"
    tail=$(echo "$file" | sed 's/.*hist__/hist__/')

    # Doppelte hist__ korrigieren
    tail=$(echo "$tail" | sed 's/hist__hist__/hist__/')

    # Neues Format
    new="${era}_${tail}"

    echo "🔧 $file  →  $new"
    mv "$f" "$dir/$new"
done

echo
echo "Fertig!"

