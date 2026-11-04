#!/bin/bash

# Itera su tutti i file che contengono .mat seguito da uno spazio e altri caratteri
for f in *.mat\ *; do
    # Verifica che il file esista effettivamente
    if [ -f "$f" ]; then
        # Taglia tutto a partire da .mat in poi e riaggiunge .mat
        newname="${f%.mat*}.mat"
        
        echo "Rinomino: '$f' -> '$newname'"
        mv -- "$f" "$newname"
    fi
done