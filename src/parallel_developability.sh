#!/bin/bash

# # Directories
STRUCTURE_DIR="/storage/antibody_data/PairedStructures/natural_models/OAS_models/OAS_structures_sample_10000/"
SASA_DIR="/storage/antibody_data/PairedStructures/natural_models/OAS_models/OAS_sasa"

# AWK script
AWK_SCRIPT="structural_properties_sasa_efficient.awk"
OUTPUT_CSV="natural_developability.csv"

# Directories
# STRUCTURE_DIR="structures_thera"
# SASA_DIR="sasa_thera"

# # AWK script
# AWK_SCRIPT="structural_properties_sasa_efficient.awk"

# OUTPUT_CSV="thera_developability_new.csv"

mkdir -p "$SASA_DIR"

# First, write header once
echo "Filename,ChargeAsymmetry,NCharged,HydrophobicSASA,AromaticSASA,PSH,PPC,PNC,PSH_CDR,PPC_CDR,PNC_CDR,HCDR1_seq,HCDR1_length,LCDR1_seq,LCDR1_length,HCDR2_seq,HCDR2_length,LCDR2_seq,LCDR2_length,HCDR3_seq,HCDR3_length,LCDR3_seq,LCDR3_length,H_seq,L_seq" > "$OUTPUT_CSV"

# Run jobs in parallel
find "$STRUCTURE_DIR" -name "*.pdb" | parallel -j 50 '
    pdbfile={}
    filename=$(basename "$pdbfile" .pdb)

    # Step 1: Run freesasa
    freesasa --shrake-rupley --format=rsa --depth=residue "$pdbfile" > "'"$SASA_DIR"'/$filename.sasa"

    # Step 2: Run AWK and append results
    awk -v filename="$filename" \
        -v pdb="$pdbfile" \
        -f "'"$AWK_SCRIPT"'" \
        "'"$SASA_DIR"'/$filename.sasa"
' >> "$OUTPUT_CSV"
