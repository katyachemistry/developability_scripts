#!/usr/bin/awk -f
#
# Faster, functionally identical rewrite of parse_sasa.awk
#
# Usage:
#   awk -v pdb=structure.pdb -f parse_sasa_fast.awk sasa_file.sasa
#
# Notes on optimizations (kept behavior identical):
# - Avoid repeated full-array scans by building per-residue atom lists once.
# - Use residue->atom-index mappings (res_t_index, res_p_index, res_n_index).
# - When computing pairwise residue distances, iterate only over atoms that belong
#   to the two residues being compared (no global scans filttoms.
# - Preserve cutoffs, SASA thresholds, and output formatting exactly.ering by equality).
# - Avoid rebuilding mappings inside loops; update maps only when adding a
#


function aa3to1(aa3,  map) {
    # standard 3-to-1 mapping
    map["ALA"]="A"; map["ARG"]="R"; map["ASN"]="N"; map["ASP"]="D";
    map["CYS"]="C"; map["GLN"]="Q"; map["GLU"]="E"; map["GLY"]="G";
    map["HIS"]="H"; map["ILE"]="I"; map["LEU"]="L"; map["LYS"]="K";
    map["MET"]="M"; map["PHE"]="F"; map["PRO"]="P"; map["SER"]="S";
    map["THR"]="T"; map["TRP"]="W"; map["TYR"]="Y"; map["VAL"]="V";

    aa3 = toupper(aa3)
    return (aa3 in map ? map[aa3] : "X")
}

BEGIN {
    # CDR ranges for IMGT-numbered structures
    hcdr = 0; lcdr = 0

    # Heavy chain IMGT CDRs
    hcdr++; cdrHstart[hcdr] = 26; cdrHend[hcdr] = 39   # CDR1
    hcdr++; cdrHstart[hcdr] = 55; cdrHend[hcdr] = 66   # CDR2
    hcdr++; cdrHstart[hcdr] = 104; cdrHend[hcdr] = 118 # CDR3

    # Light chain IMGT CDRs
    lcdr++; cdrLstart[lcdr] = 26; cdrLend[lcdr] = 39   # CDR1
    lcdr++; cdrLstart[lcdr] = 55; cdrLend[lcdr] = 66   # CDR2
    lcdr++; cdrLstart[lcdr] = 104; cdrLend[lcdr] = 118 # CDR3

    # charges and KD scale (same as original)
    charge["LYS"]=1; charge["ARG"]=1; charge["HIS"]=0.1
    charge["ASP"]=-1; charge["GLU"]=-1

    kd["ILE"]=4.5; kd["VAL"]=4.2; kd["LEU"]=3.8; kd["PHE"]=2.8; kd["CYS"]=2.5;
    kd["MET"]=1.9; kd["ALA"]=1.8; kd["GLY"]=-0.4; kd["THR"]=-0.7; kd["SER"]=-0.8;
    kd["TRP"]=-0.9; kd["TYR"]=-1.3; kd["PRO"]=-1.6; kd["HIS"]=-3.2;
    kd["GLU"]=-3.5; kd["GLN"]=-3.5; kd["ASP"]=-3.5; kd["ASN"]=-3.5;
    kd["LYS"]=-3.9; kd["ARG"]=-4.5
    minKD = -4.5; maxKD = 4.5; range = maxKD - minKD

    # --- PASS 1: read PDB and build per-residue atom lists
    pcount = 0; ncount = 0; tcount = 0

    # collect CDR sequences directly while reading PDB
    HCDR1_seq = HCDR2_seq = HCDR3_seq = ""
    LCDR1_seq = LCDR2_seq = LCDR3_seq = ""
    L_seq = ""
    H_seq = ""

    # track last residue added to avoid duplicate atoms
    delete seen_residue_cdr

    while ((getline line < pdb) > 0) {
        if (substr(line,1,4) != "ATOM") continue

        # --- Residue info ---
        res   = substr(line,18,3)
        gsub(/^ +| +$/, "", res)

        chain = substr(line,22,1)

        resnum = substr(line,23,4)   # numeric part (e.g., 111, 112)
        gsub(/^ +| +$/, "", resnum)

        icode = substr(line,27,1)    # insertion code (can be space or letter)
        gsub(/^ +| +$/, "", icode)

        # full residue ID including insertion code, used for keys
        num = resnum icode
        reskey = chain ":" num

        atom = substr(line,13,4)
        x = (substr(line,31,8) + 0)
        y = (substr(line,39,8) + 0)
        z = (substr(line,47,8) + 0)

        # --- numeric part for CDR range comparisons ---
        if (resnum ~ /^[0-9]+/) num_digits = resnum + 0
        else num_digits = 0   # fallback, should not happen

        # debug printing
#        print "  Current residue:", res, "num:", resnum, "icode:", icode, "full num:", num > "/dev/stderr"

        # --- CDR SEQUENCE COLLECTION (ONLY ON FIRST ATOM OF RESIDUE) ---
        if (!(reskey in seen_residue_cdr)) {

        # heavy chain
        if (chain == "H") {

            if (num_digits >= 27 && num_digits <= 38) {
                HCDR1_seq = HCDR1_seq aa3to1(res)
#                print "  Added to HCDR1_seq: " aa3to1(res) > "/dev/stderr"
            }
            else if (num_digits >= 56 && num_digits <= 65) {
                HCDR2_seq = HCDR2_seq aa3to1(res)
#                print "  Added to HCDR2_seq: " aa3to1(res) > "/dev/stderr"
            }
            else if (num_digits >= 105 && num_digits <= 117) {
                HCDR3_seq = HCDR3_seq aa3to1(res)
#                print "  Added to HCDR3_seq: " aa3to1(res) > "/dev/stderr"
            }

            H_seq = H_seq aa3to1(res)
        }

        # light chain
        else if (chain == "L") {

            if (num_digits >= 27 && num_digits <= 38) {
                LCDR1_seq = LCDR1_seq aa3to1(res)
#                print "  Added to LCDR1_seq: " aa3to1(res) > "/dev/stderr"
            }
            else if (num_digits >= 56 && num_digits <= 65) {
                LCDR2_seq = LCDR2_seq aa3to1(res)
#                print "  Added to LCDR2_seq: " aa3to1(res) > "/dev/stderr"
            }
            else if (num_digits >= 105 && num_digits <= 117) {
                LCDR3_seq = LCDR3_seq aa3to1(res)
#                print "  Added to LCDR3_seq: " aa3to1(res) > "/dev/stderr"
            }

            L_seq = L_seq aa3to1(res)
        }

            # mark this residue as processed
            seen_residue_cdr[reskey] = 1
        }
        

        # store heavy atoms (exclude hydrogens)
        if (atom !~ /^ H/) {
            tcount++
            tx[tcount]=x; ty[tcount]=y; tz[tcount]=z
            tres[tcount]=res
            tkey[tcount]=reskey

            # append this heavy-atom index to residue->heavy-atom list
            res_t_count[reskey]++
            res_t_index[reskey, res_t_count[reskey]] = tcount
        }

        # charged positive side-chain atoms (for salt-bridges and PPC)
        if ((res=="LYS" && atom ~ /NZ/) || (res=="ARG" && atom ~ /NH[12]/) || (res=="HIS" && atom ~ /ND1/)) {
            pcount++
            px[pcount]=x; py[pcount]=y; pz[pcount]=z
            p_res[pcount]=res
            pkey[pcount]=reskey

            # append p-atom index to residue's positive-atom list
            res_p_count[reskey]++
            res_p_index[reskey, res_p_count[reskey]] = pcount

            # mark residue as potentially positive charged
            allpos[reskey] = 1
            # set atomic residue-level charge (may be set to 0 later if salt-bridge)
            charge_atom[reskey] = charge[res]
        }

        # charged negative side-chain atoms
        if ((res=="ASP" && atom ~ /OD[12]/) || (res=="GLU" && atom ~ /OE[12]/)) {
            ncount++
            nx[ncount]=x; ny[ncount]=y; nz[ncount]=z
            n_res[ncount]=res
            nkey[ncount]=reskey

            # append n-atom index to residue's negative-atom list
            res_n_count[reskey]++
            res_n_index[reskey, res_n_count[reskey]] = ncount

            allneg[reskey] = 1
            charge_atom[reskey] = charge[res]
        }
    }
    close(pdb)

    # --- identify CDR vicinity heavy atoms (4 Å radius around CDR/anchor residues)
    cdrvicinity_cutoff2 = 4.0 * 4.0

    # collect all heavy atoms that belong to CDR residues (and anchors)
    # for IMGT-numbered residues, we use string keys instead of numeric loops

    # heavy chain
    for (r = 1; r <= hcdr; r++) {
        for (i = 1; i <= tcount; i++) {
            if (tres[i] != "") {
                reskey_i = tkey[i]
                num_i = substr(reskey_i, 3)      # remove chain letter + colon
                # extract numeric prefix for comparison
                match(num_i, /^[0-9]+/)
                num_i_digits = substr(num_i, RSTART, RLENGTH) + 0
                start_digits = cdrHstart[r]
                end_digits   = cdrHend[r]
                sub(/[^0-9].*$/, "", num_i_digits)   # e.g., 111A -> 111
                if (num_i_digits >= start_digits && num_i_digits <= end_digits) {
                    cdr_res[reskey_i] = 1
                }
            }
        }
    }


    # Convert sets to counts
    HCDR1_len = length(HCDR1_seq)
    HCDR2_len = length(HCDR2_seq)
    HCDR3_len = length(HCDR3_seq)
    LCDR1_len = length(LCDR1_seq)
    LCDR2_len = length(LCDR2_seq)
    LCDR3_len = length(LCDR3_seq)

    # light chain
    for (r = 1; r <= lcdr; r++) {
        for (i = 1; i <= tcount; i++) {
            if (tres[i] != "") {
                reskey_i = tkey[i]
                num_i = substr(reskey_i, 3)
                # extract numeric prefix
                match(num_i, /^[0-9]+/)
                num_i_digits = substr(num_i, RSTART, RLENGTH) + 0
                start_digits = cdrLstart[r]
                end_digits   = cdrLend[r]
                sub(/[^0-9].*$/, "", num_i_digits)
                if (num_i_digits >= start_digits && num_i_digits <= end_digits) {
                    cdr_res[reskey_i] = 1
                }
            }
        }
    }
    
    # detect salt bridges (3.2 Å cutoff -> 10.24 squared)
    sb_cutoff2 = 3.2 * 3.2
    for (i=1; i<=pcount; i++) {
        for (j=1; j<=ncount; j++) {
            dx=px[i]-nx[j]; dy=py[i]-ny[j]; dz=pz[i]-nz[j]
            if (dx*dx+dy*dy+dz*dz <= sb_cutoff2) {
                sbpos[pkey[i]] = 1
                sbneg[nkey[j]] = 1
                charge_atom[pkey[i]] = 0
                charge_atom[nkey[j]] = 0
            }
        }
    }

    # precompute normalized hydrophobicities
    for (r in kd) norm[r] = 1 + (kd[r] - minKD)/range  # maps to 1..2

    # cutoff for patch distance squared (7.5 Å)
    cutoff2 = 7.5 * 7.5

    # counters for exposed atoms/residues
    te_count = 0
    tc_cdr_count_exp = 0
    hydSASA = 0; aromSASA = 0
}

# --- PASS 2: parse SASA file (sasa_file.sasa)
# Keep original filters for lines: ignore REM, END, CHAIN, TOTAL
!/^REM/ && !/^END/ && !/^CHAIN/ && !/^TOTAL/ {
    # $8 is relative SASA (we only consider numeric values and >=5.00)
    rel = ($8 + 0)
    if (rel != $8) next
    if (rel < 5.00) next

    res = $2
    chain = $3
    num = $4
    gsub(/^ +| +$/, "", num)
    out = chain ":" num

    # charged exposed residues (skip if salt-bridged)
    if ((res=="LYS" || res=="ARG" || res=="HIS") && !(out in sbpos)) {
        exppos[out] = res
    }
    if ((res=="ASP" || res=="GLU") && !(out in sbneg)) {
        expneg[out] = res
    }

    # hydrophobic/aromatic SASA accumulation uses column $7 (absolute ASA)
    absASA = ($7 + 0)
    if (res=="PHE" || res=="TYR" || res=="TRP" || res=="HIS") aromSASA += absASA
    if (res=="ALA" || res=="VAL" || res=="LEU" || res=="ILE" ||
        res=="MET" || res=="PRO" || res=="PHE" || res=="TRP") hydSASA += absASA

    # if residue type has hydrophobicity entry, treat as exposed for PSH computations
    if (res in kd) {
        exposed[out] = res

        # populate exposed heavy-atom index list for this residue using res_t_index built earlier
        if (res_t_count[out] > 0) {
            exp_res_count[out] = 0
            for (k = 1; k <= res_t_count[out]; k++) {
                tidx = res_t_index[out, k]
                exp_res_count[out]++
                exp_res_index[out, exp_res_count[out]] = tidx
            }
        } else {
            # no heavy atoms found in PDB for this residue: leave count zero
            exp_res_count[out] = 0
        }
    }

    # precompute which residues are in CDR vicinity by heavy-atom proximity
    for (i = 1; i <= tcount; i++) {
        key_i = tkey[i]
        
        # only consider surface-exposed residues for vicinity
#        if (!(key_i in exposed)) continue

        for (j = 1; j <= tcount; j++) {
            key_j = tkey[j]
            
            # only compare to residues in CDRs (cdr_res has string keys)
            if (!(key_j in cdr_res)) continue

            # compute squared distance between heavy atoms
            dx = tx[i] - tx[j]
            dy = ty[i] - ty[j]
            dz = tz[i] - tz[j]
            d2 = dx*dx + dy*dy + dz*dz

            if (d2 <= cdrvicinity_cutoff2) {
                cdr_vicinity[key_i] = 1
                break  # no need to check other CDR atoms once a match is found
            }
        }
    }

    # check if residue is in CDR vicinity (CDR + 4 Å neighborhood)
    inCDR_vicinity = 0
    if ((chain ":" num) in cdr_res || (chain ":" num) in cdr_vicinity) {
        inCDR_vicinity = 1
    }

    if (inCDR_vicinity) {
        if ((res=="LYS" || res=="ARG" || res=="HIS") && !(out in sbpos)) exppos_cdr[out] = res
        if ((res=="ASP" || res=="GLU") && !(out in sbneg)) expneg_cdr[out] = res
        if (res in kd) exposed_cdr[out] = res

        # store exposed CDR-vicinity heavy atoms
        if (res in kd) {
            if (res_t_count[out] > 0) {
                exp_cdr_res_count[out] = 0
                for (k = 1; k <= res_t_count[out]; k++) {
                    tidx = res_t_index[out, k]
                    exp_cdr_res_count[out]++
                    exp_cdr_res_index[out, exp_cdr_res_count[out]] = tidx
                }
            } else {
                exp_cdr_res_count[out] = 0
            }
        }
    }

}

END {
    # --- compute nchg (number of charged residues INCLUDING buried and salt-bridged)
    nchg = 0
    for (k in allpos) nchg++
    for (k in allneg) nchg++

    # --- Charge asymmetry (VH - VL)
    VH_charge = 0; VL_charge = 0
    for (k in exppos) {
        if (k in sbpos) continue
        ch = charge[exppos[k]]
        split(k, tmp, ":"); chain = tmp[1]
        if (chain == "H") VH_charge += ch
        else if (chain == "L") VL_charge += ch
    }
    for (k in expneg) {
        if (k in sbneg) continue
        ch = charge[expneg[k]]
        split(k, tmp, ":"); chain = tmp[1]
        if (chain == "H") VH_charge += ch
        else if (chain == "L") VL_charge += ch
    }
    charge_asymmetry = VH_charge * VL_charge

    # --- Hydrophobic Patch Index (PSHmetric) over all exposed residues (non-CDR)
    ec = 0
    for (k in exposed) {
        if (exp_res_count[k] > 0) {
            ec++
            ekey[ec] = k
            eres[ec] = exposed[k]
        }
    }

    psh = 0
    for (i = 1; i < ec; i++) {
        for (j = i+1; j <= ec; j++) {
            r1 = ekey[i]; r2 = ekey[j]
            mind2 = 1e18
            for (a = 1; a <= exp_res_count[r1]; a++) {
                ai = exp_res_index[r1, a]
                ax = tx[ai]; ay = ty[ai]; az = tz[ai]
                for (b = 1; b <= exp_res_count[r2]; b++) {
                    bi = exp_res_index[r2, b]
                    dx = ax - tx[bi]; dy = ay - ty[bi]; dz = az - tz[bi]
                    d2 = dx*dx + dy*dy + dz*dz
                    if (d2 < mind2) mind2 = d2
                }
            }
            if (mind2 <= cutoff2) {
                # if salt-bridged, use Gly KD instead of native residue KD
                if ((r1 in sbpos) || (r1 in sbneg)) h1 = norm["GLY"]; else h1 = norm[eres[i]]
                if ((r2 in sbpos) || (r2 in sbneg)) h2 = norm["GLY"]; else h2 = norm[eres[j]]
                psh += h1 * h2 / mind2
            }
        }
    }

    # --- Hydrophobic Patch Index for CDR residues
    ec = 0
    for (k in exposed_cdr) {
        if (exp_cdr_res_count[k] > 0) {
            ec++
            ekey[ec] = k
            eres[ec] = exposed_cdr[k]
        }
    }

    psh_cdr = 0
    for (i = 1; i < ec; i++) {
        for (j = i+1; j <= ec; j++) {
            r1 = ekey[i]; r2 = ekey[j]
            mind2 = 1e18
            for (a = 1; a <= exp_cdr_res_count[r1]; a++) {
                ai = exp_cdr_res_index[r1, a]
                ax = tx[ai]; ay = ty[ai]; az = tz[ai]
                for (b = 1; b <= exp_cdr_res_count[r2]; b++) {
                    bi = exp_cdr_res_index[r2, b]
                    dx = ax - tx[bi]; dy = ay - ty[bi]; dz = az - tz[bi]
                    d2 = dx*dx + dy*dy + dz*dz
                    if (d2 < mind2) mind2 = d2
                }
            }
            if (mind2 <= cutoff2) {
                if ((r1 in sbpos) || (r1 in sbneg)) h1 = norm["GLY"]; else h1 = norm[eres[i]]
                if ((r2 in sbpos) || (r2 in sbneg)) h2 = norm["GLY"]; else h2 = norm[eres[j]]
                psh_cdr += h1 * h2 / mind2
            }
        }
    }

    # --- Positive Patch Compactness (PPC): GLOBAL ---
    # Build list of exposed, positive, NON–salt-bridged residues
    pc = 0
    for (k in exppos) {
        if (k in sbpos) continue
        if (exp_res_count[k] == 0) continue
        pc++
        pkey2[pc] = k
    }

    ppc = 0
    for (i = 1; i < pc; i++) {
        for (j = i+1; j <= pc; j++) {
            r1 = pkey2[i]
            r2 = pkey2[j]

            mind2 = 1e18
            for (a = 1; a <= exp_res_count[r1]; a++) {
                ai = exp_res_index[r1, a]
                ax = tx[ai]; ay = ty[ai]; az = tz[ai]
                for (b = 1; b <= exp_res_count[r2]; b++) {
                    bi = exp_res_index[r2, b]
                    dx = ax - tx[bi]
                    dy = ay - ty[bi]
                    dz = az - tz[bi]
                    d2 = dx*dx + dy*dy + dz*dz
                    if (d2 < mind2) mind2 = d2
                }
            }

            if (mind2 <= cutoff2) {
                # both residues are positive and non–salt-bridged
                q1 = charge_atom[r1]
                q2 = charge_atom[r2]
                ppc += q1 * q2 / mind2
            }
        }
    }

    # --- Positive Patch Compactness (PPC): CDR VICINITY ---
    pc = 0
    for (k in exppos_cdr) {
        if (k in sbpos) continue
        if (exp_cdr_res_count[k] == 0) continue
        pc++
        pkey2[pc] = k
    }

    ppc_cdr = 0
    for (i = 1; i < pc; i++) {
        for (j = i+1; j <= pc; j++) {
            r1 = pkey2[i]
            r2 = pkey2[j]

            mind2 = 1e18
            for (a = 1; a <= exp_cdr_res_count[r1]; a++) {
                ai = exp_cdr_res_index[r1, a]
                ax = tx[ai]; ay = ty[ai]; az = tz[ai]
                for (b = 1; b <= exp_cdr_res_count[r2]; b++) {
                    bi = exp_cdr_res_index[r2, b]
                    dx = ax - tx[bi]
                    dy = ay - ty[bi]
                    dz = az - tz[bi]
                    d2 = dx*dx + dy*dy + dz*dz
                    if (d2 < mind2) mind2 = d2
                }
            }

            if (mind2 <= cutoff2) {
                q1 = charge_atom[r1]
                q2 = charge_atom[r2]
                ppc_cdr += q1 * q2 / mind2
            }
        }
    }


    # --- Negative Patch Compactness (PNC) global and CDR ---
    nc = 0
    for (k in expneg) if (!(k in sbneg)) {
        nc++
        nkey2[nc] = k
        nres[nc]  = expneg[k]
    }

    pnc = 0
    for (i = 1; i < nc; i++) {
        for (j = i+1; j <= nc; j++) {
            r1 = nkey2[i]; r2 = nkey2[j]
            mind2 = 1e18
            for (a = 1; a <= exp_res_count[r1]; a++) {
                ai = exp_res_index[r1, a]
                ax = tx[ai]; ay = ty[ai]; az = tz[ai]
                for (b = 1; b <= exp_res_count[r2]; b++) {
                    bi = exp_res_index[r2, b]
                    dx = ax - tx[bi]; dy = ay - ty[bi]; dz = az - tz[bi]
                    d2 = dx*dx + dy*dy + dz*dz
                    if (d2 < mind2) mind2 = d2
                }
            }
            if (mind2 <= cutoff2) {
                q1 = charge_atom[r1]; q2 = charge_atom[r2]
                pnc += q1 * q2 / mind2
            }
        }
    }

    # CDR PNC
    nc = 0
    for (k in expneg_cdr) if (!(k in sbneg)) {
        nc++
        nkey2[nc] = k
        nres[nc]  = expneg_cdr[k]
    }

    pnc_cdr = 0
    for (i = 1; i < nc; i++) {
        for (j = i+1; j <= nc; j++) {
            r1 = nkey2[i]; r2 = nkey2[j]
            mind2 = 1e18
            for (a = 1; a <= exp_cdr_res_count[r1]; a++) {
                ai = exp_cdr_res_index[r1, a]
                ax = tx[ai]; ay = ty[ai]; az = tz[ai]
                for (b = 1; b <= exp_cdr_res_count[r2]; b++) {
                    bi = exp_cdr_res_index[r2, b]
                    dx = ax - tx[bi]; dy = ay - ty[bi]; dz = az - tz[bi]
                    d2 = dx*dx + dy*dy + dz*dz
                    if (d2 < mind2) mind2 = d2
                }
            }
            if (mind2 <= cutoff2) {
                q1 = charge_atom[r1]; q2 = charge_atom[r2]
                pnc_cdr += q1 * q2 / mind2
            }
        }
    }

    # --- Debug: count exposed CDR residues
    exposed_cdr_count = 0
    for (k in exposed_cdr) exposed_cdr_count++
    exposed_count = 0
    for (k in exposed) exposed_count++

    printf "%s,%.4f,%d,%.2f,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%s,%d,%s,%d,%s,%d,%s,%d,%s,%d,%s,%d,%s,%s\n", \
        filename, \
        charge_asymmetry, \
        nchg, \
        hydSASA, \
        aromSASA, \
        psh, ppc, pnc, psh_cdr, ppc_cdr, pnc_cdr, \
        HCDR1_seq, HCDR1_len, \
        LCDR1_seq, LCDR1_len, \
        HCDR2_seq, HCDR2_len, \
        LCDR2_seq, LCDR2_len, \
        HCDR3_seq, HCDR3_len, \
        LCDR3_seq, LCDR3_len, \
        H_seq, L_seq

}
