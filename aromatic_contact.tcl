# ============================================================
# Aromatic Ring Centroid Distance Analysis
# ============================================================
#
# Software: VMD
# Language: Tcl
#
# Description:
# Calculates mass-weighted centroid-to-centroid distances
# between a protein TYR aromatic ring and two ligand rings
# (C6 and C5) throughout an MD trajectory.
#
# Method:
# 1. Calculate the mass-weighted center of each aromatic ring.
# 2. Calculate the Euclidean distance between the TYR ring
#    centroid and each ligand ring centroid.
# 3. Define an interaction when the centroid distance is
#    <= 7.0 Angstrom.
# 4. Calculate interaction occupancy over the trajectory.
# 5. Calculate average, minimum, and maximum distances using
#    only interacting frames.
#
# ============================================================


# ============================================================
# USER-DEFINED PARAMETERS
# Edit ONLY this section for a different system.
# ============================================================

# Protein aromatic residue
set tyr_resid 1
set tyr_resname TYR

# Protein TYR aromatic ring atoms
set tyr_ring_atoms "CG CD1 CE1 CZ CE2 CD2"

# Ligand residue name
set ligand_resname LIG

# Ligand C6 ring atoms
set c6_ring_atoms "N1 C N C3 C2 C1"

# Ligand C5 ring atoms
set c5_ring_atoms "C1 C2 C12 C13 N4"

# Interaction cutoff in Angstrom
set cutoff 7.0

# Output file
set output_file "aromatic_distance.dat"


# ============================================================
# OUTPUT FILE
# ============================================================

set outfile [open $output_file w]

# Write header
puts $outfile "Frame\tC6ring_Dist\tC6ring_Int\tC5ring_Dist\tC5ring_Int"


# ============================================================
# TRAJECTORY INFORMATION
# ============================================================

set nf [molinfo top get numframes]


# ============================================================
# INITIALIZE COUNTERS
# ============================================================

# Interaction occupancy counters
set count_c6 0
set count_c5 0

# Sum of interacting distances
set sum_c6 0
set sum_c5 0

# Initialize minimum and maximum distances
set min_c6 9999
set max_c6 0

set min_c5 9999
set max_c5 0


# ============================================================
# TRAJECTORY LOOP
# ============================================================

for {set i 0} {$i < $nf} {incr i} {

    animate goto $i


    # --------------------------------------------------------
    # Atom selections
    # --------------------------------------------------------

    # Protein TYR aromatic ring
    set tyr_ring [atomselect top \
        "resid $tyr_resid and resname $tyr_resname and name $tyr_ring_atoms"]

    # Ligand C6 ring
    set c6ring [atomselect top \
        "resname $ligand_resname and name $c6_ring_atoms"]

    # Ligand C5 ring
    set c5ring [atomselect top \
        "resname $ligand_resname and name $c5_ring_atoms"]


    # --------------------------------------------------------
    # Calculate mass-weighted ring centers
    # --------------------------------------------------------

    set tyr_cent [measure center $tyr_ring weight mass]
    set c6_cent  [measure center $c6ring weight mass]
    set c5_cent  [measure center $c5ring weight mass]


    # --------------------------------------------------------
    # Calculate centroid-to-centroid distances
    # --------------------------------------------------------

    set dist_c6 [veclength [vecsub $tyr_cent $c6_cent]]
    set dist_c5 [veclength [vecsub $tyr_cent $c5_cent]]


    # ========================================================
    # C6 RING INTERACTION
    # ========================================================

    if {$dist_c6 <= $cutoff} {

        set interact_c6 "YES"

        # Count interacting frames
        incr count_c6

        # Sum interacting distances
        set sum_c6 [expr {$sum_c6 + $dist_c6}]

        # Minimum interacting distance
        if {$dist_c6 < $min_c6} {
            set min_c6 $dist_c6
        }

        # Maximum interacting distance
        if {$dist_c6 > $max_c6} {
            set max_c6 $dist_c6
        }

    } else {

        set interact_c6 "NO"
    }


    # ========================================================
    # C5 RING INTERACTION
    # ========================================================

    if {$dist_c5 <= $cutoff} {

        set interact_c5 "YES"

        # Count interacting frames
        incr count_c5

        # Sum interacting distances
        set sum_c5 [expr {$sum_c5 + $dist_c5}]

        # Minimum interacting distance
        if {$dist_c5 < $min_c5} {
            set min_c5 $dist_c5
        }

        # Maximum interacting distance
        if {$dist_c5 > $max_c5} {
            set max_c5 $dist_c5
        }

    } else {

        set interact_c5 "NO"
    }


    # ========================================================
    # WRITE FRAME-WISE RESULTS
    # ========================================================

    puts $outfile "$i\t$dist_c6\t$interact_c6\t$dist_c5\t$interact_c5"


    # ========================================================
    # DELETE ATOM SELECTIONS
    # ========================================================

    $tyr_ring delete
    $c6ring delete
    $c5ring delete
}


# ============================================================
# CALCULATE OCCUPANCY
# ============================================================

set occ_c6 [expr {($count_c6 * 100.0) / $nf}]
set occ_c5 [expr {($count_c5 * 100.0) / $nf}]


# ============================================================
# CALCULATE AVERAGE INTERACTING DISTANCES
# ============================================================

if {$count_c6 > 0} {
    set avg_c6 [expr {$sum_c6 / $count_c6}]
} else {
    set avg_c6 0
}

if {$count_c5 > 0} {
    set avg_c5 [expr {$sum_c5 / $count_c5}]
} else {
    set avg_c5 0
}


# ============================================================
# WRITE SUMMARY TO OUTPUT FILE
# ============================================================

puts $outfile "\n===== Occupancy ====="

puts $outfile "C6ring Occupancy = $occ_c6 %"
puts $outfile "C5ring Occupancy = $occ_c5 %"


puts $outfile "\n===== Average Interacting Distances ====="

puts $outfile \
    "Average TYR$tyr_resid-C6ring Distance = $avg_c6 Angstrom"

puts $outfile \
    "Average TYR$tyr_resid-C5ring Distance = $avg_c5 Angstrom"


puts $outfile "\n===== Minimum and Maximum Interacting Distances ====="

puts $outfile \
    "Minimum TYR$tyr_resid-C6ring Distance = $min_c6 Angstrom"

puts $outfile \
    "Maximum TYR$tyr_resid-C6ring Distance = $max_c6 Angstrom"

puts $outfile \
    "Minimum TYR$tyr_resid-C5ring Distance = $min_c5 Angstrom"

puts $outfile \
    "Maximum TYR$tyr_resid-C5ring Distance = $max_c5 Angstrom"


# ============================================================
# CLOSE OUTPUT FILE
# ============================================================

close $outfile


# ============================================================
# PRINT SUMMARY TO VMD CONSOLE
# ============================================================

puts ""
puts "=============================================="
puts "Finished centroid distance calculation"
puts "=============================================="

puts "TYR residue       : $tyr_resname $tyr_resid"
puts "Ligand residue    : $ligand_resname"
puts "Interaction cutoff: $cutoff Angstrom"
puts "Total frames      : $nf"

puts ""
puts "C6ring Occupancy = $occ_c6 %"
puts "C5ring Occupancy = $occ_c5 %"

puts ""
puts "Average TYR$tyr_resid-C6ring Distance = $avg_c6 Angstrom"
puts "Average TYR$tyr_resid-C5ring Distance = $avg_c5 Angstrom"

puts ""
puts "Minimum TYR$tyr_resid-C6ring Distance = $min_c6 Angstrom"
puts "Maximum TYR$tyr_resid-C6ring Distance = $max_c6 Angstrom"

puts ""
puts "Minimum TYR$tyr_resid-C5ring Distance = $min_c5 Angstrom"
puts "Maximum TYR$tyr_resid-C5ring Distance = $max_c5 Angstrom"

puts ""
puts "Output file: $output_file"
puts "=============================================="
