# Output file
set outfile [open aromatic_distance.dat w]

# Write header
puts $outfile "Frame\tTYR89_C6_Dist\tTYR89_C6_Int\tTYR89_C5_Dist\tTYR89_C5_Int\tTYR92_C6_Dist\tTYR92_C6_Int\tTYR92_C5_Dist\tTYR92_C5_Int"

# Total frames
set nf [molinfo top get numframes]

# =========================
# Counters for occupancy
# =========================

set count_89_c6 0
set count_89_c5 0

set count_92_c6 0
set count_92_c5 0

# =========================
# Sum of interacting distances
# =========================

set sum_89_c6 0
set sum_89_c5 0

set sum_92_c6 0
set sum_92_c5 0

# =========================
# Initialize min/max distances
# =========================

set min_89_c6 9999
set max_89_c6 0

set min_89_c5 9999
set max_89_c5 0

set min_92_c6 9999
set max_92_c6 0

set min_92_c5 9999
set max_92_c5 0

# =========================
# Loop over trajectory
# =========================

for {set i 0} {$i < $nf} {incr i} {

    animate goto $i

    # TYR89 aromatic ring
    set tyr89_ring [atomselect top "resid 89 and resname TYR and name CE1 CD1 CG CD2 CE2 CZ"]

    # TYR92 aromatic ring
    set tyr92_ring [atomselect top "resid 92 and resname TYR and name CE1 CD1 CG CD2 CE2 CZ"]

    # Ligand C6 ring
    set c6ring [atomselect top "resname RXT and name NAP CAT CAS CAR NAM CAE"]

    # Ligand C5 ring
    set c5ring [atomselect top "resname RXT and name CAT CAS NAN CAD CAC"]

    # =========================
    # Calculate centroids
    # =========================

    set tyr89_cent [measure center $tyr89_ring weight mass]
    set tyr92_cent [measure center $tyr92_ring weight mass]

    set c6_cent [measure center $c6ring weight mass]
    set c5_cent [measure center $c5ring weight mass]

    # =========================
    # Distances
    # =========================

    set dist_89_c6 [veclength [vecsub $tyr89_cent $c6_cent]]
    set dist_89_c5 [veclength [vecsub $tyr89_cent $c5_cent]]

    set dist_92_c6 [veclength [vecsub $tyr92_cent $c6_cent]]
    set dist_92_c5 [veclength [vecsub $tyr92_cent $c5_cent]]

    # =========================
    # TYR89 - C6
    # =========================

    if {$dist_89_c6 <= 7.0} {

        set int_89_c6 "YES"

        incr count_89_c6

        set sum_89_c6 [expr {$sum_89_c6 + $dist_89_c6}]

        if {$dist_89_c6 < $min_89_c6} {
            set min_89_c6 $dist_89_c6
        }

        if {$dist_89_c6 > $max_89_c6} {
            set max_89_c6 $dist_89_c6
        }

    } else {

        set int_89_c6 "NO"
    }

    # =========================
    # TYR89 - C5
    # =========================

    if {$dist_89_c5 <= 7.0} {

        set int_89_c5 "YES"

        incr count_89_c5

        set sum_89_c5 [expr {$sum_89_c5 + $dist_89_c5}]

        if {$dist_89_c5 < $min_89_c5} {
            set min_89_c5 $dist_89_c5
        }

        if {$dist_89_c5 > $max_89_c5} {
            set max_89_c5 $dist_89_c5
        }

    } else {

        set int_89_c5 "NO"
    }

    # =========================
    # TYR92 - C6
    # =========================

    if {$dist_92_c6 <= 7.0} {

        set int_92_c6 "YES"

        incr count_92_c6

        set sum_92_c6 [expr {$sum_92_c6 + $dist_92_c6}]

        if {$dist_92_c6 < $min_92_c6} {
            set min_92_c6 $dist_92_c6
        }

        if {$dist_92_c6 > $max_92_c6} {
            set max_92_c6 $dist_92_c6
        }

    } else {

        set int_92_c6 "NO"
    }

    # =========================
    # TYR92 - C5
    # =========================

    if {$dist_92_c5 <= 7.0} {

        set int_92_c5 "YES"

        incr count_92_c5

        set sum_92_c5 [expr {$sum_92_c5 + $dist_92_c5}]

        if {$dist_92_c5 < $min_92_c5} {
            set min_92_c5 $dist_92_c5
        }

        if {$dist_92_c5 > $max_92_c5} {
            set max_92_c5 $dist_92_c5
        }

    } else {

        set int_92_c5 "NO"
    }

    # =========================
    # Write output
    # =========================

    puts $outfile "$i\t$dist_89_c6\t$int_89_c6\t$dist_89_c5\t$int_89_c5\t$dist_92_c6\t$int_92_c6\t$dist_92_c5\t$int_92_c5"

    # Delete selections
    $tyr89_ring delete
    $tyr92_ring delete
    $c6ring delete
    $c5ring delete
}

# =========================
# Occupancy calculations
# =========================

set occ_89_c6 [expr {($count_89_c6*100.0)/$nf}]
set occ_89_c5 [expr {($count_89_c5*100.0)/$nf}]

set occ_92_c6 [expr {($count_92_c6*100.0)/$nf}]
set occ_92_c5 [expr {($count_92_c5*100.0)/$nf}]

# =========================
# Average distances
# =========================

if {$count_89_c6 > 0} {
    set avg_89_c6 [expr {$sum_89_c6 / $count_89_c6}]
} else {
    set avg_89_c6 0
}

if {$count_89_c5 > 0} {
    set avg_89_c5 [expr {$sum_89_c5 / $count_89_c5}]
} else {
    set avg_89_c5 0
}

if {$count_92_c6 > 0} {
    set avg_92_c6 [expr {$sum_92_c6 / $count_92_c6}]
} else {
    set avg_92_c6 0
}

if {$count_92_c5 > 0} {
    set avg_92_c5 [expr {$sum_92_c5 / $count_92_c5}]
} else {
    set avg_92_c5 0
}

# =========================
# Write results
# =========================

puts $outfile "\n===== Occupancy ====="

puts $outfile "TYR89-C6ring Occupancy = $occ_89_c6 %"
puts $outfile "TYR89-C5ring Occupancy = $occ_89_c5 %"

puts $outfile "TYR92-C6ring Occupancy = $occ_92_c6 %"
puts $outfile "TYR92-C5ring Occupancy = $occ_92_c5 %"

puts $outfile "\n===== Average Interacting Distances ====="

puts $outfile "Average TYR89-C6ring Distance = $avg_89_c6 Å"
puts $outfile "Average TYR89-C5ring Distance = $avg_89_c5 Å"

puts $outfile "Average TYR92-C6ring Distance = $avg_92_c6 Å"
puts $outfile "Average TYR92-C5ring Distance = $avg_92_c5 Å"

puts $outfile "\n===== Minimum and Maximum Interacting Distances ====="

puts $outfile "Minimum TYR89-C6ring Distance = $min_89_c6 Å"
puts $outfile "Maximum TYR89-C6ring Distance = $max_89_c6 Å"

puts $outfile "Minimum TYR89-C5ring Distance = $min_89_c5 Å"
puts $outfile "Maximum TYR89-C5ring Distance = $max_89_c5 Å"

puts $outfile "Minimum TYR92-C6ring Distance = $min_92_c6 Å"
puts $outfile "Maximum TYR92-C6ring Distance = $max_92_c6 Å"

puts $outfile "Minimum TYR92-C5ring Distance = $min_92_c5 Å"
puts $outfile "Maximum TYR92-C5ring Distance = $max_92_c5 Å"

close $outfile

puts "Finished centroid distance calculation"
