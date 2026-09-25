# Aromatic Ring Centroid Distance and Occupancy Analysis
# Trajectory must already be loaded in VMD.
# Cutoff: 7.0 Angstrom
# Run: vmd -dispdev text -e aromatic_distance.tcl

# User settings
set output_file "aromatic_distance.dat"
set protein_resid_1 89
set protein_resname_1 "TYR"
set protein_resid_2 92
set protein_resname_2 "TYR"
set protein_ring_atoms "CE1 CD1 CG CD2 CE2 CZ"
set ligand_resname "LIG"
set ligand_ring_1_atoms "NAP CAT CAS CAR NAM CAE"
set ligand_ring_2_atoms "CAT CAS NAN CAD CAC"
set cutoff 7.0

set outfile [open $output_file w]

puts $outfile "Frame\tProtein1_Ligand1_Dist\tProtein1_Ligand1_Int\tProtein1_Ligand2_Dist\tProtein1_Ligand2_Int\tProtein2_Ligand1_Dist\tProtein2_Ligand1_Int\tProtein2_Ligand2_Dist\tProtein2_Ligand2_Int"

set nf [molinfo top get numframes]

# Initialize counters
set count_1_1 0
set count_1_2 0
set count_2_1 0
set count_2_2 0

set sum_1_1 0
set sum_1_2 0
set sum_2_1 0
set sum_2_2 0

set min_1_1 9999
set max_1_1 0
set min_1_2 9999
set max_1_2 0
set min_2_1 9999
set max_2_1 0
set min_2_2 9999
set max_2_2 0

# Analyze trajectory
for {set i 0} {$i < $nf} {incr i} {

    animate goto $i

    set protein_ring_1 [atomselect top "resid $protein_resid_1 and resname $protein_resname_1 and name $protein_ring_atoms"]
    set protein_ring_2 [atomselect top "resid $protein_resid_2 and resname $protein_resname_2 and name $protein_ring_atoms"]
    set ligand_ring_1 [atomselect top "resname $ligand_resname and name $ligand_ring_1_atoms"]
    set ligand_ring_2 [atomselect top "resname $ligand_resname and name $ligand_ring_2_atoms"]

    set protein_cent_1 [measure center $protein_ring_1 weight mass]
    set protein_cent_2 [measure center $protein_ring_2 weight mass]
    set ligand_cent_1 [measure center $ligand_ring_1 weight mass]
    set ligand_cent_2 [measure center $ligand_ring_2 weight mass]

    set dist_1_1 [veclength [vecsub $protein_cent_1 $ligand_cent_1]]
    set dist_1_2 [veclength [vecsub $protein_cent_1 $ligand_cent_2]]
    set dist_2_1 [veclength [vecsub $protein_cent_2 $ligand_cent_1]]
    set dist_2_2 [veclength [vecsub $protein_cent_2 $ligand_cent_2]]

    if {$dist_1_1 <= $cutoff} {
        set int_1_1 "YES"
        incr count_1_1
        set sum_1_1 [expr {$sum_1_1 + $dist_1_1}]
        if {$dist_1_1 < $min_1_1} {set min_1_1 $dist_1_1}
        if {$dist_1_1 > $max_1_1} {set max_1_1 $dist_1_1}
    } else {
        set int_1_1 "NO"
    }

    if {$dist_1_2 <= $cutoff} {
        set int_1_2 "YES"
        incr count_1_2
        set sum_1_2 [expr {$sum_1_2 + $dist_1_2}]
        if {$dist_1_2 < $min_1_2} {set min_1_2 $dist_1_2}
        if {$dist_1_2 > $max_1_2} {set max_1_2 $dist_1_2}
    } else {
        set int_1_2 "NO"
    }

    if {$dist_2_1 <= $cutoff} {
        set int_2_1 "YES"
        incr count_2_1
        set sum_2_1 [expr {$sum_2_1 + $dist_2_1}]
        if {$dist_2_1 < $min_2_1} {set min_2_1 $dist_2_1}
        if {$dist_2_1 > $max_2_1} {set max_2_1 $dist_2_1}
    } else {
        set int_2_1 "NO"
    }

    if {$dist_2_2 <= $cutoff} {
        set int_2_2 "YES"
        incr count_2_2
        set sum_2_2 [expr {$sum_2_2 + $dist_2_2}]
        if {$dist_2_2 < $min_2_2} {set min_2_2 $dist_2_2}
        if {$dist_2_2 > $max_2_2} {set max_2_2 $dist_2_2}
    } else {
        set int_2_2 "NO"
    }

    puts $outfile "$i\t$dist_1_1\t$int_1_1\t$dist_1_2\t$int_1_2\t$dist_2_1\t$int_2_1\t$dist_2_2\t$int_2_2"

    $protein_ring_1 delete
    $protein_ring_2 delete
    $ligand_ring_1 delete
    $ligand_ring_2 delete
}

# Occupancy
set occ_1_1 [expr {($count_1_1 * 100.0) / $nf}]
set occ_1_2 [expr {($count_1_2 * 100.0) / $nf}]
set occ_2_1 [expr {($count_2_1 * 100.0) / $nf}]
set occ_2_2 [expr {($count_2_2 * 100.0) / $nf}]

# Average distances
if {$count_1_1 > 0} {set avg_1_1 [expr {$sum_1_1 / $count_1_1}]} else {set avg_1_1 0}
if {$count_1_2 > 0} {set avg_1_2 [expr {$sum_1_2 / $count_1_2}]} else {set avg_1_2 0}
if {$count_2_1 > 0} {set avg_2_1 [expr {$sum_2_1 / $count_2_1}]} else {set avg_2_1 0}
if {$count_2_2 > 0} {set avg_2_2 [expr {$sum_2_2 / $count_2_2}]} else {set avg_2_2 0}

# Results
puts $outfile "\n===== Occupancy ====="
puts $outfile "Protein1-Ligand1 Occupancy = $occ_1_1 %"
puts $outfile "Protein1-Ligand2 Occupancy = $occ_1_2 %"
puts $outfile "Protein2-Ligand1 Occupancy = $occ_2_1 %"
puts $outfile "Protein2-Ligand2 Occupancy = $occ_2_2 %"

puts $outfile "\n===== Average Interacting Distances ====="
puts $outfile "Average Protein1-Ligand1 Distance = $avg_1_1 Angstrom"
puts $outfile "Average Protein1-Ligand2 Distance = $avg_1_2 Angstrom"
puts $outfile "Average Protein2-Ligand1 Distance = $avg_2_1 Angstrom"
puts $outfile "Average Protein2-Ligand2 Distance = $avg_2_2 Angstrom"

puts $outfile "\n===== Minimum and Maximum Interacting Distances ====="
puts $outfile "Minimum Protein1-Ligand1 Distance = $min_1_1 Angstrom"
puts $outfile "Maximum Protein1-Ligand1 Distance = $max_1_1 Angstrom"
puts $outfile "Minimum Protein1-Ligand2 Distance = $min_1_2 Angstrom"
puts $outfile "Maximum Protein1-Ligand2 Distance = $max_1_2 Angstrom"
puts $outfile "Minimum Protein2-Ligand1 Distance = $min_2_1 Angstrom"
puts $outfile "Maximum Protein2-Ligand1 Distance = $max_2_1 Angstrom"
puts $outfile "Minimum Protein2-Ligand2 Distance = $min_2_2 Angstrom"
puts $outfile "Maximum Protein2-Ligand2 Distance = $max_2_2 Angstrom"

close $outfile

puts "Finished aromatic centroid distance calculation."
puts "Results saved to: $output_file"
