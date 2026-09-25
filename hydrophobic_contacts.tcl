# Hydrophobic Contact Occupancy Analysis
# Run: vmd -dispdev text -e hydrophobic.tcl

# User settings
set topology "top_file.prmtop"
set trajectory "traj_file.dcd"
set ligand_selection "resname LIG and name CAR CAS NAT C5 C6 C4 N3 C2 N1"
set protein_selection "protein and resname ALA VAL LEU ILE MET PHE TRP PRO"
set cutoff 4.5
set output_file "hydrophobic_occupancy.dat"

puts "Loading topology..."
mol new $topology

puts "Loading trajectory..."
mol addfile $trajectory waitfor all

puts "Trajectory loaded."

set nframes [molinfo top get numframes]
puts "Number of frames = $nframes"

# Check ligand selection
set lig_test [atomselect top $ligand_selection]

puts "Ligand atoms found: [$lig_test num]"

if {[$lig_test num] == 0} {
    puts "ERROR: No ligand atoms found."
    $lig_test delete
    exit
}

$lig_test delete

# Contact counter
array set counts {}

puts "Starting analysis..."

for {set f 0} {$f < $nframes} {incr f} {

    animate goto $f

    set lig [atomselect top $ligand_selection]
    set prot [atomselect top $protein_selection]

    set contacts [measure contacts $cutoff $lig $prot]

    set prot_idx [lsort -unique [lindex $contacts 1]]

    # Store residues contacted in this frame
    array unset frame_res
    array set frame_res {}

    foreach idx $prot_idx {

        set atm [atomselect top "index $idx"]

        set resid [lindex [$atm get resid] 0]
        set resname [lindex [$atm get resname] 0]
        set key "${resid}_${resname}"

        set frame_res($key) 1

        $atm delete
    }

    # Count each residue once per frame
    foreach key [array names frame_res] {

        if {[info exists counts($key)]} {
            incr counts($key)
        } else {
            set counts($key) 1
        }
    }

    $lig delete
    $prot delete

    if {$f % 1000 == 0} {
        puts "Frame $f / $nframes"
    }
}

puts "Writing results..."

set outfile [open $output_file w]

puts $outfile "Residue\tFrames\tOccupancy(%)"

foreach key [lsort [array names counts]] {

    set occ [expr {100.0 * $counts($key) / $nframes}]

    puts $outfile "$key\t$counts($key)\t[format %.2f $occ]"
}

close $outfile

puts "=================================="
puts "Done."
puts "Results written to:"
puts "$output_file"
puts "=================================="

exit
