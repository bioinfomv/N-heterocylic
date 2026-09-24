# ============================================================
# Hydrophobic Contact Occupancy Analysis
# ============================================================
#
# Purpose:
# Calculate the frame occupancy of hydrophobic protein residues
# contacting selected ligand atoms during an MD trajectory.
#
# Calculation method:
#   1. Select specified ligand atoms.
#   2. Select hydrophobic protein residues.
#   3. Identify contacts using VMD "measure contacts".
#   4. Identify protein residues involved in contacts.
#   5. Count each residue only once per frame.
#   6. Calculate occupancy:
#
#      Occupancy (%) =
#      (Number of contacted frames / Total frames) × 100
#
# Requirements:
#   - VMD
#   - Protein topology file
#   - MD trajectory
#
# Run:
#   vmd -dispdev text -e hydrophobic_occupancy.tcl
#
# ============================================================


# ============================================================
# USER INPUT
# Modify only this section
# ============================================================

# Topology file
set topology "YOUR_TOPOLOGY_FILE.prmtop"

# Trajectory file
set trajectory "YOUR_TRAJECTORY_FILE.dcd"

# Output file
set output_file "hydrophobic_occupancy.dat"

# Contact cutoff in Angstrom
set cutoff 4.5

# Ligand name
set ligand_resname "YOUR_LIGAND_RESNAME"

# Ligand atoms to analyze
# Enter atom names separated by spaces
set ligand_atoms "ATOM1 ATOM2 ATOM3"

# Hydrophobic protein residues
set hydrophobic_residues "ALA VAL LEU ILE MET PHE TRP PRO"


# ============================================================
# LOAD TOPOLOGY
# ============================================================

puts "=========================================="
puts "Hydrophobic Contact Occupancy Analysis"
puts "=========================================="

puts ""
puts "Loading topology..."

mol new $topology


# ============================================================
# LOAD TRAJECTORY
# ============================================================

puts "Loading trajectory..."

mol addfile $trajectory waitfor all

puts "Trajectory loaded."


# ============================================================
# NUMBER OF FRAMES
# ============================================================

set nframes [molinfo top get numframes]

puts "Number of frames = $nframes"

if {$nframes == 0} {
    puts "ERROR: No frames found in the trajectory."
    exit
}


# ============================================================
# CHECK LIGAND SELECTION
# ============================================================

set lig_test [atomselect top \
    "resname $ligand_resname and name $ligand_atoms"]

puts ""
puts "Ligand residue name : $ligand_resname"
puts "Ligand atoms        : $ligand_atoms"
puts "Ligand atoms found  : [$lig_test num]"

if {[$lig_test num] == 0} {

    puts ""
    puts "ERROR: No ligand atoms found."
    puts "Please check the ligand residue name and atom names."

    $lig_test delete
    exit
}

$lig_test delete


# ============================================================
# CONTACT COUNTER
# ============================================================

array set counts {}


# ============================================================
# START ANALYSIS
# ============================================================

puts ""
puts "Contact cutoff = $cutoff Angstrom"
puts "Starting analysis..."
puts ""


# ============================================================
# FRAME-BY-FRAME ANALYSIS
# ============================================================

for {set f 0} {$f < $nframes} {incr f} {

    # Move to current frame
    animate goto $f


    # --------------------------------------------------------
    # Select ligand atoms
    # --------------------------------------------------------

    set lig [atomselect top \
        "resname $ligand_resname and name $ligand_atoms"]


    # --------------------------------------------------------
    # Select hydrophobic protein residues
    # --------------------------------------------------------

    set prot [atomselect top \
        "protein and resname $hydrophobic_residues"]


    # --------------------------------------------------------
    # Calculate contacts
    # --------------------------------------------------------

    set contacts [measure contacts $cutoff $lig $prot]


    # --------------------------------------------------------
    # Get protein atom indices involved in contacts
    # --------------------------------------------------------

    set prot_idx [lsort -unique [lindex $contacts 1]]


    # --------------------------------------------------------
    # Store residues contacted in this frame
    # --------------------------------------------------------

    array unset frame_res
    array set frame_res {}


    foreach idx $prot_idx {

        set atm [atomselect top "index $idx"]

        set resid   [lindex [$atm get resid] 0]
        set resname [lindex [$atm get resname] 0]

        set key "${resid}_${resname}"

        set frame_res($key) 1

        $atm delete
    }


    # --------------------------------------------------------
    # Count each residue only once per frame
    # --------------------------------------------------------

    foreach key [array names frame_res] {

        if {[info exists counts($key)]} {
            incr counts($key)
        } else {
            set counts($key) 1
        }
    }


    # --------------------------------------------------------
    # Delete selections
    # --------------------------------------------------------

    $lig delete
    $prot delete


    # --------------------------------------------------------
    # Progress message
    # --------------------------------------------------------

    if {$f % 1000 == 0} {
        puts "Frame $f / $nframes"
    }
}


# ============================================================
# WRITE RESULTS
# ============================================================

puts ""
puts "Writing results..."

set outfile [open $output_file w]

puts $outfile "Residue\tFrames\tOccupancy(%)"


foreach key [lsort [array names counts]] {

    set occ [expr {100.0 * $counts($key) / $nframes}]

    puts $outfile \
        "$key\t$counts($key)\t[format %.2f $occ]"
}


close $outfile


# ============================================================
# FINISHED
# ============================================================

puts ""
puts "=========================================="
puts "Analysis completed successfully."
puts "=========================================="
puts "Total frames analyzed : $nframes"
puts "Contact cutoff        : $cutoff Angstrom"
puts "Output file           : $output_file"
puts "=========================================="

exit
