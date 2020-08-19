
# 2. Load MD Frames - 	source Load_MD.tcl
menu files on
display resetview
mol new {test.dcd} type {dcd} first 0 last -1 step 1 waitfor 1
animate style Loop
mol addrep 0
display resetview
mol addfile {14_protein.pdb} type {pdb} first 0 last -1 step 1 waitfor 1 0
animate style Loop
display resetview
menu files off

# 3. Measure RMSD - 	source measureRMSD.tcl
set outfile [open rmsd.dat w]
set nf [molinfo top get numframes]
set frame0 [atomselect top "protein and backbone and noh" frame 0]
set sel [atomselect top "protein and backbone and noh"]
for {set i 1} {$i <= $nf} {incr i} {
	$sel frame $i
	$sel move [measure fit $sel $frame0]
	puts $outfile "[measure rmsd $sel $frame0]"
}
close $outfile

# 4. Measure Center - 	source measureCenter.tcl
set outfile [open center.dat w]
set nf [molinfo top get numframes]
set frame0 [atomselect top protein frame 0]
set sel [atomselect top protein]
for {set i 1} {$i <= $nf} {incr i} {
	$sel frame $i
	puts $outfile "[measure center $sel]"
}
close $outfile

# 5. Measure min Max - 	source Measure_minMax.tcl
set outfile [open minMax.dat w]
set nf [molinfo top get numframes]
set frame0 [atomselect top protein frame 0]
set sel [atomselect top protein]
for {set i 1} {$i <= $nf} {incr i} {
	$sel frame $i
	puts $outfile "[measure minmax $sel]"
}
close $outfile

# 6. Write coordinate files - 	source allPDBs.tcl
mkdir ReceptorPDBQTs
cd ReceptorPDBQTs
for {set i 1} {$i < 23} {incr i} { 
         [atomselect top all frame $i] writepdb $i.pdb 
 } 
 
 #7. Exit
 exit