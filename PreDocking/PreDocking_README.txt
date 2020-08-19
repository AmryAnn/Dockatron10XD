PreDocking Procedures

Part 0

Copy the Dockatron10XD package and rename the folder after your own project.
Download and install all necessary software.

Part 1

1. Move your trajectory file (.dcd) and input coordinate file (pdb) into MD_data_from_VMD
2. Create 'LigandPDBs' and 'ReceptorPDBs' directories in /Dockatron10XD/PreDocking/
3. Move your Ligand PDB files to vmd /Dockatron10XD/PreDocking/LigandPDBs/
4. Edit MD_data_from_VMD.tcl so the dcd filename (line 5) and the pdb filename (line 9) match your filenames.
5. From the command line (WSL/Bash), move to /Dockatron10XD/PreDocking/ and enter the following command:
	
	vmd -dispdev text -eofexit < MD_data_from_VMD.tcl > MD_data_from_VMD.log
	
Expected Output:
LigandPDBs/		-containing your ligand PDB files
Step1A.log		-VMD logging output
center.dat		-center xyz coordinates of protein in each frame of trajectory
minMax.dat		-minimum and maximum xyz coordinates of protein in each frame of trajectory
rmsd.dat		-rmsd values for each frame in the trajectory
ReceptorPDBs/	-containing the receptor pdb files created from running Step1A.tcl



Part 2

1. From /Dockatron10XD/PreDocking/, enter the following command into the command line:

	python Step2.py --rmsd_dat 'full/path/to/rmsd.dat' --minMax_dat 'full/path/to/minMax.dat' --center_dat 'full/path/to/center.dat' --extra <integer value>

2. Move all .conf files to /Docking_Job/
3. Move all ligand pdbqt files to /Docking_Job/
4. Move all receptor pdbqt files to /Docking_Job/

Expected Output:
RMSD Graph				-pops up in separate window
rand10 					-randomly selected frame numbers (printed to command line - aid in 'clustered' structure selection)
.conf files				-configuration files for each receptor pdb
Receptor PDBQT files 	-in ReceptorPDBQTs/ directory
Ligand PDBQT files		-in LigandPDBQTs/ directory
