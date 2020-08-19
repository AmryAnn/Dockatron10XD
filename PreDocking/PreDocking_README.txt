PreDocking Procedures

Part 0

Copy the Dockatron10XD package and rename the folder after your own project.
Download and install all necessary software.

Part 1

1. Move your trajectory file (.dcd) and input coordinate file (pdb) into MD_data_from_VMD
2. Create 'LigandPDBQTs' and 'ReceptorPDBQTs' directories in /Dockatron10XD/PreDocking/
3. Move your Ligand PDBQT files to /Dockatron10XD/PreDocking/LigandPDBQTs/
4. Edit MD_data_from_VMD.tcl so the dcd filename (line 5) and the pdb filename (line 9) match your filenames.
5. From the command line (WSL/Bash), move to /Dockatron10XD/PreDocking/ and enter the following command:
	
	vmd -dispdev text -eofexit < MD_data_from_VMD.tcl > MD_data_from_VMD.log
	
6. From /ReceptorPDBQTs/, convert your Receptor PDB files to PDBQT files using Open Babel or AutoDockTools
7. Move your Receptor PDB files to archive
	
Expected Output:
LigandPDBQTs/		-containing your ligand PDBQT files
Step1A.log		-VMD logging output
center.dat		-center xyz coordinates of protein in each frame of trajectory
minMax.dat		-minimum and maximum xyz coordinates of protein in each frame of trajectory
rmsd.dat		-rmsd values for each frame in the trajectory
ReceptorPDBQTs/	-containing the receptor pdb files created from running Step1A.tcl


Part 2

1. From /Dockatron10XD/PreDocking/, enter the following command into the command line:

	python Step2.py --rmsd_dat 'full/path/to/rmsd.dat' --minMax_dat 'full/path/to/minMax.dat' --center_dat 'full/path/to/center.dat' --extra <integer value>

2. Copy all .conf files to /Docking_Job/
3. Copy all ligand pdbqt files to /Docking_Job/
4. Copy all receptor pdbqt files to /Docking_Job/

Expected Output:
RMSD Graph				-pops up in separate window
rand10 					-randomly selected frame numbers (printed to command line - aid in 'clustered' structure selection)
.conf files				-configuration files for each receptor pdb
Receptor PDBQT files 	-in /ReceptorPDBQTs/ and in /Docking_Job/ directories
Ligand PDBQT files		-in /LigandPDBQTs/ and in /Docking_Job/ directories
