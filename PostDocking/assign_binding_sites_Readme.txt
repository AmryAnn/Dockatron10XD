Assign Binding Site 
1. The script is assign_binding_pocket.py
	The command line parameters are:
		--l Ligand directory 
		--r Receptor directory

2. The ligand directory contains the ligand pdb files named 
		Ligand_<Num>_<Rec>_Pose.pdb *
		* These can be in sub directories

3. The receptor directory contains the receptor conformation pdb files named
		<Rec>.pdb
		and a pocket_residues.txt file

4. The pocket_residues files is in the following format:
	<pocket num>,<start residue number>-<end residue number>,*
		* repeat residue ranges as needed

5. Temporary output files are created for inspection in the receptor directory
	<Rec>_sidechains.txt		a file with the side chains extracted from the pdb file
	<Rec>_binding_sites.txt		a file with the atoms that are part of the binding sites
					binding site number is added as last column with \t
	rec_pocket_centers.txt		a file with computed binding centers for each pocket
					the script will check this file against the receptor files in the path
					and skip computations for the binding sites if the record is in this file
	

6. The final output is in the ligand directory
	ligand_binding_sites.txt	a tab separated file with header and columns <Ligand Pose Receptor BindSite>

7. A log file is saved in the directory that the script is stored.
		
	