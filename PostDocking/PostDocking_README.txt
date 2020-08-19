PostDocking Procedurs

Part 3
1. Generate pocket_residus.txt
2. Assign Binding Sites
	a. The script is assign_binding_sites.py
		The command line parameters are:
			--l Ligand directory 
			--r Receptor directory
	
	b. The ligand directory contains the ligand pdb files named 
			Ligand_<Num>_<Rec>_Pose.pdb *
			* These can be in sub directories
	
	c. The receptor directory contains the receptor conformation pdb files named
			<Rec>.pdb
			and a pocket_residues.txt file
	
	f. The pocket_residues files is in the following format:
		<pocket num>,<start residue number>-<end residue number>,*
			* repeat residue ranges as needed
	
	g. Temporary output files are created for inspection in the receptor directory
		<Rec>_sidechains.txt		a file with the side chains extracted from the pdb file
		<Rec>_binding_sites.txt		a file with the atoms that are part of the binding sites
						binding site number is added as last column with \t
		rec_pocket_centers.txt		a file with computed binding centers for each pocket
						the script will check this file against the receptor files in the path
						and skip computations for the binding sites if the record is in this file
		
	
	h. The final output is in the ligand directory
		ligand_binding_sites.txt	a tab separated file with header and columns <Ligand Pose Receptor BindSite>
	
	i. A log file is saved in the directory that the script is stored.	
	
	Example command:
	python assign_binding_sites.py --l 'C:\Users\aarid\Set3_Results\ligands\' --r 'C:\Users\aarid\Set3_Results\receptors'
	
3. Ligand Binding Sites Analysis
	a. The script is Ligand_Sites_Analysis.py 
		The command line parameters are:
			--ligand_sites	'<path to ligand_binding_sites.txt>'
			--ligand_ID		'ligand identifier'
	b. Two windows containing graphs will pop up consecutively -- save those!
	c. A loge file containing the same information that is printed on the console is saved in the directory that the script is stored.
	
	Example command:	
	python Ligand_Sites_Analysis.py --ligand_sites 'C:\Users\aarid\Set3_Results\ligand_binding_sites5_5000.txt' --ligand_ID 'L23'