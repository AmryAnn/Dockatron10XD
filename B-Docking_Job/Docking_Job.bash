#!/bin/bash/
# Create docking set directories and move vina configuration files and receptor pdbqt 
#files into the corresponding directories
for x in ./*.conf; do mkdir "${x%.*}" && mv "$x" "${x%.*}";done
for x in ./*.pdbqt; do mv "$x" "${x%.*}";done
# Move ligand pdbqts into each receptor directory
for dir in */; do cp Ligand_01_Mol1.pdbqt "$dir";done
for dir in */; do cp Ligand_02_39401.pdbqt "$dir";done
for dir in */; do cp Ligand_03_35787.pdbqt "$dir";done
for dir in */; do cp Ligand_04_10926.pdbqt "$dir";done
for dir in */; do cp Ligand_05_Mol2.pdbqt "$dir";done