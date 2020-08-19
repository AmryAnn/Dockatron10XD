__author__ = "jwhittle"

import csv
from pathlib import Path
import numpy as np
import argparse
import logging



# LPATH = "C:\BestPocket\PCSK9_Sample_Set_1\PDB_ligand_poses"
# RPATH = "C:\BestPocket\PCSK9_Sample_Set_1\\receptor_pdbs"

def compute_Distance(x1, x2):
    """ Compute the distance between two np arrays
        param: x1, x2 - np arrays of the same size
        return: np array (x,y,z)
    """
    return np.sqrt(np.sum((x2-x1)**2))



def find_LigandCenter(ligFile):
    """ Find the center of a ligand from the pdb file
        param: ligFile - Path object from pathlib
        return: np array (x,y,z)
    """
    coords = np.empty((0, 3), float)
    try:
        with open(ligFile) as pdbFile:
            for line in pdbFile:
                if line.startswith("HETATM"):
                    atomCoords = list(np.float_(line.split()[5:8]))
                    coords = np.append(coords, np.array([atomCoords]), axis=0)
                elif line.startswith("ATOM"):
                    atomCoords = list(np.float_(line.split()[5:8]))
                    coords = np.append(coords, np.array([atomCoords]), axis=0)
                else:
                    continue
    except FileNotFoundError:
        logging.warning(f"Ligand file not found - {ligFile}")
    centroid = coords.mean(axis=0)
    return centroid


def read_LigandAffinity(affFile):
    """ Read the binding affinities and other info from a VINA log file
        param: affFile - Path object from pathlib
        return: list [pose, affinity, rmsd lb, rmsd ub]
    """
    poseLines = []
    try:
        with open(affFile) as vinaFile:
            for line in vinaFile:
                if line.startswith("--"):
                    nextLine = vinaFile.readline()
                    while not nextLine.startswith("Writing"):
                        poseLines.append(nextLine.strip("\n").split())
                        nextLine = vinaFile.readline()
    except FileNotFoundError:
        logging.warning(f"Vina log file not found - {affFile}")
    return poseLines


def create_LigandPoseDict(ligPath):
    """ Step through all ligand pdb files and
            create dictionary of centers and binding info for each
        param: ligPath - Path object from pathlib
        return: dict - indexed by ligName
    """
    ligandDict = {}
    dockData = {}
    ligFiles = list(ligPath.glob('**/*.pdb'))
    logging.info(f"Found {len(ligFiles)} ligand files")

    for lFile in ligFiles:
        ligName = "".join(("L", "_".join(Path(lFile).stem.split("\\")[-1].split("_")[1:4])))
        dockName = "_".join(ligName.split("_")[:2])
        aFile = Path(lFile).parents[0].joinpath("".join((dockName[1:], "_log.txt")))

        ligCentroid = find_LigandCenter(lFile)
        if dockName not in dockData.keys():  # Check for already read ligand affinities
            dockData[dockName] = read_LigandAffinity(aFile)
        poseValues = [ligCentroid]
        if len(dockData[dockName]):
            poseNum = int(ligName[-1]) - 1
            poseValues.extend(dockData[dockName][poseNum])
        else:
            poseValues.extend(["NA", "NA", "NA"])
        ligandDict[ligName] = poseValues
    logging.info(f"Computed {len(ligandDict)} ligand centers")
    print(f"Computed {len(ligandDict)} ligand centers")
    return ligandDict


def read_PocketRanges(recPath):
    """ Read pocket_residues.txt file and
            create list of pocket residues for each pocket
        param: recPath - Path object from pathlib
        return: dict - indexed by pocketNum
    """
    pocketRange = {}
    pocketResFile = recPath.joinpath("pocket_residues.txt")
    try:
        with open(pocketResFile) as pFile:
            pocketLines = pFile.readlines()
            for pLine in pocketLines:
                pocketNum = int(pLine.split(",")[0])
                pRanges = pLine.split(",")[1:]
                newPocket = []
                for rng in pRanges:
                    start = int(rng.split("-")[0])
                    end = int(rng.split("-")[1])
                    newPocket.extend(list(range(start, end+1)))
                pocketRange[pocketNum] = newPocket
    except FileNotFoundError:
        logging.error(f"Pocket residues file not found - {pocketResFile}")
        exit()
    logging.info(f"Found {len(pocketRange.keys())} receptor binding sites")
    return pocketRange


def remove_Backbone(recFile):
    """ Read receptor pdb file and remove hydrogens and backbone atoms
        param: recFile - Path object from pathlib
        return: list - side chain atoms
    """
    hydrogens = ['H1', 'H2', 'HA', 'HB', 'HD', 'HE', 'HG', 'HH', 'HZ']
    backbone = ['H', 'N', 'C', 'O']
    sideChains = []
    try:
        with open(recFile) as pdbFile:
            for line in pdbFile:
                if line.startswith("END"):
                    sideChains.append(line)
                    break
                atom = line.split()[2]
                discard = len([a for a in hydrogens if a in atom]) + \
                          len([a for a in backbone if atom == a])
                if not discard:
                    sideChains.append(line)
            logging.debug(f"Found {len(sideChains)} side chains")
    except FileNotFoundError:
        logging.warning(f"Receptor file not found - {recFile}")

    scFileName = RPATH.joinpath("".join((recFile.stem, "_sidechains.txt")))
    try:
        with open(scFileName, "w") as scFile:
            for atom in sideChains:
                scFile.write("{}".format(atom))
    except FileNotFoundError:
        logging.warning(f"Unable to write to {scFileName}")
    return sideChains


def assign_AtomsToPockets(recFile, sideChain, pocketPos):
    """ Read receptor pdb file and remove backbone atoms
        param: recFile - Path object from pathlib
               sideChain - list of sidechain atoms in the receptor
               pocketPos - list of residues associated with each pocket
        return: list - atoms (pdb format) in the pockets with pocket number added
    """
    pocketAtoms = []
    for aCount, atom in enumerate(sideChain):
        pocketNum = [pnum for pnum, positions in pocketPos.items() if aCount in positions]
        if pocketNum:
            atom = "\t".join((atom[:-1], str(pocketNum[0])))
            pocketAtoms.append(atom)
    logging.debug(f"Found {len(pocketAtoms)} binding pocket atoms")

    siteFile = recFile.parents[0].joinpath("".join((recFile.stem, "_pocket_residues.txt")))
    try:
        with open(siteFile, "w") as siteFile:
            for atom in pocketAtoms:
                siteFile.write("{}\n".format(atom))
    except FileNotFoundError:
        logging.warning(f"Unable to write to {siteFile}")
    return pocketAtoms


def find_PocketCenters(pockets):
    """ Create dictionary of the center of each pocket in a receptor
        param: pockets - list of pocket atoms (pdb format)
        return: dict - np array (x,y,z) indexed by pocketNum
    """
    pockets.sort(key=lambda pos: int(pos.rsplit("\t", 1)[1]))
    centroids = {}
    curID = 1
    coords = np.empty((0,3), float)
    for line in pockets:
        newID = int(line.rsplit("\t", 1)[1])
        if newID != curID:
            centroids[curID] = coords.mean(axis=0)
            curID = newID
            coords = np.empty((0, 3), float)
        atomCoords = list(np.float_(line.split()[6:9]))
        coords = np.append(coords, np.array([atomCoords]), axis=0)
    centroids[curID] = coords.mean(axis=0)
    return centroids


def read_PreviousPocketCenters(pCentersFile):
    """ Read file from previous run with already computed receptor pocket centers
        param: pCenterFile - Path object from pathlib
        return: dict - np array (x,y,z) indexed by receptorNum, pocketNum
    """
    pocketDict = {}
    try:
        with open(pCentersFile) as ctrFile:
            curRec = -1
            for line in ctrFile:
                ctrData = line.split(",")
                recNum = int(ctrData[0])
                coords = np.array(list(np.float_(ctrData[2:5])))
                if recNum != curRec:
                    pocketDict[recNum] = {int(ctrData[1]): coords}
                    curRec = recNum
                else:
                    pocketDict[recNum][int(ctrData[1])] = coords
            logging.info(f"Found {len(pocketDict)} receptors with already calculated binding site centers")
    except FileNotFoundError:
        logging.warning(f"Previous pocket centers file not found - {pCentersFile}")
    return pocketDict


def create_RecptorPocketsDict(recPath, pocketPos):
    """ Step through all receptor files (pdb) and
            create dictionary with their pocket centers
        param: recPath - Path object from pathlib
               pocketPos - list of residues for each binding pocket of the receptor
        return: dict - np array (x,y,z) indexed by receptorNum, pocketNum
    """
    recFiles = list(recPath.glob('**/*.pdb'))
    logging.info(f"Found {len(recFiles)} receptor files")
    print("Found", len(recFiles), "receptor files")
    numPockets = len(pocketPos)
    pocketCentersFile = recPath.joinpath("rec_pocket_centers.txt")
    savedPocketDict = {}
    savedRecIDs = []
    if pocketCentersFile.exists():
        savedPocketDict = read_PreviousPocketCenters(pocketCentersFile)
        savedRecIDs = savedPocketDict.keys()

    rPocketDict = {}
    for rFile in recFiles:
        logging.debug(f"Processing {rFile}")
        recName = int(rFile.stem)
        if recName in savedRecIDs and len(savedPocketDict[recName]) == numPockets:
            rPocketDict[recName] = savedPocketDict[recName]
        else:
            sideChains = remove_Backbone(rFile)
            pocketResidues = assign_AtomsToPockets(rFile, sideChains, pocketPos)
            rPocketDict[recName] = find_PocketCenters(pocketResidues)

        try:
            with open(str(pocketCentersFile), "w") as ctrFile:
                for recNum in rPocketDict:
                    for pNum in rPocketDict[recNum]:
                        coords = [str(val) for val in rPocketDict[recNum][pNum]]
                        pocketStr = ",".join((str(pNum), ",".join(coords[:3])))
                        ctrFile.write("{0},{1}\n".format(recNum, pocketStr))
        except FileNotFoundError:
            logging.warning(f"Unable to write to {pocketCentersFile}")
    return rPocketDict


def assign_PocketsToLigands(ligCenters, recCenters):
    """ Assign each ligand to a binding pocket of each receptor
        param: ligCenters - dict of ligand centers
               recCenters - dict of receptor centers
        return: dict - binding site number (pocketNum) indexed by ligandID
    """
    ligSiteDict = {}
    for recID in recCenters:
        ligMatch = [lig for lig in ligCenters if int(lig.split("_")[1]) == recID]
        for ligID in ligMatch:
            distances = [compute_Distance(ligCenters[ligID][0], recCenters[recID][siteID])
                         for siteID in recCenters[recID].keys()]
            ligSiteDict[ligID] = distances.index(min(distances)) + 1

    logging.info(f"Assigned {len(ligSiteDict)} ligands to binding sites")
    print(f"Assigned {len(ligSiteDict)} ligands to binding sites")
    assigned = set(ligSiteDict.keys())
    diff = [L_id for L_id in ligCenters.keys() if L_id not in assigned]
    logging.warning(f"Ligand IDs not assigned to binding pockets {diff}")
    return ligSiteDict


def write_LigandBindingSiteFile(ligPath, lSites, lValues):
    """ Write the binding pocket,affinity and binding info for each ligand
        param: ligPath - Path object from pathlib
               lSites - dict of assigned binding sites by ligandID
               lValues - dict of binding affinity and info by ligandID
        return: dict - binding site number (pocketNum) indexed by ligandID
    """
    bindFile = ligPath.joinpath("ligand_binding_sites.txt")
    try:
        with open(bindFile, "w", newline="") as lsFile:
            ligWriter = csv.writer(lsFile, delimiter="\t")
            header = ["Lig", "Rec", "Pose", "Site", "Energy", "RMSD LB", "RMSD UB"]
            ligWriter.writerow(header)
            for ligID in lSites:
                ligInfo = ligID.split("_")
                ligInfo.extend([lSites[ligID]])
                ligInfo.extend(lValues[ligID][1:])
                ligWriter.writerow(ligInfo)
            logging.info(f"Saved ligand binding sites to {bindFile}")
    except FileNotFoundError:
        logging.error(f"Unable to write to {bindFile}")
        exit()
    return


if __name__ == "__main__":
    # Add command line argument parser
    parser = argparse.ArgumentParser(description="Assign a receptor binding site to ligands from pdb docking data.\n\n"
                                                 "Ligand and receptor files must be in pdb format.\n"
                                                 "Ligand file names in format -- Ligand_<lignum>_<recnum>_<posenum>.pdb\n"
                                                 "Receptor file names in the format -- <recnum>.pdb\n"
                                                 "Ligand path should include Vina log files with binding affinities.\n"
                                                 "Receptor path should include pocket_residues.txt file.\n"
                                                 "See readme.txt for more information on files and naming expectations.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    required.add_argument("--l", type=str, required=True,
                        help="The path for the ligand pdb and Vina log files")
    required.add_argument("--r", type=str, required=True,
                            help="The path for the receptor pdb files")
    optional.add_argument("--v", type=str, default="N", help="Verbose logging (Y/N)")
    parser._action_groups.append(optional)


    # Add logger to create a log file
    logging.basicConfig(filename="assign_binding_sites.log", format="%(asctime)s - %(levelname)s - %(message)s",
                        level=logging.NOTSET)
    logger = logging.getLogger(__name__)
    logging.info("*****Begin run of assign_binding_sites.py*****")

    # Assign command line arguments
    args = parser.parse_args()
    LPATH = Path(args.l)
    if not LPATH.exists() or not LPATH.is_dir():
        logging.error(f"Ligand path {LPATH} is not a valid directory")
        exit()
    RPATH = Path(args.r)
    if not RPATH.exists() or not RPATH.is_dir():
        logging.error(f"Receptor path {RPATH} is not a valid directory")
        exit()

    # Set logging level (not working on Win10, Python 3.7)
    if args.v.capitalize() in ("Y", "YES"):
        level = logging.INFO
        logging.info("Running in verbose mode.")
    else:
        level = logging.WARNING
    logger.setLevel(level)
    for handler in logger.handlers:
        handler.setlevel(level)


    # Call main functions of the script
    ligandValues = create_LigandPoseDict(LPATH)
    pocketPositions = read_PocketRanges(RPATH)
    receptorCenters = create_RecptorPocketsDict(RPATH, pocketPositions)
    ligandSites = assign_PocketsToLigands(ligandValues, receptorCenters)
    write_LigandBindingSiteFile(LPATH, ligandSites, ligandValues)
    logging.info("Finished assigning ligand binding sites")
    print("Finished assigning ligand binding sites")
