# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:52:39 2020

@author: aarid"""


import pandas as pd
import random 
from matplotlib import pyplot as plt
from pathlib import Path
import argparse
import logging


""" Convert MD frame number to nanoseconds, print 2 list of random frame 
numbers, and plot the RMSD"""

def rmsd(rmsd_file):
    rmsd_new = open('rmsd.txt', 'w')
    #write header to new file
    rmsd_head = 'rmsd\n'
    rmsd_new.write(rmsd_head)
    
    #write rmsd values by line to new file.
    with open(rmsd_file) as rmsd_dat:
        for Line in rmsd_dat: 
            rmsd_line = '%s' %Line
            rmsd_new.write(rmsd_line)
            #print(Line)
        rmsd_new.close()
        
        df = pd.read_csv('rmsd.txt')
        df['picoseconds'] = df.index*25
        df['nanoseconds'] = (df.index*25)/1000
        df['Frame'] = range(1, len(df)+1)
        rand10 = random.sample(range(1, len(df.Frame)+1), k=10)
        #rand20 = random.sample(range(len(df.index)+1), k=20)
        print('\n', df.head(23), '\n')
        logging.info(f'\n{df}')
        print('Frames for clustered docking:', rand10)
        logging.info(f'Frames for clustered docking: {rand10}')
        #print(rand20)
        
        plt.plot(df.nanoseconds, df.rmsd)
        plt.xlabel('Time (ns)')
        plt.ylabel('RMSD')
        plt.show()
        return df, rand10 #, rand20
    

""" Calculate the size of the docking box. """

def box_size(minMax_file, extra):
    f_min = open('NewMin.txt', 'w')
    f_max = open('NewMax.txt', 'w')
    header_min = 'minX'+'\t'+'minY'+'\t'+'minZ'+'\n'
    header_max = 'maxX'+'\t'+'maxY'+'\t'+'maxZ'+'\n'
    f_min.write(header_min)
    f_max.write(header_max)
    with open(minMax_file) as minMax_dat:
        LineNumber = 0
        for Line in minMax_dat:
            if LineNumber >= 0:
                Line = Line.strip()
                Line = Line.split()
                #print(Line[1])
                fixd0 = str(Line[0]).replace("{", "")
                fixd2 = str(Line[2]).replace("}", "")
                fixd3 = str(Line[3]).replace("{", "")
                fixd5 = str(Line[5]).replace("}", "")
                NewMin = []
                NewMax = []
                NewMin.append(fixd0)
                NewMin.append(Line[1])
                NewMin.append(fixd2)
                NewMin0 = str(NewMin[0])
                NewMin1 = str(NewMin[1])
                NewMin2 = str(NewMin[2])
                NewMin = NewMin0+'\t'+NewMin1+'\t'+NewMin2
                #print(NewMin)
                f_min.write(NewMin+'\n')
                NewMax.append(fixd3)
                NewMax.append(Line[4])
                NewMax.append(fixd5)
                NewMax0 = str(NewMax[0])
                NewMax1 = str(NewMax[1])
                NewMax2 = str(NewMax[2])
                NewMax = NewMax0+'\t'+NewMax1+'\t'+NewMax2
                #print(NewMax)
                f_max.write(NewMax+'\n')
                LineNumber = LineNumber + 1
            LineNumber = LineNumber + 1            
        f_min.close()
        f_max.close()
        """ Determine global Min and Max Values, find the difference for each
        coordinate value and add buffer distance (extra)"""    
        min_data = pd.read_csv("NewMin.txt", sep='\t')
        max_data = pd.read_csv("NewMax.txt", sep='\t')
        extra = extra    
        box = open('box_size.txt', 'w')    
        #print(min_data.head(25))
        #print(max_data.head(25))    
        minX = min_data.minX.min()
        minY = min_data.minY.min()
        minZ = min_data.minZ.min()
        #print(minX, minY, minZ)    
        maxX = max_data.maxX.max()
        maxY = max_data.maxY.max()
        maxZ = max_data.maxZ.max()
        #print(maxX, maxY, maxZ)    
        recX = maxX - minX
        recY = maxY - minY
        recZ = maxZ - minZ
        #print(recX, recY, recZ)    
        sizeX = recX + extra
        sizeY = recY + extra
        sizeZ = recZ + extra
        #print(sizeX, sizeY, sizeZ)    
        sizeX = str(sizeX)
        sizeY = str(sizeY)
        sizeZ = str(sizeZ)
        sizeBox = sizeX+'\t'+sizeY+'\t'+sizeZ
        box.write(sizeBox)
        print('Size of docking box:',sizeBox, '\n')   
        logging.info(f'Size of docking box: {sizeBox}')
        box.close()
        return sizeBox


""" Write Vina configuration files for each receptor PDB using 'center.dat'
as input"""

def Vina_configuration_files(center_file, num_modes, cpu, exhaustiveness):
    oc = open('center.tsv', 'w')
    LineNumber = 0
    with open(center_file) as center_dat:
        for Line in center_dat:
            if LineNumber >= 0:
                Line = Line.strip()
                Line = Line.split()
                #print(Line)
                Line0 = str(Line[0])
                Line1 = str(Line[1])
                Line2 = str(Line[2])
                Line = Line0+'\t'+Line1+'\t'+Line2
                oc.write(Line+'\n')
            LineNumber = LineNumber + 1  
        oc.close()
        
        """ Add 'Index' column to center.tsv """
        t = pd.read_csv('center.tsv', sep='\t')
        t['Index'] = range(1, len(t)+1)
        #print(t.head())
        t.to_csv('centerD.tsv', sep='\t')
        
        """ Write the configuration files"""
        fcd=open('centerD.tsv', 'r')
        b=open('box_size.txt', 'r')
        
        LineNumber = 0
        for Lineb in b:
            Lineb = Lineb.split()
            #print(Lineb[0])
        for Line in fcd:
            if LineNumber > 0:
                Line = Line.strip('\n')
                RecList = Line.split()
                #print(RecList)
                oTitle = "%s.conf" %LineNumber
                ocd = open(oTitle, 'w')
                intro = "#Configuration file"
                ocd.write(intro+"\n\n")
                oReceptor = "receptor = %s.pdbqt" %RecList[4]
                ocd.write(oReceptor+"\n\n")
                oCenter_x = "center_x = %s" %RecList[1]
                ocd.write(oCenter_x+"\n")
                oCenter_y = "center_y = %s" %RecList[2]
                ocd.write(oCenter_y+"\n")
                oCenter_z = "center_z = %s" %RecList[3]
                ocd.write(oCenter_z+"\n\n")
                oSize_x = "size_x = %s" %Lineb[0]
                oSize_y = "size_y = %s" %Lineb[1]
                oSize_z = "size_z = %s" %Lineb[2]
                oSize = oSize_x+"\n"+oSize_y+"\n"+oSize_z+"\n\n"
                ocd.write(oSize)
                oModes = "num_modes = %s" %num_modes
                oCpu = "cpu = %s" %cpu
                oExh = "exhaustiveness = %s" %exhaustiveness
                ocd.write(oModes+"\n"+oCpu+"\n"+oExh+"\n")
            LineNumber = LineNumber + 1
        fcd.close()
        b.close()
        ocd.close()
        return oTitle


""" Command line arguments and logfile """

if __name__ == '__main__':
    # Add command line argument parser
    parser = argparse.ArgumentParser(description="Run Step 2 of Dockatron10XD.\n\n"
                                                 "Process RMSD data.\n"
                                                 "Calculate the size of docking box.\n"
                                                 "Write Vina configuration files.", \
                                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--rmsd_dat', type=str, required=True,\
                            help='The path to rmsd.dat')
    parser.add_argument('--minMax_dat', type=str, required=True,\
                            help='The path to minMax.dat')
    parser.add_argument('--center_dat', type=str, required=True,\
                            help='The path to center.dat')
    parser.add_argument('--extra', type=int, required=True, \
                            help='Buffer distance for docking box (angstroms)')
    parser.add_argument('--num_modes', type=int, required=True, \
                            help='Number of different poses desired for each ligand-receptor pair.')
    parser.add_argument('--cpu', type=int, required=True, \
                            help='Number of cores to use for each docking set.')
    parser.add_argument('--exhaust', type=int, required=True, \
                            help='Level of exhaustiveness to use in Vina')
    
    # Add logger to create a log file
    logging.basicConfig(filename="Step2.log", format="%(asctime)s - %(levelname)s - %(message)s",
                        level=logging.NOTSET)
    logging.getLogger('matplotlib.font_manager').disabled = True
    logger = logging.getLogger(__name__)
    logging.info("*****Begin run of Dockatron10XD Step 2*****")
    
    # Assign command line arguments
    args = parser.parse_args()
    RMSD_PATH = Path(args.rmsd_dat)
    if not RMSD_PATH.exists() or not RMSD_PATH.is_file():
        logging.error(f"rmsd.dat path {RMSD_PATH} is not a valid file or directory")
        exit()
    MINMAX_PATH = Path(args.minMax_dat)
    if not MINMAX_PATH.exists() or not MINMAX_PATH.is_file():
        logging.error(f"minMax.dat path {MINMAX_PATH} is not a valid file or directory")
        exit()
    CENTER_PATH = Path(args.center_dat) 
    if not CENTER_PATH.exists() or not CENTER_PATH.is_file():
        logging.error(f"center.dat path {RMSD_PATH} is not a valid file or directory")
        exit()
    EXTRA = args.extra
    if EXTRA:
        logging.info(f'{EXTRA} buffer angstroms added to each side of docking box.')
        print('\n', EXTRA, 'buffer angstroms added to each side of docking box.')
    MODES = args.num_modes
    if MODES:
        logging.info(f'Output will include {MODES} ligand poses per ligand-receptor pair.')
        print('Output will include', MODES, 'ligand poses per ligand-receptor pair.')
    CPU = args.cpu
    if CPU:
        logging.info(f'Set {CPU} cores per docking calculation.')
        print(CPU, 'cores per docking calculation.')
    EXH = args.exhaust
    if EXH:
        logging.info(f'Docking ligand to receptor with exhaustiveness of {EXH}.')
        print('Docking ligand to receptor with exhaustiveness of', EXH)
    # Set logging level

    # Function calls
    rmsd(RMSD_PATH)
    box_size(MINMAX_PATH, EXTRA)
    Vina_configuration_files(CENTER_PATH, MODES, CPU, EXH)
    logging.info("Finished Dockatron10XD Part 2")
    print("Finished Dockatron10XD Part 2 \n")
   

